

function [E, x, En] = GDsolver(y,E,EH,A,C,R,x,nIt,tolType,tol, deb)

%GDSOLVER   Performs a GD-based pseudoinverse reconstruction for pose-dependent B0 fields
%   [E]=CGSOLVER(Y,E,EH,{A},{C},{R},{X},{NIT},{TOLTYPE},{TOL},{DEB})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {A} is a preconditioner structure
%   * {C} is a constrain structure
%   * {R} is a regularizer structure
%   * {X} is an initial condition for the array to be reconstructed
%   * {NIT} is the maximum number of iterations
%   * {TOL} is the tolerance
%   * {DEB} indicates whether to print information about convergence

if nargin<4;A=[];end
if nargin<5;C=[];end
if nargin<6;R=[];end
if nargin<8 || isempty(nIt);nIt=300;end
if nargin<9 || isempty(tolType);tolType='Energy';end % Residuals used for now 
if nargin<10 || isempty(tol);tol=1e-5;end
if nargin<11 || isempty(deb);deb=0;end

%%% INITIALISE
gpu = isa(x,'gpuArray');
Ncr = size(E.Db.cr); Ncr(end+1:6) =1;
NX = size(x);
NStates = E.NMs; 
NSegments = E.NSe; %inlcudes the outer segment (shutter)

if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'BlSzCoef');BlSzCoef = E.Db.Optim.BlSzCoef; else BlSzCoef=5 ;end % Block size for coefficient gradient calculation
if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'BlSzStates');BlSzStates = E.Db.Optim.BlSzStates; else BlSzStates=E.bS(1) ;end% Block size for coefficient gradient calculation
if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'meanC');meanC = E.Db.Optim.meanC; else meanC=0 ;end
if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'deTaylor');deTaylor = E.Db.Optim.deTaylor(1); else deTaylor=0 ;end%Note E.Db.Optim.deTaylor(2) contains information about outer iterations (handled in solveXTB_ext)
if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'Mbg');Mbg = E.Db.Optim.Mbg; else Mbg=1:NX(3);end
if isfield(E.Db,'Optim') && isfield(E.Db.Optim, 'alphaList');alphaList=E.Dl.Optim.alphaList;else alphaList=10.^(-1* [-3:1.3:8]); range = alphaList(1)/alphaList(end); end
if deTaylor<2; alphaBasis = ones([NStates 1]); else;alphaBasis = ones([1 1]);end % To refine the step size for certain motion states
x = constrain(x,C);% Constrain should already be applied in CGsolver

E = transformationConversion (E);

%%% SET ASIDE VARIABLES
oS = E.oS; bS = E.bS; dS = E.dS;
if isfield(EH,'Mbe') && ~isempty(EH.Mbe);Mbe = EH.Mbe;end %Energy computation coming from whole FOV
EH.Mbe = zeros([1 1 NX(3)],'like',real(x));EH.Mbe(Mbg)=1; %Energy computation in accordance with gradient support 
if isfield(EH,'We') && ~isempty(EH.We); We = EH.We; EH.We(:)=1;end %No shot weighting for coefficient estimation

%ENERGY PER STATES
Re = computeEnergy(y,x,E,[],EH,[],[],[],[],1)/numel(y);
ReSeg =  res_segments(Re, E.NMs, E.nEc);%ReSeg = accumarray(E.mSt, Re);%E.mSt not changed in compressMotion! might be an issue with segments/states offset (shutter)
if deTaylor>=2; ReSeg = multDimSum(ReSeg);end
En = ReSeg;

for n=1:max(nIt)    
    cr = E.Db.cr; % Set coeffients aside so E does not have to be duplicated
    res = y-encode(x,E);
    
    %%%GRADIENT INITIALISATION
    Ng = Ncr; Ng(5) = NSegments;
    grad =  zeros(Ng, 'like', real(E.Db.cr)); if gpu; grad = gpuArray(grad);end % Make sure g is not complex or it will raise errors
    
    for a=1:BlSzStates:NStates %Chunk states otwerwise is has size [Nx Nstates] and this runs out of memmory depending on NX and NStates
        E.oS(1) = a; E.bS(1) = BlSzStates; E.dS(1) = min(a+BlSzStates-1,NStates);
        vA=a:min(a+E.bS(1)-1,NStates);
        
        %%%DECODE (without summing step)
        notsum = 1; % if 0, decode step will sum over motion states
        [~,EHres] =  decode(res,EH,E,notsum); % motion states in the 5th dimension
        EHres = dynInd(EHres, Mbg,3);%take out part interested in
        
        for ii = 1:BlSzCoef:Ncr(6) %number of basis coefficients 
            vI = ii:min(ii+BlSzCoef-1 , Ncr(6));
            part1 = + 1i*2*pi*E.Db.TE * bsxfun( @times , conj(x(:)), E.Db.B(:,vI)); % + since field convention different than MRes project
            part1 = reshape(part1, [NX 1 1 length(vI)]);
            part1 = dynInd(part1, Mbg,3);
            
            dim_sum = 1:4; % sum over coils and space - not over motion states!
            temp = multDimSum( bsxfun( @times, part1, EHres), dim_sum )/numel(y); % has all the motion states in the 5th dimension
            grad = dynInd( grad, {vA,vI}, 5:6, real(temp) ) ;
        end
    end
    if deTaylor==2 
        grad = multDimSum( bsxfun(@times, grad, permute(E.Dl.f(E.Dl.Tr),[1:5 7 6])),5);
    elseif deTaylor==3
        grad = multDimSum(grad, 5);
    end
    E.oS(1)=oS(1); E.bS(1)=bS(1); E.dS(1)=dS(1); %store original oS, bS and dS
    
    
    %%%LINE SEARCH
    flag_all= zeros(size(alphaBasis));
    StepSize = zeros(size(alphaBasis));
    
    for alpha = alphaList
        %%% Possible update
        if deTaylor==2
            update = - alpha* bsxfun(@times, permute( alphaBasis,[2:5 1]), dynInd( grad ,':', 5)); %should only consist of 1 element in 5th dimension
            E.Db.cr =  cr + multDimSum(bsxfun(@times, permute(E.Dl.f(E.Dl.Tr),[1:5 7 6]),update),7);
        elseif deTaylor==3
            update = - alpha* bsxfun(@times, permute( alphaBasis,[2:5 1]), dynInd( grad , ':', 5)); %exlcude the shutter
            E.Db.cr =  cr + update ; 
        else
            update = - alpha* bsxfun(@times, permute( alphaBasis,[2:5 1]), dynInd( grad , 1:NStates, 5)); %exlcude the shutter
            E.Db.cr =  cr + update ; 
        end
        ReOption = computeEnergy(y,x,E,[],EH,[],[],[],[],1)/numel(y);
        ReSegOption = res_segments(ReOption, E.NMs, E.nEc);
        if deTaylor>=2; ReSegOption = multDimSum(ReSegOption);end
        
        %%% Flag energy reduction
        flag = (ReSegOption < ReSeg); % those segments have lower energy for that alpha
        ReSeg(flag) = ReSegOption(flag); % allow smaller step sizes after flag to reduce residual
        flag_all(flag) = 1; %if a step size if found
        
        %%% Store step size per state
        StepSize (flag==1) = alpha * alphaBasis(flag==1) ; % allow some states to have bigger step size
                
        if all(flag_all); break;end % step size for every state found
    end
    alphaBasis(~flag_all) =  alphaBasis(~flag_all) / range ; %refine step size for states where none found
    
    %%%UPDATE COEFFICIENTS
    if deTaylor==2
        update = - bsxfun(@times, permute(StepSize, [2:5 1]), dynInd(grad, ':' , 5)); %should only consist of 1 element in 5th dimension
        E.Db.cr =  cr + multDimSum(bsxfun(@times, permute(E.Dl.f(E.Dl.Tr),[1:5 7 6]),update),7);
    elseif deTaylor==3
        update = - bsxfun(@times, permute( StepSize,[2:5 1]), dynInd( grad , ':', 5)); 
        E.Db.cr =  cr + update ; 
    else
        update = - bsxfun(@times, permute(StepSize, [2:5 1]), dynInd(grad, 1:NStates , 5)); %exlcude the shutter
        E.Db.cr =  cr + update ; 
    end
    
    if ~all(flag_all) && deb>=2; fprintf('   GD iteration %d: step size zero for %d/%d\n', n, multDimSum(flag_all~=1), numel(flag));end
    if all(~flag_all);break;end%don't return since still have to substract mean field
    %need to return after break since update has to be zero
    if any(isnan(update),'all'); error('GD iteration %d: Update contains NaN\n', n); end
    
    %%%UPDATE RESIDUALS
    Re = computeEnergy(y,x,E,[],EH,[],[],[],[],1)/numel(y);
    ReSeg = res_segments(Re, E.NMs, E.nEc);
    if deTaylor>=2; ReSeg = multDimSum(ReSeg);end
    En = cat(2,En, ReSeg); %keep track of energy per state
    
figure(207)
plot(En')
pause(0.01)

    %%%CHECK CONVERGENCE
    if all( ReSeg<tol ) ;break;end % reached convergence 
end

if n==max(nIt) && deb==2;fprintf('   GD solver terminated without reaching convergence\n');end

%%% SUBSTRACT MEAN FIELD ACROSS MOTION STATES
if meanC
    fprintf('Taking out mean coefficients.\n');
    c_mean = multDimMea(E.Db.cr, 5 );
    E.Db.cr = bsxfun(@minus, E.Db.cr , c_mean);
    x = bsxfun(@times,x,dephaseBasis( E.Db.B, c_mean, size(x),E.Db.TE));
end

%%% SUBSTRACT LINEAR VARIATION W.R.T. PITCH AND ROLL
if deTaylor>=1 && deTaylor <3 && isfield(E,'Dl') && any(E.Dl.f(E.Dl.Tr)>0, 'all')
    fprintf('Taking out Taylor terms.\n');
                    
    C = permute (E.Db.cr, [5:6 1:4]);
    T = permute(E.Dl.f(E.Dl.Tr),[5:6 1:4]); %in radians
    Z = T\C; %LS of C = T*Z
    %%% Substract linear variation from c
    E.Db.cr = ipermute( C - T*Z, [5:6 1:4]);
    %%% Add to the linear model
    Dtemp = zeros(prod(NX) ,2, 'like', E.Dl.D); 
    BlSzC = 5; %For coefficients
    for jj= 1:BlSzC:Ncr(6)%matrix multiplication might become big - possible chunck into blocks
        vI=jj:min(jj+BlSzC-1,Ncr(6));
        Dtemp = Dtemp + dynInd(E.Db.B, vI,2) * dynInd(Z.',vI,1 );
    end
    E.Dl.D = E.Dl.D + reshape(Dtemp, [NX 1 1 2]);
        
end
if deTaylor==2;E.Db.cr = 0*E.Db.cr;fprintf('Setting coefficients to zero.\n');end%Remove residuals and only estimate the basis functions of the linear model

%%% RESTORE PARAMETERS
if exist('Mbe','var'); EH.Mbe = Mbe;end%currently not needed as EH not returned
if exist('We','var');  EH.We = We;end

end


