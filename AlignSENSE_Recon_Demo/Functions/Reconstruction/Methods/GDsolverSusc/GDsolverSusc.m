
function [E, En] = GDsolverSusc(y,E,EH,A,C,R,x,nIt, tol, deb)

%GDSOLVERSUSC   Performs a GD-based pseudoinverse reconstruction for linear
%model of the susceptibility model causing pose-dependent sptially varying phase
%   [E, En]=CGSOLVER(Y,E,EH,{X},{NIT},{TOL},{DEB})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {X} is an initial condition for the array to be reconstructed
%   * {NIT} is the maximum number of iterations
%   * {TOL} is the tolerance
%   * {DEB} indicates whether to print information about convergence

if nargin<5 || isempty(nIt);nIt=300;end
if nargin<6 || isempty(tol);tol=1e-5;end
if nargin<7 || isempty(deb);deb=0;end

%%% INITIALISE
%General parameters
gpu = isa(x,'gpuArray');
NX = size(x);
NStates  = E.NMs; 
NSegments  = E.NSe; %inlcudes the outer segment (shutter)

%(F)ISTA parameters
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'lineSearch');   lineSearch=E.Ds.Optim.lineSearch;else lineSearch=1; end%line search flag
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'alphaList');   alphaList=E.Ds.Optim.alphaList;else alphaList=10.^(-1* [-0:0.7:6]); end%line search step size
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'FISTA');       FISTA=E.Ds.Optim.FISTA;else FISTA=1; end
if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg, 'lambda_sparse');   lambda_sparse=E.Ds.Reg.lambda_sparse;else lambda_sparse=0; end
th = 1; temp = E.Ds.chi;
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'chiUpLim');  chiUpLim=E.Ds.Optim.chiUpLim;else chiUpLim=inf; end%constrain update TODO include in constrainSusc.m
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'chiLim');    chiLim=E.Ds.Optim.chiLim;else chiLim=[-inf inf]; end%constrain update TODO include in constrainSusc.m
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'RMSProp'); RMSProp=E.Ds.Optim.RMSProp;else RMSProp=0; end
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'beta');    beta=E.Ds.Optim.beta;else beta=0.9; end
if isfield(E.Ds,'Optim') && isfield(E.Ds.Optim, 'Mbg');     Mbg = E.Ds.Optim.Mbg; else Mbg=1:NX(3);end
assert(~(FISTA==1 && RMSProp==1), 'GDsolverSusc:: FISTA and RMSProp should not be activated together.\n');

%Energy
if isfield(EH,'Mbe');Mbe = EH.Mbe;end %In contrast to motion estimation, use whole readout dimension - NOT USED FOR NOW
EH.Mbe = zeros([1 1 NX(3)],'like',real(x));EH.Mbe(Mbg)=1;
if isfield(EH,'We')&& ~isempty(EH.We)
    We = EH.We; %set aside
    EH.Mbe= bsxfun(@times, We, EH.Mbe ); %Make sure these weights are also used in computeEnergy.m
end
if isfield(E.Ds.Reg ,'W');  W = E.Ds.Reg.W; else W = ones(size(x), 'like', real(x));end % Weighted regularisation
W = padArrayND(W,E.Ds.padDim, 1, 0,'both');
EnPrevious = energySusc(y, x, E, EH,  W);
En = EnPrevious; % keep track of energy
    
%%% ITERATE
for n=1:max(nIt)    
    chi = E.Ds.chi; % Set original value aside so no need to duplicate E
    
    %%%GRADIENT
    grad = gradSusc(y,x,E,EH,W);
    [~,grad] = constrainSusc(E,grad); %Constrain gradient
    if isfield(E.Ds,'C') && isfield(E.Ds.C,'Ma'); grad = bsxfun(@times, E.Ds.C.Ma, grad);end
    %plotND([],permute(grad,[1:3 6 4 5]),[],[],1,{[],2},[],[],[],[],E.Ds.Geom.APhiRecOrig*Tpermute(E.Ds.Geom.permuteHist{1})*Tpermute(E.Ds.Geom.permuteHist{2}));
    
    %%%GRADIENT SCALING
    if RMSProp
        if ~exist('v','var'); v=zeros(size(E.Ds.chi));end
        v = beta*v  + (1-beta)*grad.^2;
        grad = grad ./ sqrt(v+1e-3);%grad(grad>0.1) = grad(grad>0.1) ./ sqrt(v(grad>0.1)+1e-3);
    end
    
    %%% (F)ISTA LINE SEARCH
    temp_old = temp; %FISTA parameters
    if FISTA; th_old = th; else th_old=1; end
    th = (1+sqrt(1+4*th_old^2))/2;
     
   for alpha = alphaList
        %%% (F)ISTA Possible update 
        update = - alpha* grad; %only selecting one Taylor term (all others to zero)
        idx = abs(update)>chiUpLim; update(idx)= chiUpLim .* sign(update(idx));% Constrain update

        temp = chi+update;
        tempMean = temp; %multDimMea(temp);
        %plotND([], temp, [],[],0,[],[],{'temp'},[],[],89);
        if lambda_sparse>0
            temp = cat(4, FiniteDiff(temp,1), FiniteDiff(temp,2),FiniteDiff(temp,3));
            %plotND([], temp, [],[],0,[],[],{'f diff'},[],[],90);
            temp = SoftThresh( temp , alpha * lambda_sparse );%/ numel(temp);
            %plotND([], temp, [],[],0,[],[],{'f diff thresholded'},[],[],91);
            temp = tempMean + FiniteDiff(dynInd(temp,1,4),1,0)+FiniteDiff(dynInd(temp,2,4),2,0)+FiniteDiff(dynInd(temp,3,4),3,0);
        end
        %plotND([], temp, [],[],0,[],[],{'temp new'},[],[],92);
        E.Ds.chi = temp + (th_old - 1)/th * (temp - temp_old) ;
        idxLow = (E.Ds.chi)>chiLim(2);  E.Ds.chi(idxLow)= chiLim(2);% .* sign( E.Ds.chi(idx));% Constrain new estimate
        idxUp = (E.Ds.chi)<chiLim(1);  E.Ds.chi(idxUp)= chiLim(1);% 
        E=convertChiToTaylorMaps(E);%Update of chi must be converted to update D since this is called in encode
        
        %%% Check possible energy reduction 
        EnNew = energySusc(y,x,E,EH,W,[],0);         
        if  ~lineSearch || (EnNew < EnPrevious)
            EnPrevious = EnNew;
            break;
        end 

        if alpha == alphaList(end) 
            if deb >=2; fprintf('   GD iteration %d: step size zero\n   GDsolver aborted\n', n);end
            E.Ds.chi = chi; % Restore with original value
            E=convertChiToTaylorMaps(E);
            return %don't rescale step size and leave GDsolver
        end
    end
    
    %%%CHECK CONVERGENCE
    En = cat(2, En, EnNew); 
    if EnNew<tol;break;end % reached convergence   
    if n==max(nIt) && deb==2;fprintf('   GD solver terminated without reaching convergence\n');end

    if deb
        h=figure(307); 
        set(h,'color','w');
        plot(1:length(En), En);
        xlabel('Number of iterations [#]');
        ylabel('Residuals [a.u.]');
        title('Convergenve plot for GDsolverSusc');
    end

end

%%%CONSTRAIN
% %if isfield(E.Ds.C,'H'); H = E.Ds.C.H; E.Ds.C.H = [];end
% [E,E.Ds.chi] = constrainSusc(E);
% %if exist('H','var');E.Ds.C.H = H;end

%%%SET BACK PARAMETERS
if exist('Mbe', 'var');EH.Mbe=Mbe;end
if exist('We', 'var');EH.We=We;end

end