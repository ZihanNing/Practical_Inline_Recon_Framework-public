
function [E, En] = GDsolverTaylor(y,E,EH,A,C,R,x,nIt,tol, deb)

%GDSOLVERTAYLOR  Performs a GD-based pseudoinverse reconstruction for the Taylor model describing the pose-dependent B0 as proposed in 
% "Data-driven motion-corrected brain MRI incorporating pose
% dependent B0 fields", Y Brackenier, L Cordero-Grande, R Tomi-Tricot, T
% Wilkinson, P Bridgen, A Price, S J Malik, E De Vita and J V Hajnal, 2022..
%   [E, En]=GDSOLVERTAYLOR(Y,E,EH,{A},{C},{R},X,{NIT},{TOL},{DEB})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {A} is a preconditioner structure - not implemented yet
%   * {C} is a constrain structure - not implemented yet
%   * {R} is a regularizer structure - not implemented yet
%   * X is an current estimate of the array to be reconstructed
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
NSegments  = E.NSe; %Inlcudes the outer segment (shutter)

%(F)ISTA parameters
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'alphaList');   alphaList=E.Dl.Optim.alphaList;else alphaList=10.^(-1* [-0:0.7:6]); end%Line search step size
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'FISTA');       FISTA=E.Dl.Optim.FISTA;else FISTA=1; end%Flag to use FISTA
if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg, 'lambda_sparse');   lambda_sparse=E.Dl.Reg.lambda_sparse;else lambda_sparse=0; end%Hyper-parmeter for L1 loss
th = 1; temp = E.Dl.D;
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'DupLim');  DupLim=E.Dl.C.DupLim;else DupLim=inf; end%Maximum update 
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'DLim');    DLim=E.Dl.C.DLim;else DLim=inf; end%Minimum update 
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'RMSProp'); RMSProp=E.Dl.Optim.RMSProp;else RMSProp=0; end%Flag to use RMSEProp
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'beta');    beta=E.Dl.Optim.beta;else beta=0.9; end%Hyper-parmeter regarding momentum gradient calculation
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'Mbg');     Mbg = E.Dl.Optim.Mbg; else Mbg=1:NX(3);end%Are in readout to consider for gradient calculation
assert(~(FISTA==1 && RMSProp==1), 'GDsolverTaylor:: FISTA and RMSProp should not be activated together.\n');

%Energy
if isfield(EH,'Mbe');Mbe = EH.Mbe;end 
EH.Mbe = zeros([1 1 NX(3)],'like',real(x));EH.Mbe(Mbg)=1;
if isfield(EH,'We')&& ~isempty(EH.We)
    We = EH.We; %set aside
    EH.Mbe= bsxfun(@times, We, EH.Mbe ); %Make sure these weights are also used in computeEnergy.m
end
if isfield(E.Dl.Reg ,'W');  W = E.Dl.Reg.W; else W = ones(size(x), 'like', real(x));end % Weighted regularisation
EnPrevious = energyTaylor(y, x, E, EH,  W);
En = EnPrevious; %Keep track of energy

 updateSeparate=0;
 ff = ones([ones(1,5),size(E.Dl.D,6)]);
 if updateSeparate; ff = 0*repmat(ff,[ones(1,6),size(E.Dl.D,6)]);for q=1:size(ff,7);ff = dynInd(ff,[{q,q}],[6,7],1);end; end
    
%%% ITERATE
for n=1:max(nIt)    
    D = E.Dl.D; %Set original value aside so no need to duplicate E
    
    %%%GRADIENT
    grad = gradTaylor(y,x,E,EH,W);
    if isfield(E.Dl.C , 'Ma'); grad = bsxfun(@times, E.Dl.C.Ma, grad);end
    
    %%%GRADIENT SCALING
    if RMSProp
        if ~exist('v','var'); v=zeros(size(E.Dl.D));end
        v = beta*v  + (1-beta)*grad.^2;
        grad = grad ./ sqrt(v+1e-3);
    end
    
    %%% (F)ISTA LINE SEARCH
    temp_old = temp; %FISTA parameters
    if FISTA;th_old = th; else th_old=1; end
    th = (1+sqrt(1+4*th_old^2))/2;
     
    for q =1:size(ff,7)
       D = E.Dl.D;
       for alpha = alphaList
            %%% (F)ISTA Possible update 
            update = - alpha* bsxfun(@times, dynInd(ff,q,7),grad); %only selecting one Taylor term (all others to zero)
            idx = abs(update)>DupLim; update(idx)= DupLim .* sign(update(idx));% Constrain update
            temp = SoftThresh( bsxfun(@times,ones(size(W)),D + update) , alpha * lambda_sparse );%/ numel(E.Dl.D)
            %temp = bsxfun(@times, conj(W),temp);
    
            E.Dl.D = temp + (th_old - 1)/th * (temp - temp_old) ;
            idx = abs(E.Dl.D)>DLim;  E.Dl.D(idx)= DLim .* sign( E.Dl.D(idx));% Constrain new estimate

            %%% Check possible energy reduction 
            EnNew = energyTaylor(y,x,E,EH,W,[],0);         
            if  EnNew < EnPrevious
                EnPrevious = EnNew;
                break;
            end 

            if alpha == alphaList(end) 
                if deb >=2; fprintf('   GD iteration %d: step size zero\n   GDsolver aborted\n', n);end
                E.Dl.D = D; % Restore with original value
                return %don't rescale step size and leave GDsolverTaylor
            end
        end
     end 
    
    %%%CHECK CONVERGENCE
    En = cat(2, En, EnNew); 
    if EnNew<tol;break;end %Reached convergence   
    if n==max(nIt) && deb==2;fprintf('   GD solver terminated without reaching convergence\n');end

    if deb
        h=figure(307); 
        set(h,'color','w');
        plot(1:length(En), En);
        xlabel('Number of iterations [#]');
        ylabel('Residuals [a.u.]');
        title('Convergenve plot for GDsolverTaylor');
    end

end

%%%CONSTRAIN
[E,E.Dl.D] = constrainTaylor(E);

%%%SET BACK PARAMETERS
if exist('Mbe', 'var');EH.Mbe=Mbe;end
if exist('We', 'var');EH.We=We;end

end