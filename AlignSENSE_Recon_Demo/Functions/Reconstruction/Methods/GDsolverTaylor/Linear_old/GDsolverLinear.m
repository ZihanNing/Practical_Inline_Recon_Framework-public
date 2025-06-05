
function [E, En] = GDsolverLinear(y,E,EH,A,C,R,x,nIt, tol, deb)

%GDSOLVER   Performs a GD-based pseudoinverse reconstruction for linear model of the pose-dependent B0 fields
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
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'alphaList');   alphaList=E.Dl.Optim.alphaList;else alphaList=10.^(-1* [-0:0.7:6]); end%line search step size
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'FISTA');       FISTA=E.Dl.Optim.FISTA;else FISTA=1; end
if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg, 'lambda_sparse');   lambda_sparse=E.Dl.Reg.lambda_sparse;else lambda_sparse=0; end
th = 1; temp = E.Dl.D;
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'DupLim');  DupLim=E.Dl.Optim.DupLim;else DupLim=inf; end%constrain update TODO include in constrainLinear.m
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'DLim');    DLim=E.Dl.Optim.DLim;else DLim=inf; end%constrain update TODO include in constrainLinear.m
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'RMSProp'); RMSProp=E.Dl.Optim.RMSProp;else RMSProp=0; end
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'beta');    beta=E.Dl.Optim.beta;else beta=0.9; end
if isfield(E.Dl,'Optim') && isfield(E.Dl.Optim, 'Mbg');     Mbg = E.Dl.Optim.Mbg; else Mbg=1:NX(3);end
assert(~(FISTA==1 && RMSProp==1), 'FISTA and RMSProp should not be activated together.\n');
%visSegment(dynInd(x, Mbg,3), [],0)

%Energy
if isfield(EH,'Mbe');Mbe = EH.Mbe;end %In contrast to motion estimation, use whole readout dimension - NOT USED FOR NOW
EH.Mbe = zeros([1 1 NX(3)],'like',real(x));EH.Mbe(Mbg)=1;
if isfield(EH,'We')&& ~isempty(EH.We)
    We = EH.We; %set aside
    EH.Mbe= bsxfun(@times, We, EH.Mbe ); %Make sure these weights are also used in computeEnergy.m
end
if isfield(E.Dl.Reg ,'W');  W = E.Dl.Reg.W; else W = ones(size(x), 'like', real(x));end % Weighted regularisation
EnPrevious = energyLinear(y, x, E, EH,  W);
En = EnPrevious; % keep track of energy

%%% ITERATE
for n=1:max(nIt)    
    D = E.Dl.D; % Set original value aside so no need to duplicate E

    %%%GRADIENT
    grad = gradLinear(y,x,E,EH,W);
    [~,grad] = constrainLinear(E,grad); %Constrain gradient
    if isfield(E.Dl.C , 'Ma'); grad = bsxfun(@times, E.Dl.C.Ma, grad);end
        
    %plot_(abs(x), [],[],dynInd(grad,1,6),dynInd(grad,2,6),[-10 10],[],1,[],[],789);        
    
    %%%GRADIENT SCALING
    if RMSProp
        if ~exist('v','var'); v=zeros(size(E.Dl.D));end
        v = beta*v  + (1-beta)*grad.^2;
        grad = grad ./ sqrt(v+1e-3);%grad(grad>0.1) = grad(grad>0.1) ./ sqrt(v(grad>0.1)+1e-3);
    end
    
    %%% (F)ISTA LINE SEARCH
    temp_old = temp; %FISTA parameters
    if FISTA;th_old = th; else th_old=1; end
    th = (1+sqrt(1+4*th_old^2))/2;
    
    for alpha = alphaList
        %%% (F)ISTA Possible update 
        update = - alpha* grad; 
        idx = abs(update)>DupLim; update(idx)= DupLim .* sign(update(idx));% Constrain update
        temp = SoftThresh( bsxfun(@times,ones(size(W)),D + update) , alpha * lambda_sparse/ numel(E.Dl.D) );
        %temp = bsxfun(@times, conj(W),temp);
        E.Dl.D = temp + (th_old - 1)/th * (temp - temp_old) ;
        idx = abs(E.Dl.D)>DLim;  E.Dl.D(idx)= DLim .* sign( E.Dl.D(idx));% Constrain new estimate
        
        %%% Check possible energy reduction 
        EnNew = energyLinear(y,x,E,EH,W);         
        if  EnNew < EnPrevious
            EnPrevious = EnNew;
            break;
        end 
        
        if alpha == alphaList(end) 
            if deb >=2; fprintf('   GD iteration %d: step size zero\n   GDsolver aborted\n', n);end
            E.Dl.D = D; % Restore with original value
            return %don't rescale step size and leave GDsolverLinear
        end
    end
        
    %%%CHECK CONVERGENCE
    En = cat(2, En, EnNew); 
    if EnNew<tol;break;end % reached convergence   
    if n==max(nIt) && deb==2;fprintf('   GD solver terminated without reaching convergence\n');end

    if deb
        figure(307); 
        plot(1:length(En), En)
    end

end

%%%CONSTRAIN
[E,E.Dl.D] = constrainLinear(E); 
    
%%%SET BACK PARAMETERS
if exist('Mbe', 'var');EH.Mbe=Mbe;end
if exist('We', 'var');EH.We=We;end


end