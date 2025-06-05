
function [EnFi, EnRe, EnReNames] = energyLinear(y, x, E, EH, W, sep )

%ENERGYLINEAR   Computes the energy of the solution of a least squares 
%problem with different regularisations for the linaer model of pose-dependent B0 variations. 
%   [EN,W]=ENERGYLINEAR(Y,X,E,EH,{W},{SEP})
%   * Y is the measured data
%   * X is the reconstructed data
%   * E is the encoding structure
%   * EH is the decoding structure (used for weighted least squares)
%   * {W} is the weight to use for certain regularisers
%   * {SEP} whether to separate the data fidelity loss and the regularisatoin loss
%
if nargin < 5 || isempty(W); W = ones(size(x), 'like',real(x));end
if nargin < 6 || isempty(sep); sep=(nargout>1);end

%%%DATA FIDELITY 
EnFi = 0.5 *  multDimSum(computeEnergy(y,x,E,[],EH))/numel(y) ; %Mean here goes accompanied by mulDimMea in gradient

%%%REGULARISATION
EnRe = [];
EnReNames = {};

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_sparse') && ~isequal(E.Dl.Reg.lambda_sparse,0) %L1 norm 
    EnRe = cat(1,EnRe , 0.5 * E.Dl.Reg.lambda_sparse * multDimSum( abs(bsxfun(@times,ones(size(W)),E.Dl.D)))/numel(E.Dl.D) ) ;
    EnReNames = cat(1,EnReNames, 'sparse');
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_smooth') && ~isequal(E.Dl.Reg.lambda_smooth,0) %Smoothness 
    EnRe = cat(1,EnRe, 0.5 * E.Dl.Reg.lambda_smooth * ( normm( sqrt(W) .* FiniteDiff(E.Dl.D,1))/numel(E.Dl.D) + ...  % mean here goes accompanied by mulDimMea in gradient
                                                        normm( sqrt(W) .* FiniteDiff(E.Dl.D,2))/numel(E.Dl.D) + ...
                                                        normm( sqrt(W) .* FiniteDiff(E.Dl.D,3))/numel(E.Dl.D) ));
    EnReNames = cat(1,EnReNames, 'smooth');
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_ridge') && ~isequal(E.Dl.Reg.lambda_ridge,0)%Ridge 
    EnRe = cat(1,EnRe, 0.5 * E.Dl.Reg.lambda_ridge * normm(E.Dl.D)/numel(E.Dl.D));
    EnReNames = cat(1,EnReNames, 'ridge');
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_susc') && ~isequal(E.Dl.Reg.lambda_susc,0)%Implicit susceptibility model
    N = size(E.Dl.D);
    K = linearTerms(N,pi,6);
    
    temp = dynInd(K,2,6) .* fftn(fftshift(dynInd(E.Dl.D,1,6))) + dynInd(K,1,6) .* fftn(fftshift(dynInd(E.Dl.D,2,6))); %TO DO: make this flexible that different E.Dl.d work with the right Kx and Ky
    temp = sqrt(abs(W)) .* ifftshift(ifftn(temp));  %TO DO: use fftGPU

    EnRe = cat(1,EnRe, 0.5 * E.Dl.Reg.lambda_susc * normm(temp)/numel(E.Dl.D) );
    EnReNames = cat(1,EnReNames, 'susc');
end

%%%GATHER
EnFi = gather(EnFi);
EnRe = gather(EnRe);

if sep==0 && ~isempty(EnRe); EnFi = EnFi + multDimSum(EnRe,1); end

end

