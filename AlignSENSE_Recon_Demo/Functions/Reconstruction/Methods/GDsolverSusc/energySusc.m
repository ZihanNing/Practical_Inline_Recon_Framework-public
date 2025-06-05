function [EnFi, EnRe, EnReNames] = energySusc(y, x, E, EH, W, sep,deb )

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
if nargin < 7 || isempty(deb); deb=0;end

%%%DATA FIDELITY 
EnFi = 0.5 *  multDimSum(computeEnergy(y,x,E,[],EH))/numel(y) ; %Mean here goes accompanied by mulDimMea in gradient

%%%REGULARISATION
EnRe = [];
EnReNames = {};

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_sparse') && ~isequal(E.Ds.Reg.lambda_sparse,0) %L1 norm 
    temp = cat(4, FiniteDiff(E.Ds.chi,1), FiniteDiff(E.Ds.chi,2), FiniteDiff(E.Ds.chi,2)); 
    EnRe = cat(1,EnRe , 0.5 * E.Ds.Reg.lambda_sparse * multDimSum( abs(temp)) ) /numel(E.Ds.D);
    EnReNames = cat(1,EnReNames, 'sparse');
end

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_smooth') && ~isequal(E.Ds.Reg.lambda_smooth,0) %Smoothness 
    EnRe = cat(1,EnRe, 0.5 * E.Ds.Reg.lambda_smooth * ( normm( sqrt(W) .* FiniteDiff(E.Ds.chi,1))/numel(E.Ds.chi) + ...  % mean here goes accompanied by mulDimMea in gradient
                                                        normm( sqrt(W) .* FiniteDiff(E.Ds.chi,2))/numel(E.Ds.chi) + ...
                                                        normm( sqrt(W) .* FiniteDiff(E.Ds.chi,3))/numel(E.Ds.chi) ));
    EnReNames = cat(1,EnReNames, 'smooth');
end

if isfield(E.Ds,'Reg') && isfield(E.Ds.Reg,'lambda_tissueConfidence') && ~isequal(E.Ds.Reg.lambda_tissueConfidence,0)%Tissue confidence 
    EnRe = cat(1,EnRe, 0.5 * E.Ds.Reg.lambda_tissueConfidence * normm( sqrt(W).* (E.Ds.chi-E.Ds.Reg.tissueSusc) ) /numel(E.Ds.D));
    EnReNames = cat(1,EnReNames, 'ridge');
    %fprintf('energyTaylor:: Reg Ridge: %d\n', dynInd(EnRe,size(EnRe,1), 1) );  
end

%%%GATHER
EnFi = gather(EnFi);
EnRe = gather(EnRe);

if sep==0 && ~isempty(EnRe); EnFi = EnFi + multDimSum(EnRe,1); end

end

