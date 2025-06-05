
function [x, P]=multDimRidge(x,dim,dimRidge,up,P)

%MULTDIMRIDGE   
%   X=MULTDIMRIDGE(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to performe the RIDGE
%   * {DIMRIDGE} is the dimension along which to detect the peak
%   * {UP} is the upsampling factor used for peak detection.
%   * {P} are the ridge parameters to calculate phase from. If empty, they are calculated.
%   ** X is the contracted array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end
if nargin<3 || isempty(dimRidge);dimRidge=1:numDims(x);end
if nargin<4 || isempty(up);up=1;end
if nargin >4 &&~isempty(P); di=1;else di=0;end

N = size(x);

%%% Sum dimensions not to estimate linear variation
% idx = ismember(dim, dimRidge);
% dimSum = dim; dimSum(idx)=[];
dimSum = dim(~ismember(dim, dimRidge));
if ~isempty(dimSum);x = multDimSum( x, dimSum);end

%%% Detect/apply phase ramp and offset
if di==0%Detect
    [x,P] = ridgeDetectionExt(x, dimRidge, up);
else%Apply
    [x,P] = ridgeDetectionExt(x, dimRidge, up, P);
end

%%% Repeat array to match original size
Nrep = ones(size(N)); Nrep(dimSum)=N(dimSum);
x = repmat(x, Nrep);
