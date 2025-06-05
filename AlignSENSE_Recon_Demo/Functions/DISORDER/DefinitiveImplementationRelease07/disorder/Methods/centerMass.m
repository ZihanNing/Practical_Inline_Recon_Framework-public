function c=centerMass(x)

% CENTERMASS computes the center of mass of a N-D array
%   C=CENTERMASS(X) 
%   * X is the array
%   * C is the center of mass
%

gpu=isa(x,'gpuArray');
N=size(x);
ND=numDims(x);N=N(1:ND);
rGrid=generateGrid(N,gpu,N,zeros(1,ND));

mass=sum(x(:));
c=zeros(1,ND);
for n=1:ND;c(n)=multDimSum(bsxfun(@times,rGrid{n},x),1:ND)/mass;end
