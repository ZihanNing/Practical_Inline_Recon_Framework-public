function y=multDimSize(x,dim)

%MULTDIMSIZE
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end
N =size(x);
%x = N(dim);
y=[];
for n=1:length(dim);y=cat(2,y,N(dim(n)));end
