function [x,y]=multDimMin(x,dim)

%MULTDIMMIN   Takes the minimum of the elements of a multidimensional 
%array along a set of dimensions
%   X=MULTDIMMIN(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the min of the elements 
%   of the array, defaults to all
%   ** X is the contracted array
%   ** Y is a set of indexes with the minimum elements
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

if nargout==1
    for n=1:length(dim);x=min(x,[],dim(n));end
else%We assume new matlab versions
    N=size(x);
    [x,y]=min(x,[],dim,'linear');
    y=y(:);
    y=ind2subV(N,y);
end
