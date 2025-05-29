function [x,y]=multDimMax(x,dim)

%MULTDIMMAX   Takes the maximum of the elements of a multidimensional 
%array along a set of dimensions
%   X=MULTDIMMAX(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the max of the elements 
%   of the array, defaults to all
%   ** X is the contracted array
%   ** Y is a set of indexes with the maximum elements
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

if nargout==1
    for n=1:length(dim);x=max(x,[],dim(n));end
else%We assume new matlab versions
    N=size(x);
    [x,y]=max(x,[],dim,'linear');
    y=y(:);
    y=ind2subV(N,y);
end
