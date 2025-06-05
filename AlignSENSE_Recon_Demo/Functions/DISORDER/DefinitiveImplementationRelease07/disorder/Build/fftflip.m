function x=fftflip(x,dim)

%FFTFLIP   Flips a spectrum on the basis of Hermitian symmetry
%   X=FFTFLIP(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to flip, defaults to all
%   ** X is the flipped array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

x=conj(x);
for n=1:length(dim)    
    v=dim(n);
    x=flip(x,v);
    if mod(size(x,v),2)==0;x=circshift(x,1,v);end
end
