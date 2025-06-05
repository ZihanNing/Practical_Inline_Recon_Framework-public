function X=zrGPU(X,r,di)

%ZRGPU   Zeroes those rows/columns in the matrices X above a given rank R
%(in general those elements of an array along a given dimension which are
%above the rank given along another dimension or set of dimensions)
%   X=ZRGPU(X,R)
%   * X are the input matrices
%   * R is the number of rows of interest
%   * DI is the direction for zeroing, defaults to 1
%   * X are the output matrices
%

if nargin<3 || isempty(di);di=1;end

gpu=isa(X,'gpuArray');
N=size(X);N(end+1:di)=1;
NR=size(r);NR(end+1:di)=1;

assert(NR(di)==1,'The indexes for zeroing should not be defined along the dimension to zero, so the dimensionality should be zero, but it is: %d',NR(di));
inddi=1:N(di);
if gpu;inddi=gpuArray(inddi);end
perm=1:max(di,2);
if di~=2;perm([2 di])=[di 2];end
X=bsxfun(@times,X,bsxfun(@le,permute(inddi,perm),r));
