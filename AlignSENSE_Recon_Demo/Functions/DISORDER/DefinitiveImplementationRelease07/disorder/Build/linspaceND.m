function x=linspaceND(d1,d2,N)

%LINSPACEND Linearly spaced multidimensional array quite freely based on the code at 
%https://uk.mathworks.com/matlabcentral/fileexchange/22824-linearly-spaced-multidimensional-matrix-without-loop
%Note that when scalars are inputted, return array is a column vector instead of the row vector of linspace
%   X=LINSPACEND(D1,D2,{N})
%   * D1 are the lower limits of the array
%   * D2 are the upper limits of the array
%   * {N} is the number of points to be generated. It defaults to 100 for compatibility with linspace function
%   * X is the generated array
%

if ~exist('N','var') || isempty(N);N=100;end;

gpu=isa(d1,'gpuArray');
N=double(floor(N));

ND=numDims(d1);
M=size(d1);
assert(ND==numDims(d2),'d1 and d2 must have the same number of dimensions but their dimensionality is %d %d',ND,numDims(d2));
assert(all(M==size(d2)),'d1 and d2 must have the same size but their sizes are%s and%s',sprintf(' %d',M),sprintf(' %d',size(d2)));

x=single((0:N-1)/(N-1));
if gpu;x=gpuArray(x);end
perm=1:max(ND+1,2);perm(2)=ND+1;perm(ND+1)=2;
x=permute(x,perm);
x=bsxfun(@plus,d1,bsxfun(@times,x,d2-d1));
