function x=aplGPU(A,x,m)

%APLGPU   Applies the matrix A to the array x along dimension m of x
%   X=APLGPU(A,X,M)
%   * A is a matrix to apply over A.
%   * X is the array on which to apply A
%   * M is the direction along which to apply A over X
%   ** X is the transformed array
%

if isempty(x) || isempty(A);return;end

ND=numDims(x);
N=size(x);N(end+1:max(ND+1,m+1))=1;
if m~=1;x=reshape(x,[prod(N(1:m-1)) N(m) prod(N(m+1:ND))]);else x=x(:,:);end

if m==1;x=A*x;
%elseif m<ND;x=pagemtimes(x,A.');
elseif m<ND;x=matfun(@mtimes,x,A.');
else x=x*A.';
end
N(m)=size(A,1);
x=reshape(x,N);
