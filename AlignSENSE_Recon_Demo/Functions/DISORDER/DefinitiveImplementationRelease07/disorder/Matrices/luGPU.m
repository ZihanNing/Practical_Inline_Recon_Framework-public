function [X,U]=luGPU(X)

%LUGPU   Performs a LU decomposition of a set of matrices on the GPU. If one output 
%is requested, the upper, diagonal and lower triangular matrices are overwritten in 
%X. If two outputs, the method returns upper triangular and unit diagonal lower 
%triangular results. Pivoting is not implemented
%   [X,U]=LUGPU(X)
%   * X is the set of input matrices
%   * X are unit diagonal lower triangular matrices (two outputs) or the full 
%   factorization (one output)
%   * U are upper diagonal matrices
%

M=size(X);M(end+1:3)=1;
N=M(1);

X=permute(X,[3 1 2]);
for n=1:N-1
    nk=n+1:N;
    X(:,nk,n)=bsxfun(@rdivide,X(:,nk,n),X(:,n,n));             
    X(:,nk,nk)=X(:,nk,nk)-bsxfun(@times,X(:,nk,n),X(:,n,nk));
end
X=permute(X,[2 3 1]);

if nargout==2
    U=triuGPU(X);
    X=X-diagm(diagm(X)-1);
    X=trilGPU(X);
end
