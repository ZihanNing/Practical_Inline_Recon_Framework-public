function X=triBackGPU(U,X)

%TRIBACKGPU   Solves a set of triangular systems of equations by 
%backsubstitution. The code is based on
%https://uk.mathworks.com/matlabcentral/fileexchange/37459-matrix-inverse-using-lu-factorization, 
%which gives the following documentation
% Solves C = U \ B;
%          |1|         |2 2 1|
% With b = |2| and U = |0 1 4|
%          |3|         |0 0 3|
%   X=TRIBACKGPU(U,{X})
%   * U is an upper triangular matrix
%   * {X} are the coefficients to solve for (if empty, the identity matrix is used)
%   * X is the solved system of equations
%

M=size(U);M(end+1:3)=1;
N=M(1);
Ud=permute(diagm(U),[2 1 3]);
if nargin<2 || isempty(X);X=eye(N,'like',U);X=repmat(X,[1 1 M(3)]);end
X=bsxfun(@rdivide,X,Ud);
U=bsxfun(@rdivide,U,Ud);
X=permute(X,[3 2 1]);
U=permute(U,[3 1 2]);
for n=N-1:-1:1
    vn=n+1:N;
    X(:,:,n)=X(:,:,n)-sum(bsxfun(@times,U(:,n,vn),X(:,:,vn)),3);
end
X=permute(X,[3 2 1]);