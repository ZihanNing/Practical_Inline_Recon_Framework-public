function X=triForwGPU(L,X)

%TRIFORWKGPU   Solves a set of triangular systems of equations by 
%forwardsubstitution. The code is based on
%https://uk.mathworks.com/matlabcentral/fileexchange/37459-matrix-inverse-using-lu-factorization, 
%which gives the following documentation
% Solves C = L \ b
%          |1|         |1 0 0|
% With b = |2| and L = |2 1 0|
%          |3|         |3 4 1|
%   X=TRIBACKGPU(L,{X})
%   * L is a diagonally unitary lower triangular matrix
%   * {X} are the coefficients to solve for (if empty, the identity matrix is used)
%   * X is the solved system of equations
%

M=size(L);M(end+1:3)=1;
N=M(1);
if nargin<2 || isempty(X);X=eye(N,'like',L);X=repmat(X,[1 1 M(3)]);end
X=permute(X,[3 2 1]);
L=permute(L,[3 1 2]);
for n=2:N
    vn=1:n-1;
    X(:,:,n)=X(:,:,n)-sum(bsxfun(@times,L(:,n,vn),X(:,:,vn)),3);
end
X=permute(X,[3 2 1]);
