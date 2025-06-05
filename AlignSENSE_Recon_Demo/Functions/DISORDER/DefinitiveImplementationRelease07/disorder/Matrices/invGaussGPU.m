function X=invGaussGPU(X)

%INVGAUSSGPU   Computes the inverse of a set of square matrices by Gaussian 
%elimination. Pivoting is not implemented so, for stability, the matrix should 
%be positive definite
%   X=INVGAUSSGPU(X)
%   * X is the set of input matrices
%   * X are the inverses of the matrices
%

M=size(X);M(end+1:3)=1;
N=M(1);
X=permute(X,[3 1 2]);%We arrange the matrices along the first dimension
for n=1:N
    t=1./X(:,n,n);    
    X(:,:,n)=bsxfun(@times,t,X(:,:,n));
    s=-X(:,n,:);
    s(:,1,n)=0;
    X=X+bsxfun(@times,s,X(:,:,n));
    s(:,1,n)=1;
    X(:,n,:)=bsxfun(@times,t,s);
end
X=permute(X,[2 3 1]);
