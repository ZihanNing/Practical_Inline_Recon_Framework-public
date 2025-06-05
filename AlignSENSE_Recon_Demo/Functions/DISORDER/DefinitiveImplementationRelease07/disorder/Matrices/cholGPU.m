function [L,r]=cholGPU(X,tol)

%CHOLGPU   Performs multiple full rank matrix lower Cholesky 
%decompositions in parallel. Does not check for Hermitian symmetry. 
%The code is based on the extension of [1] P Courrieu, "Fast 
%Computation of Moore-Penrose Inverse Matrices," Neural Inf Proc - 
%Letters and Reviews, 8(2):25-29, 2005 for multiple pages. Probably 
%good performance only for small matrices
%   [L,R]=CHOLGPU(X,{TOL})
%   * X is the input matrix
%   * {TOL} is a given tolerance for rank deficient matrices
%   * L is the lower triangular output matrix
%   * R is an estimate for the rank
%

%dev=gpuDevice;

gpu=isa(X,'gpuArray');
if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end

N=size(X);N(end+1:3)=1;
assert(N(1)==N(2),'Cholesky only defined for squared matrices while input is %d x %d',N(1),N(2));

if nargin<2 || isempty(tol);tol=epsU*sqrt(max(diagm(X),[],2))*N(1);end%Similar to SVD
X=permute(X,[3 1 2]);%Third dimension along the first 
N=size(X);
L=zeros(N,'like',X);
tol=permute(tol,[3 1 2]);
tol2=tol.*tol;
r=ones([N(1) 1 1],'like',X);r=real(r);
N=N(2);
vnc=1:N;
L(:,:,1)=X(:,:,1);La=L(:,:,1);
for n=1:N
    vnl=1:n-1;%Maximum possible rank      
    np=indCom(L,r,3);
    npvnc=np(:,vnc,:);
    if n>1
        La=zrGPU(L(:,vnc,vnl),r-1,3);
        La=X(:,vnc,n)-sum(bsxfun(@times,La,conj(La(:,1,:))),3);
    end
    Lap=La(:,1,:);
    upTol=single(Lap>tol2);
    L(npvnc)=bsxfun(@times,1-upTol,La)+bsxfun(@times,upTol,bsxfun(@rdivide,La,max(sqrt(Lap),tol)));
    r=r+upTol;
    vnc(1)=[];
end
r=r-1;
L=permute(L,[2 3 1]);
r=permute(r,[2 3 1]);
L=matfun(@ctranspose,L);%We are interested in upper triangular convention
