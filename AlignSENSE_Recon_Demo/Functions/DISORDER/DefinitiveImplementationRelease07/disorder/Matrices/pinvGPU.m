function [X,U,S,V,r]=pinvGPU(X,tol,typ)

%PINVGPU   Estimates the pseudo-inverse of a set of matrices on the GPU. 
%Probably good performance only for small matrices. One method is based on 
%a COD as in [1] VN Katsikis, D Pappas, A Petralias, "An improved method 
%for the computation of the Moore?Penrose inverse matrix," Appl Math Comp, 
%217:9828-9834, 2011. Another version is based on [2] P Courrieu, "Fast 
%Computation of Moore-Penrose Inverse Matrices," Neural Inf Proc - Letters 
%and Reviews, 8(2):25-29, 2005, which grounds on a Cholesky factorization.
%Second method might be quicker but less accurate.
%   [X,U,S,V,R]=PINVGPU(X,TOL)
%   * X is a matrix
%   * {TOL} is the error tolerance
%   * {TYP} is the method to estimate the pseudo-inverse. Either 'Katsikis'
%   (default) or 'Courrieu'
%   * X is the pseudoinverse
%   * U is a unitary matrix of left singular vectors
%   * S is a diagonal matrix
%   * V is a unitary matrix of right singular vectors
%   * R are the estimated ranks
%

if nargin<3 || isempty(typ);typ='Katsikis';end

gpu=isa(X,'gpuArray');
if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end
%dev=gpuDevice;

N=size(X);N(end+1:3)=1;
Nmax=max(N(1:2));

if strcmp(typ,'Katsikis')
    %wait(dev);tic
    %COD
    %[U,X,V]=qrGPU(X,1);
    X=gather(X);
    U=zeros(N([1 1 3]),'like',X);
    V=zeros(N([2 2 3]),'like',X);
    if N(3)>5000
    parfor n=1:N(3);[U(:,:,n),X(:,:,n),V(:,:,n)]=qr(X(:,:,n));end%Twice quicker than GPU implementation
    else
    for n=1:N(3);[U(:,:,n),X(:,:,n),V(:,:,n)]=qr(X(:,:,n));end
    end
    if gpu;[U,X,V]=parUnaFun({U,X,V},@gpuArray);end
    %wait(dev);toc
    %TOLERANCE
    if ~exist('tol','var') || isempty(tol);tol=epsU*X(1,1,:)*Nmax;end%This is not exactly the same as for pinv. In general the tolerances will be differently scaled...
    %if ~exist('tol','var') || isempty(tol);tol=epsU*max(abs(D),[],2)*Nmax;end%with D the singular values. Think the first diagonal element of the qr decomposition will be smaller than the largest singular value

    %wait(dev);tic
    %RANK ESTIMATE AND ROW ZEROING
    r=sum(any(bsxfun(@gt,abs(X),tol),2),1);
    X=zrGPU(X,r);
    %wait(dev);toc
    
    %wait(dev);tic
    %GENERALIZED INVERSE OF RANK REVEALING UPPER TRIANGULAR FORM
    S=ginvGPU(X,r);
    %S*Xo*S-S
    %Xo*S*Xo-Xo%Right multiplication not accurate, as with pinv
    %S*Xo-(S*Xo)'
    %Xo*S-(Xo*S)'%Not accurate when doing truncation, while pinv is, key property probably
    %wait(dev);toc
    %U=zrGPU(U,r,2);

    %wait(dev);tic
    %PSEUDOINVERSE
    X=matfun(@mtimes,V,matfun(@mtimes,S,matfun(@ctranspose,U)));
    %wait(dev);toc
elseif strcmp(typ,'Courrieu')%Sometimes not well conditioning; still to be investigated
    S=[];U=[];V=[];
    %if N(1)<=N(2);A=matfun(@mtimes,X,matfun(@ctranspose,X));else A=matfun(@mtimes,matfun(@ctranspose,X),X);end
    if N(1)<=N(2);A=emtimes(X,matfun(@ctranspose,X));else A=emtimes(matfun(@ctranspose,X),X);end
    if ~exist('tol','var') || isempty(tol);tol=epsU*sqrt(max(A,[],2))*Nmax;end%Similar to SVD
    %FULL RANK CHOLESKY FACTORIZATION OF A
    if N(3)<10000;A=gather(A);end
    [L,r]=cholGPU(A,tol);
    if N(3)<10000 && gpu;[L,r]=parUnaFun({L,r},@gpuArray);end

    %ROW ZEROING
    L=zrGPU(L,r);
    
    %GENERALIZED INVERSE OF CHOLESKY FACTORIZATION
    M=ginvGPU(L,r);
    if N(1)<=N(2)
        X=emtimes(emtimes(matfun(@ctranspose,X),M),matfun(@ctranspose,M));
    else
        %X=emtimes(M,emtimes(matfun(@ctranspose,M),matfun(@ctranspose,X)));        
        X=emtimes(M,matfun(@ctranspose,emtimes(X,M)));
    end
else
    error('Undefined pseudoinversion method %s',typ);
end

