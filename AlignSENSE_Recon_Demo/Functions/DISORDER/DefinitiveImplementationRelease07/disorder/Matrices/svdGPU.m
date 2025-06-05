function [U,X,V]=svdGPU(X,tol)

%SVDGPU   Performs a SV decomposition on the GPU. The code is based on
%bidiagonalization first and then successive Givens rotations to produce
%the final diagonalization
%   [U,X,V]=SVDGPU(X,{TOL})
%   * X is a matrix
%   * {TOL} is the error tolerance
%   * U is a unitary matrix of left singular vectors
%   * X is a diagonal matrix
%   * V is a unitary matrix of right singular vectors
%

if nargin<2 || isempty(tol);tol=1e-4;end
maxNit=200;
N=size(X);
if N(1)<N(2);flipped=1;else flipped=0;end
if flipped;X=matfun(@ctranspose,X);end

%BIDIAGONALIZE
[U,X,V]=qrGPU(X,2);
%DIAGONALIZE
N=size(X);N(end+1:3)=1;
if N(2)>1
    convSVD=false([1 1 N(3)]);
    Xe=X(1:N(2),:,:);
    d=diagm(Xe);e=diagm(Xe,1);%THERE IS A PROBLEM HERE, IT IS NOT WORKING
    for n=1:maxNit
        [ei,di,Ui,Vi]=parUnaFun({e,d,U,V},@dynInd,~convSVD,3);
        [di,ei,Ui,Vi]=zbdsqrGPU(di,ei,Ui,Vi,tol);
        d=dynInd(d,~convSVD,3,di);e=dynInd(e,~convSVD,3,ei);U=dynInd(U,~convSVD,3,Ui);V=dynInd(V,~convSVD,3,Vi);
        em=max(abs(e),[],2);        
        convSVD(em<tol)=true;
        if all(convSVD(:));break;end
    end
    X(:)=0;X(1:N(2),:,:)=diagm(d);
end

if flipped
    X=matfun(@transpose,X);
    [V,U]=parUnaFun({U,V});
end
