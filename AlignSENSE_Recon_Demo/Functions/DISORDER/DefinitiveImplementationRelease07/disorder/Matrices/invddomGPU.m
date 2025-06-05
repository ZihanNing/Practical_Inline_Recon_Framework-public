function Xp=invddomGPU(X,tol,nmax)

%INVDDOMGPU   Computes the approximate inverse of strictly diagonally 
%dominant matrices by using a Neumann series.
%   X=INVDDOMGPU(X,{TOL},{NMAX},{I})
%   * X are the set of input matrices
%   * {TOL} is a tolerance for inversion convergence. It refers to the mean 
%   square error of |I-Xi*X|. Defaults to 0 to avoid potentially costly 
%   convergence checks
%   * {NMAX} are the maximum number of iterations for inversion
%   convergence. Defaults to 1. With default tol one would get an 
%   approximation of the inverse without performing matrix multiplications.
%   Note series is X^{-1}=\sum_{k=0->inf} (I-D^{-1}X)^{k}D^{-1}
%   * Xp are the inverses
%

if nargin<2 || isempty(tol);tol=0;end
if nargin<3 || isempty(nmax);nmax=1;end

N=size(X);N(end+1:3)=1;
Dp=diagm(X);
Dp=1./Dp;%We assume no elements are zero as the matrices are strictly diagonally dominant
Xp=diagm(Dp);
if nmax==0;return;end
I=eye(N(1),'like',Dp);
E=bsxfun(@minus,I,bsxfun(@times,matfun(@transpose,Dp),X));
Ea=E;
if tol~=0%If tol is 0 we rely on a fixed number of iterations
    tol2=tol*tol;
    indConv=false([1 1 N(3)]);    
    e=abs(E);
    e=gather(multDimMea(e.*e,1:2));
    indConv(e<tol2)=true;
    if all(indConv(:));return;end

    for n=1:nmax
        Nr=sum(~indConv);
        [Eai,Ei,Xi,Xpi,Dpi,indConvi]=parUnaFun({Ea,E,X,Xp,Dp,indConv},@dynInd,~indConv,3);    
        Xpi=Xpi+bsxfun(@times,Eai,Dpi);    
        Xp=dynInd(Xp,~indConv,3,Xpi);
        ei=matfun(@mtimes,Xpi,Xi)-I(:,:,1:Nr);
        ei=abs(ei);
        ei=gather(multDimMea(ei.*ei,1:2));
        indConvi(ei<tol2)=true;
        if all(indConvi(:)) || n==nmax;return;end
        Eai=matfun(@mtimes,Eai,Ei);
        Ea=dynInd(Ea,~indConv,3,Eai);indConv=dynInd(indConv,~indConv,3,indConvi);
    end
else    
    for n=1:nmax
        Xp=Xp+bsxfun(@times,Ea,Dp);%Note that with a single iteration the inverse is approximated without any matrix multiplication as diagonal multiplications can be performed without summation
        if n==nmax;return;end
        Ea=matfun(@mtimes,Ea,E);
    end
end
        
        
