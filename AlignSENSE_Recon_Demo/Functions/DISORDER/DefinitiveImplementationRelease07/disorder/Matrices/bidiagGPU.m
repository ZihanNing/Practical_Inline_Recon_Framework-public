function [U,X,V]=bidiagGPU(X,typ)

%BIDIAGGPU   Performs a upper bidiagonalization the GPU. Good performance 
%probably only for small matrices. The code is based on
%https://uk.mathworks.com/matlabcentral/fileexchange/3568-bidiag?requestedDomain=true
%Paul Godfrey
%pgodfrey@intersil.com
%Intersil Corp.
%8-12-2002
%   [U,X,V]=BIDIAGGPU(X,{TYP})
%   * X is a matrix
%   * {TYP} indicates the mode of the decomposition, if 0 the bidiagonal 
%   matrix is real (default), else if may be complex for complex input
%   * U is a unitary matrix
%   * X is a bidiagonal matrix
%   * V is a unitary matrix 
%

if ~exist('typ','var') || isempty(typ);typ=0;end
gpu=isa(X,'gpuArray');
if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end

N=size(X);
Nmin=min(N(1:2));

if isreal(X);rea=1;else rea=0;end

U=dynInd(X,1,2);U(:)=1;U=permute(U,[2 1 3]);U=diagm(U);
V=dynInd(X,1,1);V(:)=1;V=diagm(V);
Ue=U;Ve=V;
for n=1:Nmin
    v11=n:N(1);v12=n;u1=1:N(1)-n+1;
    v21=n;v22=n+1:N(2);u2=1:N(2)-n;     
    if n<=N(1)%zero a col
        [v,Ui]=zerod(1,v11,v12,Ue(u1,u1,:));
        if ~typ && ~rea;Ui(1,:,:)=bsxfun(@times,sign(conj(v(1,1,:))),Ui(1,:,:));end
        U(v11,:,:)=matfun(@mtimes,Ui,U(v11,:,:));
        X(v11,:,:)=matfun(@mtimes,Ui,X(v11,:,:));
    end    
    if n<N(2)%zero a row
        [v,Vi]=zerod(2,v21,v22,Ve(u2,u2,:));
        Vi=conj(Vi);
        if ~typ && ~rea;Vi(:,1,:)=bsxfun(@times,sign(conj(v(1,1,:))),Vi(:,1,:));end
        V(:,v22,:)=matfun(@mtimes,V(:,v22,:),Vi);
        X(:,v22,:)=matfun(@mtimes,X(:,v22,:),Vi);
    end
end
U=matfun(@ctranspose,U);
if ~typ
    X=real(X);
    x=X(1:Nmin,1:Nmin,:);
    x=sign(diagm(x));
    x(x==0)=1;
    if N(1)<=N(2)
        U=bsxfun(@times,U,x);
        X=bsxfun(@times,matfun(@transpose,x),X);
    else
        V=bsxfun(@times,x,V);
        X=bsxfun(@times,X,x);
    end
end

function [v,UVi]=zerod(d,v1,v2,UVe)  
    x=X(v1,v2,:);v=x;
    if d==2;v=permute(x,[2 1 3]);end
    inda=find(x(1,1,:)==0);
    if ~isempty(inda);x(1,1,inda)=epsU;end%what if x is a zero vector or has x(1)=0?       
    xa=abs(x);
    v(1,1,:)=v(1,1,:)+sign(x(1,1,:)).*sqrt(sum(xa.*xa,d));
    va=abs(v);
    vg=sqrt(sum(va.*va,1)/2);
    v=bsxfun(@rdivide,v,vg);  
   
    UVi=UVe-bsxfun(@times,v,matfun(@ctranspose,v));
end

end
