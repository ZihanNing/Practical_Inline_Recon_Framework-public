function [U,X,V]=qrGPU(X,typ)

%QRGPU   Performs a upper triangulation / bidiagonalization using GPU code. Good performance 
%probably only for small matrices. The code is loosely based on
%https://uk.mathworks.com/matlabcentral/fileexchange/3568-bidiag?requestedDomain=true
%Paul Godfrey
%pgodfrey@intersil.com
%Intersil Corp.
%8-12-2002
%   [U,X,V]=QRALTGPU(X,{TYP})
%   * X is a matrix
%   * {TYP} indicates the mode of the decomposition, if first index equals 0 (default) it 
%   performs QR decomposition, if it equals 1 it performs QR with minimum fill-in, so 
%   VX=UR, if 2, it performs bidiagonalization. If second index equals 0 the R diagonals 
%   (and in the case of diagonalization upper diagonal) are real (default), else it may be 
%   complex for complex input. Minimum fill-in is based on the strategy in [1] P Businger and 
%   GH Golub, "Linear Least Squares Solutions by Householder Transformations," Num Mat, 
%   7:269-276, 1965
%   * U is a unitary matrix
%   * X is a bidiagonal matrix
%   * V is a unitary matrix or a column permutation array 
%

if nargin<2 || isempty(typ);typ=zeros(1,2);end
typ(end+1:2)=0;
gpu=isa(X,'gpuArray');
if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end
%dev=gpuDevice;

N=size(X);N(end+1:3)=1;
Nmin=min(N(1:2));

if isreal(X);rea=1;else rea=0;end

if typ(1)==1
    Vp=1:N(2);
    if gpu;Vp=gpuArray(Vp);end
    Vp=repmat(Vp,[1 1 N(3)]);
end

U=dynInd(X,1,2);U(:)=1;U=permute(U,[2 1 3]);U=diagm(U);Ue=U;
V=dynInd(X,1,1);
if typ(1)==1;V(:)=0;else V(:)=1;end
V=diagm(V);
if typ(1)==2;Ve=V;end
for n=1:Nmin
    v11=n:N(1);v12=n;u1=1:N(1)-n+1;
    if n<=N(1)%zero a col          
        if typ(1)==1%Pivoting, column with the largest norm (for the elements that will form the diagonal)
            p=pivot(X(v11,n:N(2),:),1);            
            indToFlip=p~=n;
            if any(indToFlip(:))                
                [Xi,Vpi,pi]=parUnaFun({X,Vp,p},@dynInd,indToFlip,3);
                np=repmat(n,size(pi));
                npp=horzcat(np,pi);
                nppi=indCom(Xi,npp,2);
                Xi(nppi)=Xi(flip(nppi,2));
                nppi=indCom(Vpi,npp,2);
                Vpi(nppi)=Vpi(flip(nppi,2));
                X=dynInd(X,indToFlip,3,Xi);Vp=dynInd(Vp,indToFlip,3,Vpi); 
            end
        end
        [v,Ui]=zerod(X(v11,v12,:),1,Ue(u1,u1,:));%This is particularly slow, not sure why
        if ~typ(2) && ~rea;Ui(1,:,:)=bsxfun(@times,sign(conj(v(1,1,:))),Ui(1,:,:));end
        U(v11,:,:)=matfun(@mtimes,Ui,U(v11,:,:));
        X(v11,:,:)=matfun(@mtimes,Ui,X(v11,:,:));
    end
    if typ(1)==2
        v21=n;v22=n+1:N(2);u2=1:N(2)-n;        
        if n<N(2)%zero a row       
            [v,Vi]=zerod(X(v21,v22,:),2,Ve(u2,u2,:));
            Vi=conj(Vi);
            if ~typ(2) && ~rea;Vi(:,1,:)=bsxfun(@times,sign(conj(v(1,1,:))),Vi(:,1,:));end
            V(:,v22,:)=matfun(@mtimes,V(:,v22,:),Vi);
            X(:,v22,:)=matfun(@mtimes,X(:,v22,:),Vi);
        end
    end
end
U=matfun(@ctranspose,U);
if ~typ(2)
    if typ(1)==2;X=real(X);end
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
if typ(1)~=2;X=triuGPU(X);end%Zeroing elements below the main diagonal
if typ(1)==1
    nppi=indCom(V,Vp,1);%Conversion to matrix form
    V(nppi)=1;
end

function [v,UVe]=zerod(x,d,UVe)
    x(1,1,x(1,1,:)==0)=epsU;%what if x is a zero vector or has x(1)=0?
    v=x;
    if d==2;v=permute(x,[2 1 3]);end    
    xa=abs(x);
    v(1,1,:)=v(1,1,:)+sign(x(1,1,:)).*sqrt(sum(xa.*xa,d));
    va=abs(v);
    vg=sqrt(sum(va.*va,1)/2);
    v=bsxfun(@rdivide,v,vg);
    UVe=UVe-bsxfun(@times,v,matfun(@ctranspose,v));
end

function p=pivot(x,d)
    x=abs(x);
    x=sum(x.*x,d);
    [~,p]=max(x,[],3-d);
    p=p+n-1;
end

end
