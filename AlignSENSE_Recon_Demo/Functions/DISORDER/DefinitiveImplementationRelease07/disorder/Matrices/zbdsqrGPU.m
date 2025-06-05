function [d,e,U,V]=zsQRGPU(d,e,U,V,tol)

%ZSQRGPU   Diminishes the strength of bidiagonal matrices (generally as 
%part of the SVD algorithm). Good performance probably only for small 
%matrices. The code is based on an example at
%http://www.math.pitt.edu/~sussmanm/2071/lab09/index.html
%which was in turn based on
% Demmel & Kahan zero-shift QR downward sweep
%   [D,E,U,V]=ZSQRGPU(D,E,U,V,{TOL})
%   * D are the elements on the diagonal of the bidiagonal form
%   * E are the elements on the subdiagonal-1 of the bidiagonal form
%   * U is a unitary matrix
%   * V is a unitary matrix
%   * {TOL} is the error tolerance
%   * D are the elements on the diagonal of the bidiagonal form
%   * E are the elements on the subdiagonal-1 of the bidiagonal form
%   * U is a unitary matrix
%   * X is a diagonal matrix
%   * V is a unitary matrix 
%

if ~exist('tol','var') || isempty(tol);tol=1e-4;end
gpu=isa('d','gpuArray');
dev=gpuDevice;
%wait(dev);tic
N=size(d,2);Nmaxmax=N-1;
P=size(V,3);

%FIND THE SUPPORT
x=(abs(e)>=tol);
[Vma,Nin(1,1,:)]=max(x,[],2);
[~,Nin(2,1,:)]=max(flip(x,2),[],2);
Nin=single(Nin);%Single indexing may be problematic on very large (>1e8, perhaps) matrices
Nin(:,1,Vma==0)=N;
Nin(2,1,:)=N-Nin(2,1,:)+1;%Last diagonal element for which convergence has not been achieved
N=Nin(2,1,:)-Nin(1,1,:)+1;
%N(:)'

%HARD DECISION ON CONVERGED PAGES AND ASSURE THEY ARE NOT PASSED INTO AGAIN
convDia=false([1 1 P]);
if gpu;convDia=gpuArray(convDia);end
convDia(N<=1)=true;
e(convDia)=0;

%DECIDE ON FLIPPING AND COMPUTE MINIMUM EIGENVALUE OF LAST (OR FIRST) 2X2 BLOCK FOR SHIFTING
sh=d(1,1,:);sh(:)=0;fl=sh+1;%Not to flip by default
flv=repmat(sh,[2 1 1]);flv(1,1,:)=1;flv(2,1,:)=-1;%Not to flip / to flip
[di,ei,fli,flvi,Nini]=parUnaFun({d,e,fl,flv,Nin},@dynInd,~convDia,3);  
fli(indDim(di,Nini(2,1,:),2)>indDim(di,Nini(1,1,:),2))=2;%Whether to flip
flvi=indDim(flvi,fli,1); 
Ninr=indDim(Nini(:,1,:),3-fli,1);
Ninre=Ninr-(fli-1);
shi=sv2(indDim(di,Ninr-flvi,2),indDim(di,Ninr,2),indDim(ei,Ninre-flvi,2));%Shift value
sh=dynInd(sh,~convDia,3,shi);fl=dynInd(fl,~convDia,3,fli);
flv=indDim(flv,fl,1);

%SHIFTED SWEEP
Ninr=indDim(Nin(:,1,:),fl,1);
f=indDim(d,Ninr,2);
inda=f>1e-9;
if any(inda(:))
    inda=find(inda);
    dinda=indDim(d(1,:,inda),Ninr(1,1,inda),2);
    f(inda)=(abs(dinda)-sh(inda)).*(sign(dinda)+sh(inda)./dinda);
end
g=indDim(e,Ninr-fl+1,2);
indfl1=fli==1;
indfl2=fli==2;
%wait(dev);toc

indf=0:P-1;
if gpu;indf=gpuArray(indf);end
Ninr=double(Ninr);flv=double(flv);fl=double(fl);

for n=1:Nmaxmax
    %First we check what pages have not converged
    convDia(N<=n)=true;
    Nnoconv=sum(~convDia);
    if any(~convDia(:))
        convDiaF=find(~convDia);
        nc=Ninr+flv*(n-1);        
        [fi,gi,di,ei,Ni,nci,fli,flvi,Ui,Vi,indfl1i,indfl2i]=parUnaFun({f,g,d,e,N,nc,fl,flv,U,V,indfl1,indfl2},@dynInd,convDiaF,3);
        indfi=indf(1:Nnoconv);
        c=fi;s=fi;q=fi;
        
        %1D indexes
        ncie=nci-(fli-1);        
        nci2=horzcat(ncie,ncie+1);
        if any(indfl1i(:))
            indfl1i=find(indfl1i);
            nci2indfl1i=nci2(:,:,indfl1i);
            indfifl1i=indf(1:length(indfl1i));
        else
            nci2indfl1i=[];
        end
        if any(indfl2i(:));
            indfl2i=find(indfl2i);
            nci2indfl2i=nci2(:,:,indfl2i);
            indfifl2i=indf(1:length(indfl2i));
        else
            nci2indfl2i=[];
        end
        ncipflvi=nci+flvi;                
        
        %ND indexes  
        cncie=indCom(ei,ncie,2,[],indfi);       
        cnci=indCom(di,nci,2,[],indfi);
        cncipflvi=indCom(di,ncipflvi,2,[],indfi);        
        
        %wait(dev);tic
        givensrotcalc(fi,gi);
        %wait(dev);toc
        %wait(dev);tic
        if n~=1;ei(indCom(ei,ncie-flvi,2,[],indfi))=q;end
        dp=di(cnci);ep=ei(cncie);
        fi=c.*dp+s.*ep;
        ei(cncie)=c.*ep-s.*dp;
        dp=di(cncipflvi);
        gi=s.*dp;
        di(cncipflvi)=c.*dp;
        %wait(dev);toc

        %wait(dev);tic
        if ~isempty(nci2indfl1i)         
            A=Vi(:,:,indfl1i);
            cnci2ind=indCom(A,nci2indfl1i,2,[],indfifl1i);
            A(cnci2ind)=givensrotappl(c(indfl1i),s(indfl1i),A(cnci2ind));            
            Vi(:,:,indfl1i)=A;
        end
        if ~isempty(nci2indfl2i)
            A=Ui(:,:,indfl2i);
            cnci2ind=indCom(A,nci2indfl2i,2,[],indfifl2i);
            A(cnci2ind)=givensrotappl(c(indfl2i),-s(indfl2i),A(cnci2ind));
            Ui(:,:,indfl2i)=A;
        end
        %wait(dev);toc

        %wait(dev);tic
        givensrotcalc(fi,gi);
        %wait(dev);toc
        
        %wait(dev);tic
        di(cnci)=q;
        dp=di(cncipflvi);ep=ei(cncie);
        fi=c.*ep+s.*dp;
        di(cncipflvi)=c.*dp-s.*ep;
        inda=n<Ni-1;
        if any(inda(:))
            inda=find(inda);
            einda=ei(1,:,inda);
            ncienda=ncie(inda)+flvi(inda);
            cncienda=indCom(einda,ncienda,2);
            eindap=einda(cncienda);
            gi(inda)=s(inda).*eindap;
            einda(cncienda)=c(inda).*eindap;
            ei(1,:,inda)=einda;      
        end
        %wait(dev);toc
        
        %wait(dev);tic
        if ~isempty(nci2indfl1i)
            A=Ui(:,:,indfl1i);
            cnci2ind=indCom(A,nci2indfl1i,2,[],indfifl1i);
            A(cnci2ind)=givensrotappl(c(indfl1i),s(indfl1i),A(cnci2ind));
            Ui(:,:,indfl1i)=A;
        end
        if ~isempty(nci2indfl2i)
            A=Vi(:,:,indfl2i);
            cnci2ind=indCom(A,nci2indfl2i,2,[],indfifl2i);
            A(cnci2ind)=givensrotappl(c(indfl2i),-s(indfl2i),A(cnci2ind));
            Vi(:,:,indfl2i)=A;
        end
        %wait(dev);toc        
        %pause
        inda=n==Ni-1;
        if any(inda(:));ei(1,:,inda)=indDim(ei(1,:,inda),ncie(inda),2,fi(inda));end
        f=dynInd(f,convDiaF,3,fi);g=dynInd(g,convDiaF,3,gi);d=dynInd(d,convDiaF,3,di);e=dynInd(e,convDiaF,3,ei);U=dynInd(U,convDiaF,3,Ui);V=dynInd(V,convDiaF,3,Vi);
    end    
end

function G=givensrotappl(c,s,M)    
    cs=horzcat(c,s);
    NM=size(M,1);cs=repmat(cs,[NM 1 1]);
    M1=M(:,1,:);M2=M(:,2,:);
    G=horzcat(M1,-M1).*cs+horzcat(M2,M2).*cs(:,[2 1],:);
    %NM=size(M,1);c=repmat(c,[NM 1 1]);s=repmat(s,[NM 1 1]);
    %G=horzcat(M(:,1,:).*c+M(:,2,:).*s,-M(:,1,:).*s+M(:,2,:).*c);
end

function givensrotcalc(f,g)
    inda=(f==0);
    indb=abs(f)>abs(g) & ~inda;
    indc=~inda & ~indb; 
    if any(inda(:))
        inda=find(inda);
        c(inda)=0;
        s(inda)=1;
        q(inda)=g(inda);
    end
    if any(indb(:))
        indb=find(indb);
        t=g(indb)./f(indb);
        t1=sqrt(1+t.*t);
        c(indb)=1./t1;
        s(indb)=t.*c(indb);
        q(indb)=f(indb).*t1;
    end
    if any(indc(:))     
        indc=find(indc);
        t=f(indc)./g(indc);
        t1=sqrt(1+t.*t);
        s(indc)=1./t1;
        c(indc)=t.*s(indc);
        q(indc)=g(indc).*t1;
    end
end

function sh=sv2(dd,ee,ff)
    d1=dd.*dd+ff.*ff;%First diagonal element
    d2=ee.*ee;%Second diagonal element
    D=ff.*ee;%Cross diagonal elements
    D=d1.*d2-D.*D;%Determinant
    T=(d1+d2)/2;%Trace
    D=sqrt(max(T.*T-D,0));%Discriminant    
    sh=sqrt(T-D);%Smaller eigenvalue    
    sh(D==0)=0;
end 

end

%These are some lines to accelerate sweeping without shifting:
% f=d(1,1,:);oldc=f;oldc(:)=1;g=e(1,1,:);
% for n=1:N-1%One sweep approx before testing convergence externally
%     nc=n:n+1;    
%     [c,s,q]=givensrotcalc(f,g);
%     if n~=1;e(1,n-1,:)=olds.*q;end
%     f=oldc.*q;g=s.*d(1,n+1,:);h=c.*d(1,n+1,:);
%     V(:,nc,:)=givensrotappl(c,s,V(:,nc,:));
% 
%     [c,s,d(1,n,:)]=givensrotcalc(f,g);
%     if n<N-1;g=e(1,n+1,:);oldc=c;olds=s;f=h;end        
%     U(:,nc,:)=givensrotappl(c,s,U(:,nc,:));
% end
% e(1,N-1,:)=h.*s;
% d(1,N,:)=h.*c;
