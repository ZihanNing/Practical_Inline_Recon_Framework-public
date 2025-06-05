function [x,B]=TIEUnwrapping(x,voxSiz,weightFunc,typInp)

%TIEUNWRAPPING calls a phase unwrapping based on [1] Z Zhao et al. 
%"Comparative study of phase unwrapping algorithms based on solving the 
%Poisson equation," Meas Sci Technol, 31(065004):1-15, Mar 2020. Only
%alternatives based on the DCT are implemented
%   [X,B]=TIEUNWRAPPING(X,{VOXSIZ},{WEIGHTFUNC},{TYPINP})
%   * X is the complex data used to estimate the field
%   * {VOXSIZ} are the voxel sizes of the data
%   * {WEIGHTFUNC} weighting function
%   * {TYPINP} is the type of input function, either LSNONIT, TIENONIT,
%   LSIT, TIEIT respectively for LS and TIE non iterative and iterative, it
%   defaults to LSIT. Note iterative returns a multiple of 2*pi phase based
%   on [2] Z Zhao et al. "Robust 2D phase unwrapping algorithm based on the 
%   transport of intensity equation", Meas Sci Technol, 30(015201):1-8, 
%   Nov 2018.
%   * X is the estimation of the unwrapped field
%   * B is the distance between the unwrapped and the wrapped field
%

if nargin<2 || isempty(voxSiz)
    ND=numDims(x);
    voxSiz=ones(1,ND);
elseif numDims(voxSiz)==0
    ND=numDims(x);
else
    ND=length(voxSiz);
end
if nargin<3 || isempty(weightFunc);weightFunc='Magnitude';end
if nargin<4 || isempty(typInp);typInp='LSIt';end

gpu=isa(x,'gpuArray');

%UNWRAPPING PARAMETERS
w=voxSiz;w(end+1:ND)=1;
%fprintf('Used weights:%s\n',sprintf(' %.2f',w));

%Precondition
A.mi=ones(1,ND);
N=size(x);N=N(1:ND);
A.Ps=buildFilter(N*2,'2ndFiniteDiscrete',w,0,0,A.mi);
if gpu;A.Ps=gpuArray(A.Ps);end
A.Ps=real(A.Ps);
A.Ps(1)=1e9;%For inversion
A.Ps=A.Ps.^(-1);

NO=size(x);
x=resSub(x,ND+1:numDims(x));
NL=size(x,ND+1);
B=angle(x);x=abs(x);
for l=1:NL
    %xl=dynInd(x,l,ND+1);
    Bwrap=dynInd(B,l,ND+1);

%     for n=1:ND
%         NN=size(xl,n);
%         padEl=zeros(1,ND+1);padEl(n)=1;
%         Baux=diff(Bl,1,n);    
%         Baux=padarray(wrapToPi(Baux)/w(n),padEl,0,'post');    
%         Bauxabs=ones(1,'like',xl);
%         if ~isempty(strfind(weightFunc,'Magnitude'));Bauxabs=bsxfun(@times,Bauxabs,sqrt(xl.*dynInd(xl,[2:NN 1],n)));end     
%         if ~isempty(strfind(weightFunc,'Gradient'));Bauxabs=bsxfun(@times,Bauxabs,abs(cos(angle(exp(1i*B).*exp(-1i*dynInd(B,[2:NN 1],n))))));end%%%UNSURE ABOUT THIS, SEE OPTIONS IN SRNCPUNWRAPPING    
% 
%         if n==1
%             Bn=Baux;
%             EH.Mb=Bauxabs;
%         else
%             Bn=cat(ND+1,Bn,Baux);
%             EH.Mb=cat(ND+1,EH.Mb,Bauxabs);
%         end   
%     end
% 
%     E.Ff.Fd=1;EH.Fb.Fd=1;
%     E.Ff.w=w;EH.Fb.w=w;
%     E.Ff.ND=ND;EH.Fb.ND=ND;
%     E.Sm=0;
%     if E.Sm;Bn=gather(Bn);end
% 
%     %PHASE UNWRAPPING SOLVER
%     Bn=CGsolver(Bn,E,EH,A);
%     Bn=real(Bn);

    if strcmp(typInp,'TIENonIt') || strcmp(typInp,'LSNonIt')
        xl=real(filtering(inputPoisson(Bwrap,typInp,ND),A.Ps,1));        
    else   
        [B1,K1,xl,res,resn1]=unwrapRes(Bwrap,Bwrap);        
        nn=1;
        while 1
            [B1,K2,xl,res,resn2]=unwrapRes(res,Bwrap,B1);
            nn=nn+1;
            if all(K2(:)==K1(:)) || abs(resn1-resn2)<1e-5;break;end             
            K1=K2;resn1=resn2;
        end
        fprintf('Number of iterations: %d\n',nn);
    end
    x=dynInd(x,l,ND+1,xl);
end
x=reshape(x,NO);B=reshape(B,NO);
if nargout>=2;B=angle(exp(1i*(x-B)));end

    function [B,K,x,res,resn]=unwrapRes(res,Bwrap,B0)               
        B=real(filtering(inputPoisson(res,typInp,ND),A.Ps,1));
        if nargin>=3;B=B+B0;end
        B=B+mean(Bwrap(:))-mean(B(:));
        K=round((B-Bwrap)/2/pi);  %calculate integer K
        x=Bwrap+2*K*pi;
        res=wrapToPi(x-B);
        resn=sqrt(normm(res));
    end
end

function xou=inputPoisson(xin,typInp,ND)
    xou=xin;xou(:)=0;
    if strcmp(typInp,'TIENonIt') || strcmp(typInp,'TIEIt');xin=exp(1i*xin);end
    for n=1:ND
        padEl=zeros(1,ND);padEl(n)=1;
        if strcmp(typInp,'TIENonIt') || strcmp(typInp,'TIEIt');xaux=angle(diff(xin,1,n));else xaux=wrapToPi(diff(xin,1,n));end
        %xaux=wrapToPi(diff(xin,1,n));    
        xaux=padarray(xaux,padEl,0,'both');
        xaux=diff(xaux,1,n);       
        xou=xou+xaux;
    end
    if strcmp(typInp,'TIENonIt') || strcmp(typInp,'TIEIt');xou=imag(conj(xin).*xou);end
end

% function x=CGsolver(y,E,EH,A)    
%     nX=200;
%     toler=1e-3;
%     n=0;
%     r=decode(y,EH,E);n=n+1;
%     
%     NX=size(r);ND=numDims(r);
%     x=zeros(NX,'like',r);
%     z=precondition(r,A);
%     p=z;
%     rsold=multDimSum(conj(z).*r,1:ND);
%     l=true;
%     if sqrt(min(abs(rsold(:))))<1e-9;l=false;end
% 
%     %ITERATIONS
%     while l
%         %SYSTEM MATRIX
%         Ap=decode(encode(p,E),EH,E);n=n+2;
%     
%         %UPDATES
%         al=conj(rsold)./multDimSum(conj(p).*Ap,1:ND);
%         xup=bsxfun(@times,al,p);
%         x=x+xup;
%         xup=max(abs(xup(:)).^2);
%         n
%         xup
%         %visSegment(x)
%         if xup<toler || n>=nX;break;end
%         r=bsxfun(@minus,r,bsxfun(@times,al,Ap));
%         z=precondition(r,A);
% 
%         rsnew=multDimSum(conj(z).*r,1:ND);
%         if sqrt(min(abs(rsnew(:))))<1e-9;break;end
%         be=bsxfun(@times,rsnew,1./rsold);
%         p=z+bsxfun(@times,be,p);
%         rsold=rsnew;
%     end
% end
