function [x,B]=CNCGUnwrapping(x,voxSiz,weightFunc,typInp)

%CNCGUNWRAPPING calls an 2-norm phase unwrapping (Convex Norm Conjugate
%Gradient) based on [1] DC Ghiglia, LA Romero. "Minimum Lp-norm 
%two-dimensional phase  unwrapping," J Opt Soc Am A, 13(10):1999-2013, Oct 
%1999, which is pre-weighted using the magnitude data as proposed in [2] DC 
%Ghiglia and LA Romero. "Robust two-dimensional weighted and unweighted
%phase unwrapping that uses fast transforms and iterative methods," J Opt
%Soc Am A, 11(1):107-117, Jan 1994, for dynamic field mapping
%   [X,B]=CNCGUNWRAPPING(X,{VOXSIZ},{WEIGHTFUNC},{TYPINP})
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

%preconditionU
A.mi=ones(1,ND);
N=size(x);N=N(1:ND);
A.Ps=-buildFilter(N*2,'2ndFiniteDiscrete',w,0,0,A.mi);
if gpu;A.Ps=gpuArray(A.Ps);end
A.Ps=real(A.Ps);
A.Ps(1)=1e9;%For inversion
A.Ps=A.Ps.^(-1);

NO=size(x);
x=resSub(x,ND+1:numDims(x));
NL=size(x,ND+1);
B=angle(x);x=abs(x);
for l=1:NL
    xl=dynInd(x,l,ND+1);
    Bl=dynInd(B,l,ND+1);
    NNX=numel(xl);

    for n=1:ND
        NN=size(xl,n);   
        Bauxabs=ones(1,'like',xl);
        padEl=zeros(1,ND+1);padEl(n)=1;
        if ~isempty(strfind(weightFunc,'Magnitude'));Bauxabs=bsxfun(@times,Bauxabs,sqrt(xl.*dynInd(xl,[2:NN 1],n)));end     
        %if
        %~isempty(strfind(weightFunc,'Gradient'));Bauxabs=bsxfun(@times,Bauxabs,abs(cos(angle(exp(1i*B).*exp(-1i*dynInd(B,[2:NN 1],n))))));end%%%UNSURE ABOUT THIS, REPLACED AS BELOW
        if ~isempty(strfind(weightFunc,'Gradient1'));Bauxabs=Bauxabs.*padarray(abs(cos(wrapToPi(diff(Bl,1,n))/2)),padEl,0,'post');end
        if ~isempty(strfind(weightFunc,'Gradient2'));Bauxabs=Bauxabs.*padarray(min(1./(1-abs(cos(wrapToPi(diff(Bl,1,n))/2))),1e3),padEl,0,'post');end
        if ~isempty(strfind(weightFunc,'Gradient3'));Bauxabs=Bauxabs.*padarray(max(cos(wrapToPi(diff(Bl,1,n))),1e-3),padEl,0,'post');end%Does not converge
        if ~isempty(strfind(weightFunc,'Gradient4'));Bauxabs=Bauxabs.*padarray(min(1./(1-max(cos(wrapToPi(diff(Bl,1,n))),1e-3)),1e3),padEl,0,'post');end


        if n==1;EH.Mb=Bauxabs;else EH.Mb=cat(ND+1,EH.Mb,Bauxabs);end    
    end

    E.Ff.Fd=1;EH.Fb.Fd=1;
    E.Ff.w=w;EH.Fb.w=w;
    E.Ff.ND=ND;EH.Fb.ND=ND;
    E.Sm=0;

    %PHASE UNWRAPPING SOLVER    
    if strcmp(typInp,'LSNonIt')
        Bn=inputPoisson(Bl,ND,w);
        if E.Sm;Bn=gather(Bn);end
        Bn=real(CGsolver(Bn,E,EH,A));   
    else
        [B1,K1,Bn,res1,resn1]=unwrapRes(Bl,Bl);        
        nn=1;
        while 1
            [B1,K2,Bn,res2,resn2]=unwrapRes(res1,Bl,B1);            
            %fprintf('Iteration: %d / Max diff residual: %.2f / Mea diff residual: %.2f / Diff norm: %.2f / DiffPoints: %d\n',nn,max(abs(res1(:)-res2(:))),sqrt(normm(res2,res1)),abs(resn1-resn2),sum(K2(:)~=K1(:)));
            %pause
            nn=nn+1;
            if all(K2(:)==K1(:)) || abs(resn1-resn2)<1e-6*NNX || nn==100;break;end             
            K1=K2;resn1=resn2;res1=res2;
        end
        %fprintf('Number of iterations: %d\n',nn);
    end
    x=dynInd(x,l,ND+1,Bn);
end
x=reshape(x,NO);B=reshape(B,NO);
if nargout>=2;B=angle(exp(1i*(x-B)));end
    function [B,K,x,res,resn]=unwrapRes(res,Bwrap,B0)               
        B=inputPoisson(res,ND,w);
        if E.Sm;B=gather(B);end
        B=real(CGsolver(B,E,EH,A));
        if nargin>=3;B=B+B0;end
        %B=B+mean(Bwrap(:))-mean(B(:));
        %B=B+mean(abs(xl(:)).^2.*Bwrap(:))-mean((abs(xl(:)).^2).*B(:));
        K=round((B-Bwrap)/2/pi);  %calculate integer K
        x=Bwrap+2*K*pi;
        res=wrapToPi(x-B);
        resn=sqrt(normm(res));
    end
end

function x=CGsolver(y,E,EH,A)    
    nX=200;
    toler=2e-3;
    n=0;
    r=decodeU(y,EH,E);n=n+1;
    
    NX=size(r);ND=numDims(r);
    x=zeros(NX,'like',r);
    z=preconditionU(r,A);
    p=z;
    rsold=multDimSum(conj(z).*r,1:ND);
    l=true;
    if sqrt(min(abs(rsold(:))))<1e-9;l=false;end

    %ITERATIONS
    while l
        %SYSTEM MATRIX
        Ap=decodeU(encodeU(p,E),EH,E);n=n+2;
    
        %UPDATES
        al=conj(rsold)./multDimSum(conj(p).*Ap,1:ND);
        xup=bsxfun(@times,al,p);
        x=x+xup;
        xup=max(abs(xup(:)).^2);
        %visSegment(x)
        if xup<toler || n>=nX;break;end
        r=bsxfun(@minus,r,bsxfun(@times,al,Ap));
        z=preconditionU(r,A);

        rsnew=multDimSum(conj(z).*r,1:ND);
        if sqrt(min(abs(rsnew(:))))<1e-9;break;end
        be=bsxfun(@times,rsnew,1./rsold);
        p=z+bsxfun(@times,be,p);
        rsold=rsnew;
    end
end

function xou=inputPoisson(xin,ND,w)
    for n=1:ND      
        padEl=zeros(1,ND+1);padEl(n)=1;
        Baux=diff(xin,1,n);    
        Baux=padarray(wrapToPi(Baux)/w(n),padEl,0,'post');            
        if n==1;xou=Baux;else xou=cat(ND+1,xou,Baux);end
    end
end
