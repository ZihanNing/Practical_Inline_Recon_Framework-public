function x=solveXMS2D(x,y,W,T,H,Hsmooth,S,Precond,Ak,xkGrid,kGrid,nX,NzSlab,outlD,BlSz,Ma)

%SOLVEXMS2D reconstructs an image under motion
%   X=SOLVEXMS2D(X,Y,W,T,H,HSMOOTH,S,SH,PRECOND,AK,XKGRID,KGRID,NX,TOLER,NZSLAB,GPU,OUTLD,BLSZ)
%   computes the best x for a given T
%   * X is the input image
%   * Y is the measured data
%   * W is a spatial mask
%   * T are the transform parameters
%   * H is the slice profile filter
%   * HSMOOTH is the regularization kernel
%   * S are the sensitivities
%   * PRECOND is the preconditioner
%   * AK is a sampling mask in the phase encoding direction
%   * XKGRID is a grid of points in the spatio-spectral domain
%   * KGRID is a grid of points in the spectral domain
%   * NX is the maximum number of iterations of CG method
%   * NZSLAB is the number of slices of the slab
%   * OUTLD is a mask for shot rejection
%   * BLSZ is the size of blocks for gpu processing
%   ** X is the reconstructed image
%

gpu=isa(x,'gpuArray');

N=size(x);NT=size(T);
for p=1:2
    NRun(p)=ceil(NT(4+p)/BlSz(p));
    NRunRem(p)=mod(NT(4+p),BlSz(p));
    for s=1:NRun(p)
        if s~=NRun(p) || NRunRem(p)==0
            vS{p}{s}=(s-1)*BlSz(p)+1:s*BlSz(p);
        else
            vS{p}{s}=(s-1)*BlSz(p)+1:(s-1)*BlSz(p)+NRunRem(p);
        end
    end
end
 
if sum(T(:))~=0
    etDir=cell(NRun(1),NRun(2));etInv=cell(NRun(1),NRun(2));
    for s=1:NRun(1)
        for t=1:NRun(2)
            xkGridAux=xkGrid;
            xkGridAux{1}{3}=dynInd(xkGridAux{1}{3},vS{2}{t},6);
            xkGridAux{2}{2}=dynInd(xkGridAux{2}{2},vS{2}{t},6);
            etDir{s}{t}=precomputeFactors3DTransform(xkGridAux,[],kGrid,dynInd(T,{vS{1}{s},vS{2}{t}},5:6),1,0);
            etInv{s}{t}=precomputeFactors3DTransform(xkGridAux,[],kGrid,dynInd(T,{vS{1}{s},vS{2}{t}},5:6),0,0);             
        end
    end
end

%CG method
%Initialization
NS=size(S);NY=size(y);

AkOutlDisc=bsxfun(@times,Ak,outlD);

yEnd=zeros([N(1:2) NzSlab 1 1 N(3)],'like',y);
for s=1:NRun(1)
    if sum(T(:))~=0;yS=bsxfun(@times,y,dynInd(outlD,vS{1}{s},5));else yS=y;end 
    yS=bsxfun(@times,yS,dynInd(Ak,vS{1}{s},5));

    yS=ifftGPU(yS,1,[],NS(1));
    %yS=ifold(yS,1,NS(1),NY(1));
    yS=sum(bsxfun(@times,yS,conj(S)),4);
    NZR=size(yS);NZR(end+1:5)=1;
    yZR{1}{2}=zeros(NZR,'like',yS);
    yZR{2}{1}=zeros([NZR(1:2) NzSlab NZR(4:5) NZR(3)],'like',yS);
    yS=extractSlabsOld(yS,NzSlab,1,0,yZR);
    yS=filtering(yS,H);
    if s==1
        NX=size(yS); 
        [F,FH]=buildStandardDFTM(NX,0,gpu);
    end
    
    for t=1:NRun(2)
        if sum(T(:))~=0
            etS=gpuTSt(etInv);              
            yEnd=dynInd(yEnd,vS{2}{t},6,dynInd(yEnd,vS{2}{t},6)+transform3DSinc(dynInd(yS,vS{2}{t},6),etS,0,F,FH));
        else
            yEnd=dynInd(yEnd,vS{2}{t},6,dynInd(yEnd,vS{2}{t},6)+sum(dynInd(yS,vS{2}{t},6),5));
        end
    end    
end
y=yEnd;
clear yEnd yS
NZR=size(y);
NZR(end+1:6)=1;
yZR{1}{1}=zeros([NZR(1:2) NZR(6)],'like',y);
yZR{2}{2}=zeros(NZR,'like',y);
y=extractSlabsOld(y,NzSlab,0,0,yZR);
y=W.*y;

Ap=applyCG(x);
r=y-Ap; 
z=Precond.*r;
p=z; 
rsold=multDimSum(conj(z).*r,1:3);

%Iterations
for n=1:nX
    Ap=applyCG(p);
    al=conj(rsold)/multDimSum(conj(p).*Ap,1:3);
    xup=al*p;
    x=x+xup;           
    r=r-al*Ap;
    z=Precond.*r;
    rsnew=multDimSum(conj(z).*r,1:3);
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew;
    if sqrt(abs(rsnew))<1e-10;break;end
end

function x=applyCG(x)   
    if ~isempty(Ma);xM=bsxfun(@times,x,Ma);end
    xB=filtering(x,Hsmooth);
    x=extractSlabsOld(x,NzSlab,1,1,yZR);
    NX=size(x);
    xEnd=zeros(NX,'like',x);
    for s=1:NRun(1)
        if sum(T(:))~=0;NXS=[NX(1:4) length(vS{1}{s}) NX(6)];else NXS=NX;end
        xS=zeros(NXS,'like',x);
        for t=1:NRun(2)
            if sum(T(:))~=0    
                etS=gpuTSt(etDir);
                xS=dynInd(xS,vS{2}{t},6,transform3DSinc(dynInd(x,vS{2}{t},6),etS,1,F,FH));
            else
                xS=dynInd(xS,vS{2}{t},6,dynInd(x,vS{2}{t},6));
            end                
        end
        xS=filtering(xS,H); 
        xS=extractSlabsOld(xS,NzSlab,0,1,yZR);
        xS=bsxfun(@times,xS,S);
        %xS=fold(xS,1,NS(1),NY(1));
        xS=fftGPU(xS,1,[],NY(1));    
        xS=bsxfun(@times,xS,AkOutlDisc(:,:,:,:,vS{1}{s}));   
        xS=ifftGPU(xS,1,[],NS(1));       
        %xS=ifold(xS,1,NS(1),NY(1));                     
        xS=sum(bsxfun(@times,xS,conj(S)),4);
        NZR=size(xS);NZR(end+1:5)=1;
        yZR{1}{2}=zeros(NZR,'like',xS);
        yZR{2}{1}=zeros([NZR(1:2) NzSlab NZR(4:5) NZR(3)],'like',xS);        
        xS=extractSlabsOld(xS,NzSlab,1,0,yZR);
        xS=filtering(xS,H);
        for t=1:NRun(2)
            if sum(T(:))~=0 
                etS=gpuTSt(etInv);               
                xEnd=dynInd(xEnd,vS{2}{t},6,dynInd(xEnd,vS{2}{t},6)+transform3DSinc(dynInd(xS,vS{2}{t},6),etS,0,F,FH)); 
            else
                xEnd=dynInd(xEnd,vS{2}{t},6,dynInd(xEnd,vS{2}{t},6)+sum(dynInd(xS,vS{2}{t},6),5));
            end
        end        
    end
    x=xEnd;    
    x=extractSlabsOld(x,NzSlab,0,0,yZR);
    x=x+xB;
    if ~isempty(Ma);x=x+xM;end
    x=x.*W;
end

function etS=gpuTSt(etIn)
    if gpu
        etS{1}=gpuArray(etIn{s}{t}{1});
        for m=1:3
            for l=2:3;etS{l}{m}=gpuArray(etIn{s}{t}{l}{m});end
        end
    else
        etS{1}=etIn{s}{t}{1};
        for m=1:3
            for l=2:3;etS{l}{m}=etIn{s}{t}{l}{m};end
        end
    end
end

end