function [x,T,outlD]=optimizeLevelMS2D(x,y,T,S,W,Ak,xGrid,parX,debug,res,estT,outlD,alpha,toler,Ma)

%OPTIMIZELEVELMS2D gets the aligned reconstruction for a given resolution
%level
%   [X,T,OUTLD]=OPTIMIZELEVELMS2D(X,Y,T,S,W,AK,XGRID,MNORM,PARX,DEBUG,RES,ESTT,GPU,SLTH,SLOV,OUTLD,ALPHA)
%   applies the motion corrected reconstruction algorithm for MS sequences
%   at a given resolution level.
%   X is the input image
%   Y is the measured data
%   T are the transform parameters
%   S are the sensitivities
%   W is a spatial mask
%   AK is a sampling mask in the phase encoding direction
%   XGRID is a grid of points in the spatial domain
%   PARX is a structure with the reconstruction parameters
%   DEBUG is a flag to print debug information: 0-> no; 1-> basic; 2-> 
%   verbose
%   RES is the subsampling factor for this resolution level
%   ESTT is a flag that determines whether to apply motion correction
%   GPU is a flag that determines whether to use gpu processing
%   OUTLD is a mask for shot rejection
%   ALPHA is the regularization factor to treat the slice profiles
%   It returns:
%   X, the reconstructed image
%   T, the estimated tranform paramaters
%   OUTLD, the estimated shot rejection mask

gpu=isa(x,'gpuArray');

%ROI computation
ROIsec=0;
ROI=computeROI(W,ROIsec);
if parX.threeD
    if parX.correct==0;NzSlab=3;else NzSlab=parX.NzSlab;end
else
    NzSlab=1;
end

%z-Overlap
H=sliceProfile(NzSlab,parX.SlTh,parX.SlOv,parX);
NW=size(W);NW(end+1:3)=1;
if ~parX.threeD;alpha=0;end
Hsmooth=buildFilter([1 1 NW(3)],'FractionalFiniteDiscrete',[1 1 1],gpu,parX.fractionReg);
Hsmooth=alpha*abs(Hsmooth).^2;

%ROI extraction in the readout direction
[y,S,W,x,xGrid{2}]=parUnaFun({y,S,W,x,xGrid{2}},@extractROI,ROI,1,2);
if ~isempty(Ma);Ma=extractROI(Ma,ROI,1,2);end

N=size(W);
N(3)=NzSlab;

kGrid{1}=permute(single(-floor(N(1)/2):ceil(N(1)/2)-1),[2 1])/res;
kGrid{2}=single(-floor(N(2)/2):ceil(N(2)/2)-1)/res;
kGrid{3}=permute(single(-floor(N(3)/2):ceil(N(3)/2)-1),[1 3 2]);

%%Update transform terms
for m=1:3;kGrid{m}=2*pi*kGrid{m}/N(m);end    
yZ{2}{2}=single(zeros([1 1 N(3) 1 1 size(xGrid{3},3)]));


xGrid{3}=extractSlabsOld(xGrid{3},N(3),1,1,yZ);

per(1,:)=[1 3 2];per(2,:)=[2 1 3];
xkGrid=cell(2,3);
for n=1:2
    for m=1:3
        xkGrid{n}{m}=bsxfun(@times,xGrid{per(3-n,m)},kGrid{per(n,m)});%First index denotes that the k takes the first dimension of xkgrid; second index denotes the rotation        
    end
end
fact(1,:)=[1 2 3 1 1 2];fact(2,:)=[1 2 3 2 3 3];
kkGrid=cell(1,6);
for m=1:6;kkGrid{m}=bsxfun(@times,kGrid{fact(1,m)},kGrid{fact(2,m)});end

%Initialize solver
winic=1;
%winic=1e-2;
NT=size(T); 
w=winic*ones(NT(1:6));
flagw=zeros(NT(1:6));
nExtern=1000;
nX=5;
nT=1;

%Preconditioner computation
if isempty(Ma)
    Precond=(sum(real(conj(S).*S),4)+1e-9).^(-1);
else
    Precond=(sum(real(conj(S).*S),4)+Ma).^(-1);
end

y=fftGPU(y,1);
mNorm=1;%max(abs(y(:)));
%fprintf('Normalization: %.2f\n',mNorm);%4800-5300-3500-6700-4800
y=y/mNorm;
x=x/mNorm;
for n=1:nExtern       
    %if debug==2;tstart=tic;end
    xant=x;
    %Solve for x
    %Block sizes for gpu computation
    if res==1 && sum(T(:))~=0    
        BlSz(1)=1;BlSz(2)=ceil(size(T,6)/2);
    else
        BlSz(1)=floor(size(T,5)/2);BlSz(2)=size(T,6);
        if BlSz(1)==0;BlSz(1)=1;end
    end        
    x=solveXMS2D(x,y,W,T,H,Hsmooth,S,Precond,Ak,xkGrid,kGrid,nX,NzSlab,outlD,BlSz,Ma);      
    %if debug==2;telapsed=toc(tstart);fprintf('Time solving x: %.4fs\n',telapsed);end
  
    %figure
    %imshow(squeeze(outlD),[])
  

    %visReconstruction(dynInd(x,88,3),0)
    %if n==1
    %visReconstruction(dynInd(Ma,88,3),0)
    %end
    if parX.refineSoftMask
        M=buildFilter(2*size(x),'tukeyIso',0.125,gpu,1,1);%0.125 (here) - 5 (line below) converges quickly but gives zeros in artifacted slices
        Ma=abs(filtering(abs(x)*mNorm,M,1))/2;
        Ma=1./(abs(Ma).^2+0.01);%Inverse from 0.01 to 400-Direct from 0.0025 to 100
        
        %visReconstruction(dynInd(Ma,88,3))
        Precond=(sum(real(conj(S).*S),4)+Ma).^(-1);
    end
    %visReconstruction(Ma)

    %Solve for T    
    if parX.correct>0 && estT
        %if debug==2;tstart=tic;end
        BlSz=1;
        [T,w,flagw,outlD]=solveTMS2D(x,y,T,H,S,Ak,xkGrid,kkGrid,kGrid,nT,w,flagw,NzSlab,parX.outlP,parX.thplc,BlSz,winic);
        %if debug==2;telapsed=toc(tstart);fprintf('Time solving T: %.4fs\n',telapsed);end                               
    end    
    xant=x-xant;
    xant=real(xant.*conj(xant));
    xant=max(xant(:));       
    if debug>=2;fprintf('Iteration XT %04d - Error %0.2g \n',n,xant);end
    if xant<toler
        if debug>0;fprintf('Iteration XT %04d - Error %0.2g \n',n,xant);end 
        break
    end    
end
x=x*mNorm;
x=extractROI(x,ROI,0,2);
