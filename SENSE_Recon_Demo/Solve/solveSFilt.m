function [S,M]=solveSFilt(x,y,voxsiz,maskTh)

%SOLVERSFILT   Solves for sensitivities using a median filter of the ratio
%of the surface and body coils
%   [S,M]=SOLVERSFILT(X,Y,{VOXSIZ},{MASKTH}) 
%   * X is the body coil
%   * Y is the surface coil
%   * {VOXSIZ} is the voxel size
%   * {MASKTH} is a threshold of body coil intensities for mask extraction
%   ** S are the estimated sensitivities
%   ** M is a signal covariance map for regularization
%

if nargin<3 || isempty(voxsiz);voxsiz=ones(1,3);end
if nargin<4 || isempty(maskTh);maskTh=0.5;end

N=size(x);
resol=4;
fact=voxsiz/resol;%We operate at 4mm
Nf=round(N.*fact);
gpu=useGPU;

%0.25 parameter is dependent of resolution, it should be made independent
%H=buildFilter(2*Nf,'gauss',0.5,gpu,0.5,1);
H=buildFilter(2*Nf,'gauss',0.125,gpu,0.5,1);
%H=buildFilter(2*Nf,'gauss',0.5,gpu,0.5,1);
%H=buildFilter(2*Nf,'gauss',0.125,gpu,0.5,1);
%H=buildFilter(2*Nf,'tukeyIso',0.25,gpu,1,1);
%H=buildFilter(2*Nf,'tukeyIso',0.125,gpu,1,1);

%lim=0.75;
%H=buildFilter(2*Nf,'gauss',lim,gpu,1,1);
if gpu;x=gpuArray(x);end
x=resampling(x,Nf,0,2);
S=y;S(:)=0;

%parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
parS.maskTh=maskTh;%0.2;%Threshold of the body coil intensities for mask extraction%It was 0.2
parS.Otsu=[];%Binary vector of components for multilevel extraction (it picks those components with 1)

M=abs(refineMask(x,parS));

%visSegment(x,M,0)
if any(M(:)==0);
    MH=filtering(M,H,1);
    MH=real(MH);
    MH=MH-min(MH(:));
    MH=MH+0.1*prctile(MH(:),75);
    MH=MH/max(MH(:));
else
    MH=M;
end

%visSegment(MH.*1i,[],0)


for s=1:size(y,4)
    ys=dynInd(y,s,4);
    if gpu;ys=gpuArray(ys);end
    ys=resampling(ys,Nf,0,2);
    Ss=bsxfun(@rdivide,ys,x);
    Ss(isnan(Ss))=0;        
    Ss=bsxfun(@times,M,Ss)+bsxfun(@times,cdfFilt(Ss,'med',[3 3 3]),1-M);
    Ss=bsxfun(@times,filtering(Ss,H,1),1./MH);
    Ss=resampling(Ss,N,0,2);
    S=dynInd(S,s,4,gather(Ss));
end
if gpu
    y=gpuArray(y);
    S=gpuArray(S);
end
M=sum(conj(S).*y,4)./(normm(S,[],4)+1e-12);
S=gather(S);
H=resampling(H,N,4,2);
M=gather(100*filtering(M,H,1));