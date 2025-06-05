function [x,amse]=waveletSVShrinkage(x,dimsM,parR,L)

%WAVELETSVSHRINKAGE  performs SV shrinkage in a wavelet domain, following 
%[1] W Leeb, "Matrix denoising for weighted loss functions and 
%heterogeneous signals", arXiv, 1902.0947v1, 2019
%   [X,AMSE]=WAVELETSVSHRINKAGE(X,DIMSM,{PARR})
%   * X is the array on which to perform the patch-based SVD shrinkage
%   * DIMSM is the last dimension of X along which to build the rows of the
%   matrices
%   * {PARR} are the parameters to build the matrices from the pathes
%   * {L} is a structure containing the noise formation
%   ** XH is the recovered array
%   ** AMSEH is an estimate of the asymptotic mean square error
%

%DEFAULT VALUES AND DIMENSIONS
gpu=isa(x,'gpuArray');
onM=ones(1,dimsM);
if nargin<3 || isempty(parR)
    parR.Verbosity=0;%To plot debugging information
    parR.J=zeros(1,dimsM);%Number of levels for wavelet decomposition
    parR.S=ones(1,dimsM);%Number of spatial subdivisions
    parR.O=ones(1,dimsM);%Overlap factor
    parR.UsePacket=0;%To use wavelet packets
    parR.backNoise=1e-3;%Background noise for wavelet padding
    parR.wname='db4';%Wavelet filter
end
if nargin<4;L=[];end
if length(parR.J)==1;parR.J=repmat(parR.J,[1 dimsM]);end
if length(parR.S)==1;parR.S=repmat(parR.S,[1 dimsM]);end
if length(parR.O)==1;parR.O=repmat(parR.O,[1 dimsM]);end

%WPAD THE DATA FOR WAVELET DECOMPOSITION
NX=size(x);NX(end+1:dimsM)=1;ND=length(NX);
NXM=NX(1:dimsM);
scJ=2.^parR.J;%Scale corresponding to J
scJS=scJ.*parR.S.*parR.O;%Spatial and spectral subdivisions
NXMPad=scJS.*ceil(NXM./scJS);
Ma=ones(NXM,'like',real(x));
Ma=resampling(Ma,NXMPad,2);
if isempty(L)
    L{1}=Ma+parR.backNoise*(1-Ma);
else%We assume no correlations at the moment
    L{1}=reshape(L{1},NXM);
    L{1}=resampling(L{1},NXMPad,2);
    Ma=Ma.*single(L{1}>0);
    L{1}=L{1}+parR.backNoise*(1-Ma);
end

%ADD NOISE
x=resampling(x,NXMPad,2);
x=bsxfun(@times,x,Ma)+parR.backNoise*bsxfun(@times,plugNoise(x,1),1-Ma);

%OVERLAP 
NOverlaps=prod(parR.O);
xhat=x;xhat(:)=0;
what=L{1};what(:)=0;%Weights

%INFO FOR PARCELLATION
NN=prod(NX(dimsM+1:ND));
NXsplit=ones(1,2*dimsM);
NXW=NXMPad./parR.S;
NXsplit(1:2:2*dimsM)=NXW;
NXsplit(2:2:2*dimsM)=parR.S;
permX=1:2*dimsM+1;
permX(dimsM+2:2*dimsM+1)=2:2:2*dimsM;
permX(1:dimsM+1)=[1:2:2*dimsM 2*dimsM+1];

%LOCALIZE IN THE B-SPACE
if isfield(parR,'bvalDec');IJ=parR.bvalDec;else IJ=[];end

%BUILD THE WAVELET MATRICES        
[W,~,I]=buildWaveletMatrix(NXW,parR.wname,parR.J,parR.UsePacket,gpu);
for n=1:length(W)  
    perm=1:dimsM;perm([1 n])=[n 1];
    W{n}=permute(W{n}(:),perm);
end
II=zeros(length(I(:)),max(I(:)));
II=indDim(II,I(:),2,1);
if all(parR.J==0);II=[];W=[];end

for o=1:NOverlaps
    oV=ind2subV(parR.O,o)-1;   
    oV=oV.*NXW./parR.O;
    [xsh,Lsh{1},wsh]=parUnaFun({x,L{1},what},@circshift,oV);

    %PARCELLATION OF SPACE
    xsh=reshape(xsh,[NXsplit numel(x)/prod(NXsplit)]);
    [Lsh{1},wsh]=parUnaFun({Lsh{1},wsh},@reshape,[NXsplit 1]);

    [xsh,Lsh{1},wsh]=parUnaFun({xsh,Lsh{1},wsh},@permute,permX);
    [xsh,Lsh{1},wsh]=parUnaFun({xsh,Lsh{1},wsh},@resSub,dimsM+2:2*dimsM+1);

    NPatches=size(xsh,dimsM+2);
    bP=1;
    for p=1:bP:NPatches;vP=p:min(p+bP-1,NPatches);NP=length(vP);
        %EXTRACT PATCH
        [xp,Lp{1}]=parUnaFun({xsh,Lsh{1}},@dynInd,vP,dimsM+2);

        %ARRANGE COIL INFO
        Lp{1}=reshape(Lp{1}(:),[ones(1,dimsM) prod(NXW) NP]);

        %ARRANGE IN MATRIX FORM
        xp=reshape(xp,[prod(NXW) NN NP]);%Arrange in matrix form

        %CALL THE SOLVER---EXTEND TO 3D ENTRIES!
        [xp,amsep]=heterogeneousSVShrinkage(xp,{Lp,[]},{II,IJ},{W,[]},{NXW,[]},1);
        if strcmp(parR.combRule,'Invva');wamsep=1./(amsep+eps);else wamsep=1;end       
        
        %RESHAPE BACK
        xp=reshape(xp,[NXW NN NP]);
        xsh=dynInd(xsh,vP,dimsM+2,bsxfun(@times,xp,wamsep));
        wsh=dynInd(wsh,vP,dimsM+2,repmat(wamsep,[NXW 1 1]));
    end
    xsh=reshape(xsh,[NXW NN parR.S]);
    wsh=reshape(wsh,[NXW 1 parR.S]);
    [xsh,wsh]=parUnaFun({xsh,wsh},@ipermute,permX);
    xsh=reshape(xsh,[NXMPad NX(dimsM+1:ND)]);
    wsh=reshape(wsh,[NXMPad 1]);
    xhat=xhat+circshift(xsh,-oV);
    what=what+circshift(wsh,-oV);
end
x=bsxfun(@rdivide,xhat,what);

%CONSTRAIN
x=bsxfun(@times,x,Ma);
x=resampling(x,NXM,2);
