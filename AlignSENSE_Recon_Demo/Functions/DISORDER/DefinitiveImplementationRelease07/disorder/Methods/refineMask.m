
function x=refineMask(x,parS,voxsiz)

%REFINEMASK   Refines or generates a mask
%   X=REFINEMASK(X,{PARS},{VOXSIZ})
%   * X is the mask or image for which a mask is extracted
%   * {PARS} is a structure with the parameters to refine the mask
%   * {VOXSIZ} is the spacing of the image space
%   ** X is a mask
%

if nargin<2 || isempty(parS)
   parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
   parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
   parS.Otsu=[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)
   parS.nErode=0;%Erosion for masking (in mm)
   parS.nDilate=0;%Dilation for masking (in mm)
   parS.conComp=1;%Whether to get the largest connected component after erosion (when >0). If set to 2, holes in masks will be filled, provided they are not the FOV boundaries. 
else
    if ~isfield(parS,'maskNorm');parS.maskNorm=1;end
    if ~isfield(parS,'maskTh');parS.maskTh=1;end
    if ~isfield(parS,'Otsu');parS.Otsu=[0 1];end
    if ~isfield(parS,'nErode');parS.nErode=0;end
    if ~isfield(parS,'nDilate');parS.nDilate=0;end
    if ~isfield(parS,'conComp');parS.conComp=1;end
end     
if nargin<3 || isempty(voxsiz);voxsiz=ones(1,3);end

gpu=isa(x,'gpuArray');gpuIn=single(gpuDeviceCount && ~blockGPU);

if gpuIn;x=gpuArray(x);end

%FEATURE EXTRACTION FOR MASKING
xabs=abs(x);
xabs=xabs.^parS.maskNorm;
if isempty(parS.Otsu)
    parS.maskTh=parS.maskTh^parS.maskNorm;
    meanx=mean(xabs(:));
    x=single(xabs>meanx*parS.maskTh);
elseif length(unique(xabs(:)))>2    
    xabs=reshape(imquantize(gather(xabs(:)),multithresh(gather(xabs(xabs>1e-6)),length(parS.Otsu)-1)),size(xabs));
    if gpu;xabs=gpuArray(xabs);end        
    x=single(ismember(xabs,find(parS.Otsu)));%YB: because of the find(Otsu), x only consists of 1 value!
else%We already have binary mask
    x=xabs;
end
xabs=[];

NW=min(numDims(x),3);
mirr=ones(1,NW);%For Neumann boundary conditions, the most natural in this application (0 for periodic)
if any(parS.nErode>0)
    if length(parS.nErode)==1;parS.nErode=parS.nErode*ones(1,NW);end
    x=morphFourier(x,-parS.nErode(1:NW),voxsiz(1:NW),mirr);
end
if parS.conComp
      CC=bwconncomp(logical(gather(x)));
      numPixels = cellfun(@numel,CC.PixelIdxList);
      [~,idx] = max(numPixels);
      x(:)=0;x(CC.PixelIdxList{idx})=1;
      x=single(x);
      if gpu;x=gpuArray(x);end
end
if any(parS.nDilate>0)
    if length(parS.nDilate)==1;parS.nDilate=parS.nDilate*ones(1,NW);end
    x=morphFourier(x,parS.nDilate(1:NW),voxsiz(1:NW),mirr);
end

if parS.conComp==2
    x=1-x;%Invert mask
    CC=bwconncomp(logical(gather(x)));
    NN=size(x);NN(end+1:3)=1;%This won't work in 4D probably
    for n=1:length(CC.PixelIdxList)%By inverting x, the binary mask can now contain more than 1 component       
        PixelInd=ind2subV(NN,CC.PixelIdxList{n});              
        PixelIndInv=bsxfun(@minus,NN+1,PixelInd);
        if ~any(PixelInd(:)==1) && ~any(PixelIndInv(:)==1);x(CC.PixelIdxList{n})=0;end%Not array boundary areas, so by setting to zero (and later inerting), you fill mask
    end
    x=1-x;%Invert back
    x=single(x);
    if gpu;x=gpuArray(x);end
end

if ~gpu;x=gather(x);end
