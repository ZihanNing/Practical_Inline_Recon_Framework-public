function [y,MT]=extractROI(x,ROI,ext,dims,MT)

%EXTRACTROI extracts a given ROI from the image or fills the full image
%   [X,MT]=extractROI(X,ROI,EXT,{DIMS},{ROI})
%   * X is an input image
%   * ROI is the ROI to be extracted
%   * {EXT} is a flag that determines whether to extract (1) or to fill 
%   (0). It defaults to extract
%   * {DIMS} indicates the active dimensions for ROI extraction, it 
%   defaults to all
%   * {MT} are the homogeneous coordinates of the volume
%   ** Y is the output image
%   ** MT are the updated homogeneous coordinates of the volume
%

ND=min(size(ROI,1),numDims(x));
if nargin<3 || isempty(ext);ext=1;end
if nargin<4 || isempty(dims);dims=1:ND;end
if nargin<5;MT=[];end

ROIDyn=cell(1,length(dims));
for n=1:length(ROIDyn);ROIDyn{n}=ROI(dims(n),1):ROI(dims(n),2);end

if ext
    y=dynInd(x,ROIDyn,dims);
    MT=modifyGeometryROI(MT,ROI,dims);
else
    N=size(x);N(dims)=ROI(dims,3);    
    y=single(zeros(N));
    if isa(x,'gpuArray');y=gpuArray(y);end
    y=dynInd(y,ROIDyn,dims,x);
    MT=modifyGeometryROI(MT,-ROI,dims);%LIKELY TO BE WRONG! BUG!
end
