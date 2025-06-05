
function x=morphFourier(x,r,d,mirr,dist,sep)

%MORPHFOURIER   Applies morphological operations with ellipsoidal
%structuring elements in Fourier space
%   X=MORPHFOURIER(X,R,D,{MIRR},{DIST},{SEP})
%   * X is the mask or image for which to apply a given morphological
%   operation
%   * R are the radii of the structuring element along different
%   orientations (if positive dilation, if negative erosion)
%   * D is the spacing of the image
%   * {MIRR} are the boundary conditions. Defaults to 1, Neumann (0 for
%   periodic)
%   * {DIST} indicates whether to compute a distance function to the mask,
%   defaults to 0
%   * {SEP} indicates whether to use a separable distance, defaults to 0
%   * X is a transformed mask or image
%

if ~exist('dist','var') || isempty(dist);dist=0;end;
if ~exist('sep','var') || isempty(sep);sep=0;end;

gpu=isa(x,'gpuArray');

NW=length(d);
if ~exist('mirr','var') || isempty(mirr);mirr=ones(1,NW);end;

rd=abs(r)./d;rd(end+1:3)=1;%rd in voxel units
M=size(x);M(end+1:NW)=1;M=M(1:NW).*(1+mirr);%Mirrir for Neumann condition
rGrid=generateGrid(M,gpu,M,ceil((M+1)/2));
for n=length(rGrid)+1:3;rGrid{n}=0;end
%This was too gpu memory demanding
%[sGrid{1},sGrid{2},sGrid{3}]=ndgrid(rGrid{1}(:),rGrid{2}(:),rGrid{3}(:));
%sGrid=abs(bsxfun(@times,cat(4,sGrid{:}),1./permute(rd,[1 3 4 2])));
%if ~sep;idx=find(sum(sGrid.^2,4)<=1);else idx=find(dynInd(sGrid,1,4)<=1 & dynInd(sGrid,2,4)<=1 & dynInd(sGrid,3,4)<=1);end;sGrid=[];
if ~sep
    for n=1:3;rGrid{n}=(rGrid{n}/(rd(n)+eps)).^2;end%Grid normalised to rd
    idx=find(bsxfun(@plus,bsxfun(@plus,rGrid{1},rGrid{2}),rGrid{3})<=1);rGrid=[];%So here thresholding w.r.t. 1
else%This may still produce problems with gpu memory
    [sGrid{1},sGrid{2},sGrid{3}]=ndgrid(rGrid{1}(:),rGrid{2}(:),rGrid{3}(:));
    for n=1:3;sGrid{n}=(sGrid{n}/(rd(n)+eps)).^2;end
    idx=find(sGrid{1}<=1 & sGrid{2}<=1 & sGrid{3}<=1);rGrid=[];sGrid=[];
end
M(end+1:3)=1;
se=zeros(M,'like',real(x));
se(idx)=1;
NP=length(idx);idx=[];
for m=1:NW
    se=ifftshift(se,m);
    se=fftGPU(se,m);%Bring kernel to Fourier domain
end
%Have to extract half the kernel as constructed in Cosine domain. Analogous with using builFilter(2*N, ... , c=1), which results in a filter of size N.
v=cell(1,NW);
for m=1:NW
    if ~mirr(m);v{m}=1:M(m);else v{m}=1:ceil(M(m)/2);end
end
se=dynInd(se,v,1:NW);

x=filtering(x,se,1);

if ~dist
    if r<0;x=single(abs(x)>(NP-0.5));else x=single(abs(x)>0.5);end%kernel chosesn so that PSF has certain width
else
    x=(real(x)/NP).^dist;
    x=mapToZeroOne(x);
end