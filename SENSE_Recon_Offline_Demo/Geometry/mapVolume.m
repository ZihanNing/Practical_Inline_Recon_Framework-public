function [xo,MTy]=mapVolume(xi,y,MTx,MTy,addGrid,extGrid,fillVal,intType)

%MAPVOLUME   Maps the volume x into the geometry of volume y
%   [XO,MTY]=MAPVOLUME(XI,Y,MTX,MTY,{ADDGRID},{EXTGRID},{FILLVAL},{INTTYPE})
%   * XI is the volume to be mapped
%   * Y is the volume onto which to map
%   * MTX are the homogeneous coordinates of the first volume
%   * MTY are the homogeneous coordinates of the second volume
%   * {ADDGRID} serves to add a value to the grid (1 in ReconFrame format 
%   versus 0, default in nifti)
%   * {EXTGRID} extends the ROI of the second volume symmetrically by a given
%   number of voxels
%   * {FILLVAL} is the value to fill the background
%   * {INTTYPE} is the type of interpolation, defaults to 'linear'
%   ** XO is the mapped volume
%   ** MTY are the homogeneous coordinates of the mapping
%

if nargin<5 || isempty(addGrid);addGrid=0;end
if size(addGrid,2)==1;addGrid=repmat(addGrid,[1 3]);end
if size(addGrid,1)==1;addGrid=repmat(addGrid,[2 1]);end

if nargin<6 || isempty(extGrid);extGrid=0;end
if length(extGrid)==1;extGrid=repmat(extGrid,[1 3]);end
if nargin<7;fillVal=[];end
if nargin<8 || isempty(intType);intType='linear';end

gpu=isa(xi,'gpuArray');

MTT=MTx\MTy;%EQUIVALENT TO MTT
MTy=modifyGeometryROI(MTy,-extGrid(:)+1);  
Nsource=size(xi);Nsource(end+1:3)=1;Nsource=Nsource(1:3);
Ndestin=size(y);Ndestin(end+1:3)=1;Ndestin=Ndestin(1:3);y=[];clear y;

rdGrid=generateGrid(Ndestin+2*extGrid,gpu,Ndestin+2*extGrid,ones(1,3)-addGrid(2,:)+extGrid);
[dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rdGrid{1}(:),rdGrid{2}(:),rdGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;rdGrid=[];clear rdGrid;
destinGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');dGrid{4}=[];
if gpu;MTT=gpuArray(MTT);end

sdGrid=generateGrid(Nsource,gpu,Nsource,ones(1,3)-addGrid(1,:));
[sGrid{1},sGrid{2},sGrid{3}]=ndgrid(sdGrid{1}(:),sdGrid{2}(:),sdGrid{3}(:));sdGrid=[];clear sdGrid;

destinGrid=MTT*destinGrid;
for m=1:3
    dGrid{m}=reshape(destinGrid(m,:),Ndestin+2*extGrid);
    if isempty(fillVal)
        dGrid{m}(dGrid{m}<min(sGrid{m}(:)))=min(sGrid{m}(:));
        dGrid{m}(dGrid{m}>max(sGrid{m}(:)))=max(sGrid{m}(:));
    end
end;destinGrid=[];clear destinGrid;
if isempty(fillVal);fillVal=0;end

NXO=size(xi);NXO(1:3)=Ndestin+2*extGrid;NXO(end+1:4)=1;
[xi,NXI]=resSub(xi,4:max(numDims(xi),4));NXI(end+1:4)=1;
xo=zeros([NXO(1:3) NXI(4)],'like',xi);
if ~strcmp(intType,'linear') && ~strcmp(intType,'nearest');xo=gather(xo);xi=gather(xi);sGrid=cellFun(sGrid,@gather);dGrid=cellFun(dGrid,@gather);end        
for n=1:NXI(4);xo=dynInd(xo,n,4,interpn(sGrid{1},sGrid{2},sGrid{3},dynInd(xi,n,4),dGrid{1},dGrid{2},dGrid{3},intType,fillVal));end
if gpu;xo=gpuArray(xo);end
xo=reshape(xo,NXO);
