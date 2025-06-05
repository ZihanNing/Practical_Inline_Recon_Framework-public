function xo=mapVolume(xi,y,MTx,MTy,addGrid,extGrid, interpolation)

%MAPVOLUME   Maps the volume x into the geometry of volume y
%   XO=MAPVOLUME(XI,Y,MTX,MTY,{ADDGRID},{EXTGRID},{INTERPOLATION})
%   * XI is the volume to be mapped
%   * Y is the volume onto which to map
%   * MTX are the homogeneous coordinates of the first volume
%   * MTY are the homogeneous coordinates of the second volume
%   * {ADDGRID} serves to add a value to the grid (1 in ReconFrame format versus 0, default in nifti)
%   * {EXTGRID} extends the ROI of the second volume symmetrically by a given number of voxels
%   * {INTERPOLATION} indicates the interpolation used in interpn.m. Defaults to linear interpolation
%   ** XO is the mapped volume
%

if nargin<5 || isempty(addGrid);addGrid=0;end
if length(addGrid)==1;addGrid=repmat(addGrid,[1 3]);end

if nargin<6 || isempty(extGrid);extGrid=0;end
if length(extGrid)==1;extGrid=repmat(extGrid,[1 3]);end

if nargin<7 || isempty(interpolation);interpolation='linear';end%Others ('cubic' or 'spline') can be used but require more memory

gpu=isa(xi,'gpuArray');

MTT=MTx\MTy;%EQUIVALENT TO MTT; YB: MTT=inv(MTx)*MTy
Nsource=size(xi);Nsource(end+1:3)=1;Nsource=Nsource(1:3);
Ndestin=size(y);Ndestin(end+1:3)=1;Ndestin=Ndestin(1:3);

rdGrid=generateGrid(Ndestin+2*extGrid,gpu,Ndestin+2*extGrid,ones(1,3)-addGrid+extGrid);
[dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rdGrid{1}(:),rdGrid{2}(:),rdGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;%YB: dGrid{4} ones since homogeneous coordinate
destinGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');dGrid{4}=[];
if gpu;MTT=gpuArray(MTT);end

sdGrid=generateGrid(Nsource,gpu,Nsource,ones(1,3)-addGrid);
[sGrid{1},sGrid{2},sGrid{3}]=ndgrid(sdGrid{1}(:),sdGrid{2}(:),sdGrid{3}(:));

destinGrid=MTT*destinGrid;%YB:since MTT=inv(MTx)*MTy this goes from array grid y --> RAS --> back to gridx; so this is now the query grid used in interpn 
for m=1:3
    dGrid{m}=reshape(destinGrid(m,:),Ndestin+2*extGrid);
    dGrid{m}(dGrid{m}<min(sGrid{m}(:)))=min(sGrid{m}(:));
    dGrid{m}(dGrid{m}>max(sGrid{m}(:)))=max(sGrid{m}(:));
end

NXO=size(xi);NXO(1:3)=Ndestin+2*extGrid;NXO(end+1:4)=1;%YB: NXO can have in 5th... dimension
[xi,NXI]=resSub(xi,4:max(numDims(xi),4));NXI(end+1:4)=1;%YB: Put every dimension > 4  in 4th
xo=zeros([NXO(1:3) NXI(4)],'like',xi);
for n=1:NXI(4);xo=dynInd(xo,n,4,interpn(sGrid{1},sGrid{2},sGrid{3},dynInd(xi,n,4),dGrid{1},dGrid{2},dGrid{3},interpolation,0));end
xo=reshape(xo,NXO);%YB: restore higher dimensions
