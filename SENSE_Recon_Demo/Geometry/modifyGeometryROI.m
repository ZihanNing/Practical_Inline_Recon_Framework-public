function MT=modifyGeometryROI(MT,ROI,dims)

%MODIFYGEOMETRYROI   Modifies the geometry using a given ROI
%   MT=MODIFYGEOMETRYROI(MT,ROI,{DIMS})
%   * MT are the homogeneous coordinates of the volume
%   * ROI is a given ROI to be extracted
%   * {DIMS} indicates the active dimensions for ROI extraction, it 
%   defaults to all (maximum of 3)
%   ** MT are the updated homogeneous coordinates of the volume
%

if isempty(MT);return;end

ND=min(size(ROI,1),3);
if nargin<3 || isempty(dims);dims=1:ND;end

vROI=zeros(4,1);vROI(4)=1;                  
vROI(dims)=vROI(dims)+ROI(dims,1)-1;
vROI=MT*vROI;
MT(1:3,4)=vROI(1:3); 
