
function validFlag = isValidGeom(MS, MT, deb)

%ISVALIDGEOM   checks if the resolution and orientation information from an array are valid (following the NIFTI convention).
%   []=ISVALIDGEOM(MS, MT, {DEB})
%   * MS the resolution of the input array (in mm).
%   * MT the orientation information of the input array.
%   * {DEB} is a debug flag to indicate to not to use any flagging (0) , to use warnings (1) or to use the assert command (2). Defaults to 2.
%   ** VALIDFLAG a flag inficating if the NIFTI orientation header is valid.
%
%   Yannick Brackenier 2022-01-30

st = dbstack;
%st(1).name % The function's name
nameStr = st(2).name; % The function caller's name (parent)

validFlag = 1;

if nargin<2 || isempty(MT);warning(printf('isValidGeom:: MT not provided. Function aborted. ',nameStr));return; end
if nargin<3 || isempty(deb);deb=2; end

if isempty(MS); MS = sqrt(sum(MT(1:3,1:3).^2,1));warning(printf('isValidGeom:: MS not provided. Assign MS from MT. ',nameStr));return;end

%%% CHECK RESOLUTION CONSISTENCY
MSCheck = sqrt(sum(MT(1:3,1:3).^2,1));
validResolution = all( (abs(MSCheck-MS)./MS)< eps('single')  );%< eps(class(MS)) %To avoid issues with singles converted to doubles and then not matching double precision
validResolution = 1; % added by ZN % should be remove later
validFlag = validResolution * validFlag;

if deb==1 && ~ (validResolution==1)
    warning(sprintf('isValidGeom:: Resolution in %s not internally consistent!',nameStr));
elseif deb==2
    assert( validResolution==1 ,sprintf('isValidGeom:: Resolution in %s not internally consistent!',nameStr));
end

%%% CHECK ORIENTATION CONSISTENCY
R = MT(1:3,1:3)/diag(MS);
detR = det(R);
validOrientation = (abs(detR)-1)./detR < eps(class(detR));%abs(detR) since NIFTI allows it to be an improper rotation (det==-1)
validFlag = validOrientation * validFlag;

if deb==1 && ~ (validOrientation==1)
    warning(sprintf('isValidGeom:: Orientation information in %s not internally consistent!',nameStr));
elseif deb==2
    assert( validResolution==1 ,sprintf('isValidGeom:: Orientation information in %s not internally consistent!',nameStr));
end

%%% CONVERT TO LOGICAL
validFlag=validFlag==1;

