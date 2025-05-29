clc
cd ('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Projects\GeomTesting')
addpath(genpath('.'));
addpath(genpath('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Software\DISORDER\DefinitiveImplementationRelease07'));
addpath(genpath('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Software\Utilities'));
addpath(genpath('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Projects\Reconstruction'));

dataDir = 'C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Experiments\Data\2021_03_02_GeomTesting_Phantom\';

%% Load data for specified acquisition
scan = 1;
if scan ==1 
    fileName = strcat( dataDir ,'ISD-PACS\202103021142_AWP79013_4_GRE_MIDRES_DISO_16_16_LIN_LIP\s004a1001.nii');%Default
    fileNameRaw = strcat( dataDir ,'20210302_FID18142_P_GRE_MIDRES_DISO_16_16_LIN_LIP.dat');
    dimFlip = [2];
elseif scan ==2
    fileName = strcat( dataDir ,'ISD-PACS\202103021142_AWP79013_20_GRE_DISO_PARot10v2/s020a1001.nii');%PA10. Attention: The actual rotation is AP, not PA
    fileNameRaw = strcat( dataDir ,'20210302_FID18151_P_GRE_DISO_PARot10v2.dat');
    dimFlip = [3];
elseif scan ==3
    fileName = strcat( dataDir ,'ISD-PACS\202103021142_AWP79013_16_GRE_DISO_HFRot10v2/s016a1001.nii');%HF10
    fileNameRaw = strcat( dataDir ,'20210302_FID18149_P_GRE_DISO_HFRot10v2.dat');
    dimFlip = [3];
elseif scan ==4
    fileName = strcat( dataDir ,'ISD-PACS\202103021142_AWP79013_18_GRE_DISO_LRRot10v2/s018a1001.nii');%LR10 
    fileNameRaw = strcat( dataDir ,'20210302_FID18150_P_GRE_DISO_LRRot10v2.dat');
    dimFlip = [2];
end

%% Image and orientation from the NIFTI header
%Load image
nii = load_untouch_nii(fileName);
x = single(nii.img); %PE-RO-SL with PE=AP / RO = IS / SL=LR ~= Right handed!

%%% sform
sForm = diag(ones(1,4));
sForm(1,:)=nii.hdr.hist.srow_x;
sForm(2,:)=nii.hdr.hist.srow_y;
sForm(3,:)=nii.hdr.hist.srow_z;

%% Orientation from raw data
%Load raw data and try to 
TW=mapVBVD(fileNameRaw);
rec=TWIX2rec(TW);
% Attention: Siemens raw k-space is converted to image space by taking the
% FORWARD Fourier transform. This is elaborated on in a post on the IDEA
% forum (see FFTscale.png in the Documentation folder)

%% Read image and re-order dimensions
xRaw = sqrt( mean(rec.y.^2, 4)); %RO-PE-SL since taken from TWIX data 
MSRaw{1} = rec.Enc.AcqVoxelSize; %Voxel spacing

%Re-order so PE-RO-SL (Siemens convention: see IDEA manual)
xRaw = permute(xRaw,[2 1 3 4]);
MSRaw{1} = permute(MSRaw{1},[2 1 3 4]);%need to permute spacing as well

%% Inspect raw data image
plotND(abs(xRaw));

%% Orientation information from TWIX header
%%% Quaternions + translation
quaternionRaw = TW.image.slicePos(4:7,1);
tranRaw = TW.image.slicePos(1:3,1);%mm I guess/wr.t. centre of FOV

%%% Construct transformation matrix from PE-RO-SL (PRS) to Patient Coordinate System (PCS)
T_PRS2PCS = diag(ones(1,4));
T_PRS2PCS(1:3,1:3) = quaternion_rotation_matrix((quaternionRaw))*diag(MSRaw{1}); 
T_PRS2PCS(1:3,4) = tranRaw;

or = ( ((size(xRaw)+0)/2 - 0*[1 1 1]) +[0 0 1*(-.5)])' ; % account for fact that tranRaw is referred to the center of FOV, not the first element in the array
T_PRS2PCS(1:3,4) = T_PRS2PCS(1:3,4) - ( T_PRS2PCS(1:3,1:3)*or +[0;0;0]);

%%% Convert from PCS to Right Anterior Superior (convention in NIFTI header)
PCS2RAS = diag([-1 -1 1 1]);%For HeadFirstSupine. Can be generalised from TW.hdr.Dicom.tPatientPosition  
T_PRS2RAS = PCS2RAS * T_PRS2PCS;% From PE-RO-SL to RAS

y =[]; y{1} = xRaw;
MTRaw = []; MTRaw{1} = T_PRS2RAS;
writeNII(fileNameRaw(1:end-4), {'Raw'},y, MSRaw, MTRaw)

%% Check if orientation matrix is the same when extracted from NIFTI header and TWIX header
%Inspect the axis to flip between NIFTI array and xRaw: can only compare
%the orientation matrix if array has dimensions stored in the same way
dimFlip = [3]; %This changes over scans ... I don't know why but the important thing is that the raw data itself is consistent with the orientation
xPlot = xRaw; 
MTPlot = MTRaw{1};

for i=1:length(dimFlip); xPlot = flip(xPlot, dimFlip(i));end
for i=1:length(dimFlip); MTPlot(:,dimFlip(i)) = -MTPlot(:,dimFlip(i)) ;end%s_form modification (sign difference)
for i=1:length(dimFlip); MTPlot(:,4) = MTPlot(:,4) - (size(xPlot, dimFlip(i))-1) .* MTPlot(:,dimFlip(i));end
plot_(abs(x),[],[52 40 50],abs(xPlot));

%%% Compare
sForm
MTPlot

diff = MTPlot(:,4)-sForm(:,4);
MTPlot\diff

% One can see the orientation matches perfectly, but there is a slight
% offset in the translation. To better understand, we would have to look at
% what the origina of all reference frames are (e.g. including table
% position, etc.).
%One parameter that might be helpful later: TW.hdr.Dicom.lSBCSOriginPositionZ

%@Lucilio update: This offset is I think how ISDPACS deals with this
%convention. I think if you would look at the arrays, they would also be
%shifted one unit.

%% Orientation of volume? Consistent with name given
% Here we calculate the Euler Angles (EA) of the transformation matrix and
% this can be used to assess whether the FOV is tiled as provided in the
% filename

[U,~,V] = svd(T_PRS2RAS(1:3,1:3));%Only rotation
temp = U*V';%take out voxel dimensions - necessary for rotation matrix
fileNameRaw
EA = 180/pi*rotm2eul(temp,'ZYX')% since active convention / ZYX is SAR;

