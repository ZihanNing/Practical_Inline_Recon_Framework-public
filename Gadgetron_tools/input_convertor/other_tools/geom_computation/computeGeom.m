function geom = computeGeom(AcqVoxelSize,matrixSize,hdr,image_hdr,TablePosTra,computeGeom_ACS,ACSSize)
%%% This function is to compute geom based on the twix-like converted raw
%%%
%%% INPUT:
%%% AcqVoxelSize: the voxel size of the highres data
%%% matrixSize: this parameter should be align with the AcqVoxelSize
%%% hdr: twix.hdr
%%% image_hdr: twix.image
%%% TablePostTra: twix.hdr.Dicom.lGlobalTablePosTra
%%% computeGeom_ACS: 1- compute the geom of calibration based on the geom
%%% of high-res data
%%% ACSSize: the matrix size of calibration scan
%%% 
%%% OUTPUT:
%%% geom - for disorder function: rec.Par.Mine = geom;
%%%
%%% by Zihan @ king's
%%% March-2025

if nargin < 5 || isempty(TablePosTra); TablePosTra=[];end 
if nargin < 6 || isempty(computeGeom_ACS); computeGeom_ACS=0;end 
if (nargin < 7 || isempty(ACSSize)) && (computeGeom_ACS==1); computeGeom_ACS=0;end 

fprintf('Computing geometry for reconstruction array.\n');

%%% CHANGE DIMENSION TO [LIN COL PAR]
% readin as [Col Lin Par]
ND = 6; perm=1:ND; perm(1:3)=[2 1 3];
AcqVoxelSize = AcqVoxelSize(perm(1:3));
matrixSize = matrixSize(perm(1:3));
if computeGeom_ACS; ACSSize = ACSSize(perm(1:3)); end

%%% COMEPUTE GEOM
idxToUse = [ image_hdr.sort_Lin(1,1) , image_hdr.sort_Par(1) ];%Make sure we don't extract information from outside the elliptical shutter;
% might should use Lin & Par before zero-padding, need test later

% scaling
geom.Asca=diag([AcqVoxelSize 1]);

% rotation
geom.Arot = double(eye(4));
geom.Arot(1:3,1:3)=double(cat(2, image_hdr.read_dir(:,1), image_hdr.phase_dir(:,1), image_hdr.slice_dir(:,1))); % ZN: input bucket %RO-PE-SL to PCS 
Arot_tmp = geom.Arot;
geom.Arot(1:3,1) = Arot_tmp(1:3,2);geom.Arot(1:3,2) = Arot_tmp(1:3,1); 
clear Arot_tmp
    
% translation
geom.Atra=eye(4); 
geom.Atra(1:3,4) = double(image_hdr.slicePos(:,1)); 
if ~isempty(TablePosTra)
    geom.Atra(3,4) = geom.Atra(3,4) + TablePosTra; % considered the position of table (in HF direction)
end
% account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = matrixSize;%Set inf to make sure this element is replaced
orig = ( ceil((N(1:3)+1)/2) - [0 0 .5] )';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK

geom.Atra(1:3,4)= geom.Atra(1:3,4)- geom.Arot(1:3,1:3)*geom.Asca(1:3,1:3)*orig;

%COMBINED MATRIX
geom.MTT=eye(4);%YB: not sure what this does --> de-activated
geom.APhiRec=geom.MTT*geom.Atra*geom.Arot*geom.Asca;

%MOVE TO RAS
geom.patientPosition = hdr.Dicom.tPatientPosition;
geom.PCS2RAS = getPCS2RAS(geom.patientPosition);
geom.APhiRec = geom.PCS2RAS  * geom.APhiRec;% From [Lin Col Par] to RAS
[~,geom.APhiRec] = mapNIIGeom([], geom.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED
geom.APhiRecOrig = geom.APhiRec;%Backup as how was read in

%DEDUCE ACQUISITION ORDER
[geom.MPS, geom.slicePlane] = acquisitionOrder(geom.APhiRec);
geom.RO = geom.MPS(4:5);
geom.PE1 = geom.MPS(1:2);
geom.PE2 = geom.MPS(7:8); 
fprintf('(Slice orientation: %s\n', geom.slicePlane);
fprintf('    Readout: %s\n', geom.RO );
fprintf('    1st Phase encode direction: %s\n', geom.PE1 );
fprintf('    2nd Phase encode (slice) direction: %s\n\n', geom.PE2 );

geom.FoldOverDir = geom.MPS(1:2);%First PE direction
geom.FatShiftDir = geom.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

%%% GEOMETRY COMPUTATION FOR AUTOCALIBRATION ARRAY
if computeGeom_ACS && exist('ACSSize','var')
    if ~isempty(ACSSize) % ZN: to enable ACS line
        fprintf('Computing geometry for autocalibration array.\n');
        YSize = matrixSize;%Set inf to make sure this element is replaced
        [geom.ACSVoxelSize, geom.APhiACS] = mapNIIGeom([], geom.APhiRec,'resampling',[],YSize,ACSSize);%TO CHECK:AD HOC - REVERSE ENGINEERED
    end
end

%%% CONVERT BACK TO [COL LIN PAR]
perm=1:ND; perm(1:3) = [2 1 3];
[AcqVoxelSize , geom.APhiRec] = mapNIIGeom(AcqVoxelSize , geom.APhiRec,'permute', perm);
if isfield(geom,'APhiACS') && isfield(geom,'ACSVoxelSize')
    [geom.ACSVoxelSize, geom.APhiACS] = mapNIIGeom(geom.ACSVoxelSize, geom.APhiACS,'permute', perm);
end 

geom.permuteHist = [];geom.permuteHist{1} = perm(1:4);

%%% Make MPS consisitent with RO-PE-SL
MPStemp = geom.MPS ;
geom.MPS(1:2) = MPStemp(4:5);
geom.MPS(4:5) = MPStemp(1:2);