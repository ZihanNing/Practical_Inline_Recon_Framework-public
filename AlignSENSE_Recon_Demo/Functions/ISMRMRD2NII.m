
function [x, MS, MT] = ISMRMRD2NII (fileName, writeNIIFlag, disph5Info, plotFlag)

%ISMRM2NII Converts ISMRM iamge data (ISMRM-RD) to a NIFTI structure.
%   * FILENAME is the filename of the raw data file.
%   * {WRITENIIFLAG} a flag whether to write the NIFTI file in .nii extension.
%   * {DISPH5INFO} a flag whether to display the .h5 file information.
%   * {PLOTFLAG} a flag whether to display the image.
%   ** x is the data array
%   ** MS is the resolution
%   ** MT is the NIFTI orientation information
%

if nargin < 2 || isempty(writeNIIFlag);writeNIIFlag=0;end
if nargin < 3 || isempty(disph5Info);disph5Info=0;end
if nargin < 4 || isempty(plotFlag);plotFlag=0;end

%% Control
if ~exist(strcat(fileName,'.h5'), 'file');error(['File ' strcat(fileName,'.h5') ' does not exist. Please generate it using the ISMRMRD converter.']);end

%% Read
h5Info = h5info(strcat(fileName,'.h5'));
if disph5Info; h5disp(strcat(fileName,'.h5'));end

groups = h5Info.Groups(:);
numSeries = length(groups.Groups);
Names = groups(:).Name;

nameStudy = Names;

if numSeries>1
    tt = strsplit(groups.Groups(1).Name,'/');
else
    tt = strsplit(groups.Groups.Name,'/');
end
nameImage = tt{end};%'image_0';

data = h5read(strcat(fileName,'.h5'), sprintf('%s/%s/data',nameStudy,nameImage));
hdr = h5read(strcat(fileName,'.h5'), sprintf('%s/%s/header',nameStudy,nameImage));
attributes = h5read(strcat(fileName,'.h5'), sprintf('%s/%s/attributes',nameStudy,nameImage));

%% GEOMETRY COMPUTATION
%SCALING
rec.Enc.AcqVoxelSize = single( hdr.field_of_view') ./ single(hdr.matrix_size');
rec.Par.Mine.Asca = diag([ rec.Enc.AcqVoxelSize 1]);

%ROTATION
rec.Par.Mine.Arot = single(eye(4));
rec.Par.Mine.Arot(1:3,1:3)=single( cat(2, hdr.read_dir(:,1), hdr.phase_dir(:,1), hdr.slice_dir(:,1))); %RO-PE-SL to PCS

%TRANSLATION
rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = single(hdr.position(1:3,1));%RO-PE-SL to PCS

% account for fact that  meas.head.position is referred to the center of FOV, not the first element in the array
if isstruct(data)
    dataBackup = data;
    if isfield(dataBackup,'real');data = dataBackup.real;end
    if isfield(dataBackup,'imag');data = data + 1i* dataBackup.imag;end
    dataBackup = [];
end
NY = size(data); NY=NY(1:3);
orig = ((NY+0)/2 - 0 - [0 0 .5])'; orig=orig(1:3);

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)- rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;
                    
%COMBINED MATRIX
rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

%MOVE TO RAS
rec.Par.Mine.PCS2RAS = diag([-1 -1 1 1]);%For HeadFirstSupine. Can be generalised from hdr.measurementInformation.patientPosition
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From RO-PE-SL to RAS
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;

%% ASSIGN
x = data;
MS = rec.Enc.AcqVoxelSize;
MT = rec.Par.Mine.APhiRec;

%% WRITE
if writeNIIFlag
    xWrite{1}=x;
    MSWrite{1}=MS;
    MTWrite{1}=MT;
    addpath(genpath('/home/ybr19/Software/Utilities'));
    addpath(genpath('/home/ybr19/Software/DISORDER'));
    writeNIIExt(fileName, {''},xWrite,MSWrite,MTWrite);
end

%% PLOT
if plotFlag
    addpath(genpath('/home/ybr19/Software/Utilities'));
    addpath(genpath('/home/ybr19/Software/DISORDER'));
    plotND([], x, [],[], 0, {[],2},MT);
end

%%% FROM https://github.com/aTrotier/EDUC_GT_MATLAB/wiki/Demo-1-:-Bucket
% function img = read_image_h5(filename)
% if nargin < 1
%     [file, PATHNAME] = uigetfile('*.h5');
%     filename = fullfile(PATHNAME,file);
% end
% 
% S=hdf5info(filename);
% 
% if exist(filename, 'file')
%     dset = ismrmrd.Dataset(filename, 'dataset');
% else
%     error(['File ' filename ' does not exist.  Please generate it.'])
% end
% 
% % hdr = ismrmrd.xml.deserialize(dset.readxml);
% 
% disp(dset.fid.identifier)
% 
% S=hdf5info(filename);
% attributes=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(1).Name;
% dataset=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(2).Name;
% header=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(3).Name;
% 
% img=hdf5read(filename,dataset);
% img=squeeze(img);
% end