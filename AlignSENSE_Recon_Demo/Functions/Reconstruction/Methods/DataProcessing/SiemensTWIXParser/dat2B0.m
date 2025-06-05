
function [] = dat2B0(fileRef, fileName, supportReadout, writeNIIFlag, writeRAW, pathOu, isPilotTone)

%DAT2B0 processes a Siemens .dat file from an multiple echo time (TE) acquisition and calculated B0 map. 
%   []=DAT2B0(REFNAME,FILENAME,{UPPORREADOUT},{WRITENIIFLAG},{PATHOU},{ISPILOTTONE})
%   * REFNAME is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * FILENAME is the name of the raw data (.dat) of the multiple TE acquisition.
%   * {SUPPORTFOV} is the ralative range of the readout FOV to be extracted.
%   * {WRITENIIFLAG} is a flag to write the B0 map to NIFTI files (1) or to write both maps and raw images (2). Defaults to 1.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the multiple TE file is stored.
%   * {isPilotTone} indicated whether Pilot Tone was used during this acquisition. This activates a median filter across k-space to detect outliers. Poor performance, so only use when applicable.
%       

[pathOuTemp,file]=fileparts(fileName);
fprintf('Parsing TWIX object to reconstruction structure:\n  B0 mapping: %s\n',file);

if nargin<3 || isempty(supportReadout); supportReadout=[];end
if nargin<4 || isempty(writeNIIFlag); writeNIIFlag=2;end
if nargin<5 || isempty(writeRAW); writeRAW=1;end
if nargin<6 || isempty(pathOu); pathOu=[];end
if nargin<7 || isempty(isPilotTone); isPilotTone=0;end
if isempty(pathOu); pathOu = pathOuTemp;end

%%% INVERT TWIX DATA
[recB0 ,TWB0] = dat2Rec(fileName, supportReadout, 0, 1, [], isPilotTone);

%%% ADD PARAMETERS FOR AFI ACQUISITION
recB0.Par.Labels.TE = TWB0.hdr.Meas.alTE/1000;%In ms
recB0.Par.Labels.TE(recB0.Par.Labels.TE ==0)=[];
recB0.Par.Labels.SliceGaps = 0; %Used in solveB07T.m - disabled here
recB0.Alg.parS=[];recB0.Dyn.Typ2Rec=[];
recB0.Alg.parS.conComp=2;%Set to 0 if multiple objects in image
recB0.Alg.parU.UnwrapMeth='CNCG';
recB0.Alg.useTE3=length(recB0.Par.Labels.TE)>2;

%%% LOAD SENSITIVITY MAPS FOR RECONSTRUCTION
if ~exist(strcat(fileRef,'.mat'),'file'); dat2Se(fileRef);end
ss = load(strcat(fileRef,'.mat'));recS =ss.recS; ss=[];
N = size(recB0.y);
MTS = recS.Par.Mine.APhiRec; MTB0 = recB0.Par.Mine.APhiRec;
recB0.S = mapVolume (recS.S, ones(N(1:3)) , MTS, MTB0); recS=[];

%%% RECONSTRUCT IMAGES
x = SWCC(recB0.y, recB0.S);%Sensitivity Weighted Coil Combination
nD = numDims(x);
x = permute(x,[1:3 nD 4:(nD-1) ] );%Different echoes stored in 4th dimension
recB0.x=x;

%x = dynInd( x, [1 2], 4 , dynInd(x, [2 1],4) ); For B0GRE

recB0.Alg.Sensitivities = fileRef;%Store the name of the sensitivities used for reconstructing the images
recB0 = rmfield(recB0,'S');%Remove to save memory

%%% SOLVE FOR B0 FIELD
fprintf('Estimating B0 map\n');
recB0 = solveB07T(recB0);
%recB0 = rmfield(recB0,'x');%Remove to save memory

%%% SAVE RAW DATA
if writeRAW
    if exist(strcat(fileName,'.mat'),'file')
        load(fileName,'rec');
        if exist('rec','var')
            recB0.y=[];recB0.N=[];recB0.Assign=[];%Remove raw data to avoid storing twice
        end
        save(strcat(fileName,'.mat'),'recB0','-v7.3','-append');%If this file is also used for something else, dont'r overwrite .mat file    
    else
        save(strcat(fileName,'.mat'),'recB0','-v7.3');
    end
    fprintf('File saved:\n   %s\n', fileName);
end

%%% SAVE NIFTI
if writeNIIFlag
    recB0.Par.Mine.Modal=3;%Writing data to Re-F0 = B0 maps
    recB0.Plan.Suff = ''; recB0.Plan.SuffOu='';recB0.Names.pathOu = pathOu;
    recB0.Names.Name = file;
    recB0.Plan.Types{11}='B0';recB0.Plan.TypeNames{11}='B0';recB0.Dyn.Typ2Wri=zeros(1,50);recB0.Dyn.Typ2Wri(11)=1;recB0.Dyn.Typ2Rec=vertcat(recB0.Dyn.Typ2Rec,11);
    if writeNIIFlag>1;recB0.x=x;recB0.Plan.Types{10}='x';recB0.Plan.TypeNames{10}='Aq';recB0.Dyn.Typ2Wri(10)=1;recB0.Dyn.Typ2Rec=vertcat(recB0.Dyn.Typ2Rec,10);end
    recB0.Alg.OnlyJSON = 0;recB0.Fail=0;
    recB0.Par.Mine.Proce=1;recB0.Par.Mine.Nat=1;recB0.Alg.OverDec=ones(1,3);recB0.Par.Mine.Signs=[];
    writeData(recB0);
end
