
function [recB1, TWB1] = dat2B1(refName, fileName, supportReadout, writeNIIFlag, pathOu, isPilotTone)

%DAT2B1 processes a Siemens .dat file from an Actual Flip Angle (AFI) acquisition and calculated B1 and B0 maps (if multiple echoes are provided). 
%   []=DAT2B1(REFNAME,FILENAME,{UPPORREADOUT},{WRITENIIFLAG},{PATHOU},{ISPILOTTONE})
%   * REFNAME is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * FILENAME is the name of the raw data (.dat) of the AFI acquisition.
%   * {SUPPORTFOV} is the ralative range of the readout FOV to be extracted.
%   * {WRITENIIFLAG} is a flag to write B1 and B0 maps to NIFTI files (1) or to write both maps and raw images (2). Defaults to 1.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the AFI file is stored.
%   * {isPilotTone} indicated whether Pilot Tone was used during this acquisition. This activates a median filter across k-space to detect outliers. Poor performance, so only use when applicable.
%

[pathOuTemp,file]=fileparts(fileName);
fprintf('Parsing TWIX object to reconstruction structure:\n%s\n',file);

if nargin<3 || isempty(supportReadout); supportReadout=[];end
if nargin<4 || isempty(writeNIIFlag); writeNIIFlag=2;end
if nargin<5 || isempty(pathOu); pathOu=[];end
if nargin<6 || isempty(isPilotTone); isPilotTone=0;end
if isempty(pathOu); pathOu = pathOuTemp;end

%%% INVERT TWIX DATA
TWB1=mapVBVD(strcat(fileName,'.dat'));     
if length(TWB1)>1
    TWB1{2}.noise = TWB1{1}.noise;%Take noise from first one and discard all the rest
    TWB1 = TWB1{2};
end

if isPilotTone
    recB1=TWIX2rec_PilotTone(TWB1,supportReadout);
else
    recB1=TWIX2rec(TWB1,supportReadout);
end

%TODO: include dat2Rec for consistency and json writing

%%% ADD PARAMETERS FOR AFI ACQUISITION
recB1.Par.Labels.SliceGaps = 0; %Used in solveB1.m - disabled here
recB1.Par.Scan.FastImgMode='AFI'; %Hardcoded now assuming only acquiring AFI B1 maps
recB1.Par.Labels.RepetitionTime = [TWB1.hdr.MeasYaps.alTR{1}/1000  TWB1.hdr.MeasYaps.alTR{1}/1000*TWB1.hdr.MeasYaps.sWipMemBlock.alFree{11}]; %;%2needed for AFI calculation. TWB1.hdr.MeasYaps.sWipMemBlock.alFree(11) is TR2/TR1
recB1.Par.Scan.AcqVoxelSize = recB1.Enc.AcqVoxelSize;
recB1.Dyn.Typ2Rec=[];

%%% LOAD SENSITIVITY MAPS FOR RECONSTRUCTION
if ~exist(strcat(refName,'.mat'),'file'); dat2Se(refName);end
ss = load(strcat(refName,'.mat'));recS =ss.recS; ss=[];
N = size(recB1.y);
MTS = recS.Par.Mine.APhiRec;MTB1 = recB1.Par.Mine.APhiRec;
recB1.S = mapVolume (recS.S, ones(N(1:3)), MTS, MTB1);recS=[];

%%% RECONSTRUCT IMAGES
x= SWCC (recB1.y, recB1.S);%Simple PI recon
nD = numDims(x);
x = permute(x,[1:3 nD 4:(nD-1) ] );%Different echoes stored in 4th dimension
if size(x,4)>2%AFI also contains 3 echoes in first TR to estimate B0
    recB1.x = dynInd( x, [1 4],4);%First echo of each TR
    estimateB0=1;
else
    recB1.x = dynInd(x, 1:2,4);
    estimateB0=0;
end
recB1.Alg.Sensitivities = refName;%Store the name of the sensitivities used for reconstructing the images
recB1 = rmfield(recB1,'S');%Remove to save memory

%%% SOLVE FOR B1 FIELD
fprintf('Estimating B1 map\n');
recB1=solveB1(recB1);%Now normalised B1 map is concatenated in recB1.B

%%% SOLVE FOR B0
if estimateB0
    recB1.Par.Labels.TE = TWB1.hdr.Meas.alTE(1:3)/1000;%In ms
    recB1.Alg.parS=[];recB1.Alg.parS.conComp=2;%Set to 0 if multiple objects in image
    recB1.Alg.parU.UnwrapMeth='CNCG';
    recB1.Alg.useTE3=length(recB1.Par.Labels.TE)>2;
    recB1.x = dynInd(x,1:3,4);
    fprintf('Estimating B0 map\n');
    recB1=solveB07T(recB1);%Now B0 map is concatenated in recB1.B0   
end

%%% SAVE RAW DATA
recB1 = rmfield(recB1,'x');%Remove to save memory
save(strcat(fileName,'.mat'),'recB1','-v7.3');    
fprintf('File saved:\n   %s\n', fileName);

%%% WRITE METADATA TO JSON
recJSON = recB1;
%Remove all big arrays
recJSON.y=[];
recJSON.S=[];recJSON.N=[];recJSON.Assign=[];recJSON.B=[];recJSON.B0=[];
recJSON.x=[];recJSON.M=[];
if isfield(recJSON.Par.Labels,'Shim') && isfield(recJSON.Par.Labels.Shim, 'crossTerms');recJSON.Par.Labels.Shim=rmfield(recJSON.Par.Labels.Shim,'crossTerms');end

%Save JSON
savejson('',recJSON,sprintf('%s.json',fileName));%writeJSON has specific fields 
    
%%% SAVE NIFTI
if writeNIIFlag
    recB1.Par.Mine.Modal=4;%Writing data to Re-F1 = B1 maps
    recB1.Plan.Suff = ''; recB1.Plan.SuffOu='';recB1.Names.pathOu = pathOu;
    recB1.Names.Name = file;
    recB1.Plan.Types{11}='B';recB1.Plan.TypeNames{11}='B1';recB1.Dyn.Typ2Wri=zeros(1,50);recB1.Dyn.Typ2Wri(11)=1;
    if estimateB0; recB1.Plan.Types{9}='B0';recB1.Plan.TypeNames{9}='B0';recB1.Dyn.Typ2Wri(9)=1;recB1.Dyn.Typ2Rec=vertcat(recB1.Dyn.Typ2Rec,9);end
    if writeNIIFlag>1; recB1.x=x; recB1.Plan.Types{10}='x';recB1.Plan.TypeNames{10}='Aq';recB1.Dyn.Typ2Wri(10)=1;recB1.Dyn.Typ2Rec=vertcat(recB1.Dyn.Typ2Rec,10);end
    recB1.Alg.OnlyJSON = 0;recB1.Fail=0;
    recB1.Par.Mine.Proce=1;recB1.Par.Mine.Nat=1;recB1.Alg.OverDec=ones(1,3);recB1.Par.Mine.Signs=[];
    writeData(recB1);
end

end