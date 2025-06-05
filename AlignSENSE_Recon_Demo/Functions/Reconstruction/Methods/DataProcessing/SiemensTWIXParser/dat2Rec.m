

function [rec, TW] = dat2Rec(fileName, supportReadout, writeNIIFlag, writeRAWFlag, pathOu, isPilotTone, removeOversampling, noiseName,resRec)
    
%DAT2REC converts a Siemens .dat file to a reconstruction structure with all the necessary data in it. 
%   [REC , TW] = DAT2REC(FILENAME,SUPPORTREADOUT,{WRITENIIFLAG},{WRITERAWFLAG},{PATHOU},{ISPILOTTONE},{REMOVEOVERSAMPLING},{NOISENAME},{RESREC})
%   * FILENAME is the name of the raw data (.dat) of the multiple TE acquisition.
%   * {SUPPORTREADOUT} is the ralative range of the readout FOV to be extracted (between 0 and 1).
%   * {WRITENIIFLAG} s a flag to write the raw coil images to NIFTI files (1). Defaults to 0.
%   * {WRITERAWFLAG} is a flag to write the reconstruction structure to a .mat file.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the .dat file is stored.
%   * {ISPILOTTONE} indicated whether Pilot Tone was used during this acquisition. For implications, see TWIX2Rec.m.
%   * {REMOVEOVERSAMPLING} indicated whether to remove Siemens oversampling.
%   * {NOISENAME} is the name of the file from which to use noise samples for channel de-correlation.
%   * {RESREC} is the resolution at which to downsample the reconstruction structure.
%   ** REC is the reconstruction object.
%   ** TW is the untouched TWIX object.
%       

[pathOuTemp,file] = fileparts(fileName);
fprintf('=====  Converting .dat file to a reconstruction structure  ====\n   %s\n', file);

if nargin<2 || isempty(supportReadout); supportReadout=[];end
if nargin<3 || isempty(writeNIIFlag); writeNIIFlag=0;end
if nargin<4 || isempty(writeRAWFlag); writeRAWFlag=0;end
if nargin<5 || isempty(pathOu); pathOu=[];end
if nargin<6 || isempty(isPilotTone); isPilotTone=0;end
if nargin<7 || isempty(removeOversampling); removeOversampling=1;end
if nargin<8 || isempty(noiseName); noiseName='';end
if nargin<9 || isempty(resRec); resRec=[];end
if isempty(pathOu); pathOu = pathOuTemp;end

%%% CONVERT .DAT FILE TO TWIX OBJECT
fprintf('Converting .dat file to a TWIX object.\n');
TW = mapVBVD(strcat(fileName,'.dat'));

%%% CONVERT TWIX TO A REC STRUCTURE
fprintf('Converting TWIX object to a reconstruction structure.\n');
if iscell(TW);TW=TW{end};end%Take last one by default as previous ones might be adjustment volumes
rec = TWIX2rec(TW, supportReadout, isPilotTone, [], removeOversampling, noiseName, fullfile(pathOu,file), resRec);
rec.Names.Name=file;
rec.Names.pathOu=pathOu;

%%% WRITE RAW
if writeRAWFlag
    rec = gatherStruct(rec);
    save(strcat(pathOu,filesep,file,'.mat'),'rec','-v7.3');    
    fprintf('Raw file saved:\n   %s\n', strcat(pathOu,filesep,file,'.mat'));
end

%%% SAVE NIFTI
if writeNIIFlag>0
    outDir = fullfile(pathOu,'Parsing'); if exist(outDir,'dir')~=2; mkdir(outDir);end
    xW=[];xW{1} = rec.y; 
    MSW=[];MSW{1} = rec.Enc.AcqVoxelSize; 
    MTW =[];MTW{1} = rec.Par.Mine.APhiRec;
    writeNII( strcat(outDir,filesep,rec.Names.Name) , {'y'},xW, MSW, MTW);
    fprintf('NIFTI file with raw coil images saved:\n   %s\n', strcat(outDir,filesep,rec.Names.Name));
end

%%% WRITE METADATA TO JSON FILE
writeJSONFlag=1;%Hard-coded to 1
if writeJSONFlag>0
    recJSON = rec;
    %Remove all big arrays - look at TWIX2Rec.m to see which big arrays are generated
    recJSON.y=[];
    recJSON.S=[];recJSON.N=[];recJSON.Assign=[];
    recJSON.x=[];recJSON.M=[];
    recJSON.Par.Labels.Shim.correctFact_mA2au=[];
%     recJSON.ACS=[];
    if isfield(recJSON,'Alg') && isfield(recJSON.Alg,'covMatrix');recJSON.Alg.covMatrix=[];end
    if isfield(recJSON,'PT')
        recJSON.PT.yProjRO=[];
        recJSON.PT.pSliceImage=[];
        recJSON.PT.pTimeTest=[];
    end
    recJSON.NY=size(rec.y);
    if isfield(recJSON.Par.Labels,'Shim') && isfield(recJSON.Par.Labels.Shim, 'crossTerms');recJSON.Par.Labels.Shim=rmfield(recJSON.Par.Labels.Shim,'crossTerms');end
    %Save JSON
    savejson('',recJSON,sprintf('%s.json',fullfile(pathOu,file)));%writeJSON has specific fields 
    fprintf('JSON file saved:\n   %s\n', fullfile(pathOu,file));
end


end