
function recS = dat2Se(refName, refBName, supportReadout, writeNIIFlag, writeRAWFlag, pathOu, isPilotTone)

%DAT2SE processes a Siemens .dat file from a low-resolution acquisition and calculates the sensitivity maps using the ESPIRiT algorithm. 
%   [RECS]=DAT2SE(FILENAME,{UPPORREADOUT},{WRITENIIFLAG},{WRITERAWFLAG},{PATHOU},{ISPILOTTONE})
%   * REFNAME is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * {SUPPORTREADOUT} is the ralative range of the readout FOV to be extracted (between 0 and 1).
%   * {WRITENIIFLAG} is a flag to write the B0 map to NIFTI files (1) or to write both maps and raw images (2). Defaults to 1.
%   * {WRITERAWFLAG} is a flag to write the reconstruction structure to a .mat file.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the multiple TE file is stored.
%   * {ISPILOTTONE} indicated whether Pilot Tone was used during this acquisition. This activates a median filter across k-space to detect outliers. Poor performance, so only use when applicable.
%   ** REC is the reconstruction object.
%

[pathOuTemp,file,suff]=fileparts(refName);
if ~strcmp(refBName,'') && ~isempty(refBName); [~,fileB,~]=fileparts(refBName); end

fprintf('=====  Converting %s.dat file to coil sensitivities  ====\n',file);
if ~strcmp(refBName,'') && ~isempty(refBName); fprintf('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT      ===>   Using body coil from acquisition %s\n',fileB);end

if nargin<3 || isempty(supportReadout); supportReadout=[];end
if nargin<4 || isempty(writeNIIFlag); writeNIIFlag=0;end
if nargin<5 || isempty(writeRAWFlag); writeRAWFlag=1;end
if nargin<6 || isempty(pathOu); pathOu=[];end
if nargin<7 || isempty(isPilotTone); isPilotTone=0;end
if isempty(pathOu); pathOu = pathOuTemp;end
logFlag = 1;%~isempty(pathOu);

%%% INVERT TWIX DATA
writeNIIFlagRec = (writeNIIFlag>1);%Writing NIFTI of the raw coil images
writeRAWFlagRec = 1;
removeOversampling=1;
if exist(strcat(pathOuTemp,file,'.mat'),'file')
    load(strcat(pathOuTemp,file,'.mat'),'rec');
    if exist('rec','var')
        recS=rec;%Attention: when recS contains gpuArray variables and you re-use dat2Rec.m, this gives the following error Data no longer exists on the GPU.
        rec=[];
    else
        recS = dat2Rec(refName, supportReadout, writeNIIFlagRec, writeRAWFlagRec, [], isPilotTone, removeOversampling);
    end
else
    recS = dat2Rec(refName, supportReadout, writeNIIFlagRec, writeRAWFlagRec, [], isPilotTone, removeOversampling); 
end
if ~strcmp(fileB,'') && ~isempty(fileB)
     recB = dat2Rec(refBName, supportReadout, writeNIIFlagRec, writeRAWFlagRec, [], isPilotTone, removeOversampling); 
else
     recB = []; 
end

%%% ACTIVATE LOGGING
if logFlag
    if ~isempty(pathOu); logDir = strcat(pathOu,filesep,'Re-Se_Log');else;logDir='Re-Se_Log';end
    logName =  strcat( logDir,filesep, file, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~exist(logDir,'dir'); mkdir(logDir);end
    diary(logName) %eval( sprintf('diary %s ', logName))
end

%%% DEAL WITH MULTIPLE ECHOES
fprintf('Dealing with muliple echoes: Taking the first echo.\n');
nD = numDims(recS.y);
if nD>4; recS.y = permute(recS.y,[1:4 nD 5:(nD-1) ] );end%Different echoes stored in 5th dimension
recS.y = dynInd(recS.y,1,5);%Extract 1st echo to have minimal phase
%TO DO: if multiple echoes, use phase unwrapping and B0 mapping to remove
%B0 phase and only use transmit/receive phase in ESPIRiT

%%% SOLVE FOR COIL SENSITIVITIES
fprintf('Estimating coil sensitivities.\n');
recS.Names.pathOu = pathOu;
recS.Names.Name = file;
recS.Plan.Suff=''; recS.Plan.SuffOu='';

recS=solveSensit7T(recS, recB);
if logFlag;  diary off;end

%%% SAVE RAW DATA
if writeRAWFlag
    if exist(strcat(pathOuTemp,file,'.mat'),'file')
        load(strcat(pathOuTemp,file,'.mat'),'rec');
        if exist('rec','var')
            recS.y=[];recS.N=[];recS.Assign=[];%Remove raw data to avoid storing twice
        end
        save(strcat(pathOuTemp,file,'.mat'),'recS','-append');%If this file is also used for something else, dont'r overwrite .mat file    
    else
        save(strcat(pathOuTemp,file,'.mat'),'recS','-v7.3');
    end
    fprintf('Raw file saved:\n   %s\n', file);
end

%%% SAVE NIFTI
if writeNIIFlag>0
    fprintf('Writing NIFTI files.\n');
    pathOuNII = strcat(pathOu, '/Re-Se/'); if ~exist( pathOuNII,'dir');mkdir(pathOuNII);end
    fileSave = strcat( pathOuNII,file);
    xW=[];xW{1} = recS.S; 
    MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; 
    MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
    writeNII(fileSave, {'Se'},xW, MSW, MTW);
end

end