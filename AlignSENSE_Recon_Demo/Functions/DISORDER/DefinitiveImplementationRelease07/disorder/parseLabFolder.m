function prot=parseLabFolder(path,modal,rootOu,rootIn)

%PARSELABFOLDER   Parses a raw data folder using reconframe. It creates a
%set of .mat files which are placed in the subfolder ZZ-HE of the
%rootOu/path folder and a prot.txt file with the corresponding set of 
%detected modalities and filenames which is placed in the rootOu/path folder.
%   PARSELABFOLDER(PATH,{MODAL},{ROOTOU},{ROOTIN})
%   * PATH is the relative path to raw data
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%   * {ROOTIN} is the folder where the raw data has been mounted. It
%   defaults to Data/rawSource (relative to the user's home folder)
%

addpath(genpath(fileparts(mfilename('fullpath'))));

%SET DEFAULT PARAMETERS
[user,versCode,versSave,pathRF,pathSt]=versRecCode;
if ~exist('rootOu','var') || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);
if ~exist('modal','var');modal=[];end
if ~exist('rootIn','var') || isempty(rootIn);rootIn=fullfile(filesep,'home',user,'Data','rawSource');end
%DETECTING THE RAW FOLDER
pathIn=rawFolderDetection(path,rootIn);
fprintf('Parsing raw data folder %s\n',pathIn);tsta=tic;

%ADD THE RECONFRAME CODE TO THE PATH
assert(logical(exist(pathRF,'dir')),'ReconFrame code not found at %s',pathRF);
addpath(genpath(pathRF));

%INITIALIZING VARIABLES AND DESTINATION FOLDERS
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(headFolder,'dir')
    mkdir(headFolder);
else
    fprintf('The destination folder already existed, overwritting parsed content\n');
    if exist(protFile,'file');delete(protFile);end
    delete(fullfile(headFolder,'*'));
end
if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end
prot.A_Modal=[];prot.B_FileName=[];

%READ SPECIFIC STUDY INFORMATION AND DETECT CASE OF INTEREST
studFile=fullfile(pathSt,'cases.txt');
if exist(studFile,'file');warning('off','MATLAB:namelengthmaxexceeded');stud=tdfread(studFile,' ');warning('on','MATLAB:namelengthmaxexceeded');else stud=[];end
curStud.Stud='None';curStud.IssuId='None';curStud.IssuPa=0;
if ~isempty(stud)
    for n=1:size(stud.Case,1)
        if strcmp(strtrim(stud.Case(n,:)),path)
            curStud.Stud=strtrim(stud.Study(n,:));curStud.IssuId=strtrim(stud.IssueId(n,:));curStud.IssuPa=stud.IssuePar(n,:);break;
        end
    end 
    studC{1}='dHCPNeoBra';studC{3}='dHCPTodBra';studC{5}='dHCPEpiBra';
    for p=[1 3 5]
        studNeon=studiesDetection(p,[],rootOu);%Neonatal studies don't need to be in the file unless they are exceptional...
        for n=1:length(studNeon)
            if strcmp(studNeon{n},path)   
                curStud.Stud=studC{p};
                break;
            end
        end
    end
end
curStudMod=curStud;

%TRAVERSE THE RAW DATA FOLDER, EXTRACT THE LAB INFO FROM THE FILES, DUMP IT
%TO A MATLAB STRUCTURE AND WRITE THIS STRUCTURE TO FILE
files=dir(fullfile(pathIn,'*.lab'));

alreadyref=0;%Indicating that there is already a valid reference scan
for f=1:length(files)
    files(f).name=files(f).name(1:end-4);
    targetFile=files(f).name;
    sourceFile=fullfile(pathIn,sprintf('%s.lab',targetFile));
    destinFile=fullfile(headFolder,sprintf('%s.mat',targetFile));    
    if exist(sourceFile,'file')
        %fprintf('Parsing file %s\n',sourceFile);
        Par=[];
        try
            MR=MRecon(sourceFile);
        catch
            fprintf('Information for %s seems corrupted\n',sourceFile);
            continue
        end
        Param=MR.Parameter;
        %fprintf('Scan technique: %s\n',Param.Scan.Technique);
        %fprintf('Fast imaging mode: %s\n',Param.Scan.FastImgMode);
        %fprintf('Scan mode: %s\n',Param.Scan.ScanMode);
        if isempty(Param.Scan.Technique)
            fprintf('No scan technique stored for series %s\n',files(f).name);
        elseif isempty(Param.Scan.ScanType)%Introduced for ma_09012017_1359315_5_2_dhcp8mpragesenseV4
            fprintf('No scan type stored for series %s\n',files(f).name);
        elseif Param.Scan.Stacks>=3
            %fprintf('Series %s comprised of %d stacks, probably a survey\n',files(f).name,Param.Scan.Stacks) 
        %elseif isfield(Param.Labels,'RefScan') && strcmp(Param.Labels.RefScan,'CoilSurvey')
            %fprintf('Code not prepared to interpret coil surveys\n');
        elseif (strcmp(Param.Scan.Technique,'T1FFE') && strcmp(Param.Scan.FastImgMode,'None') && strcmp(Param.Scan.ScanMode,'3D') && strcmp(Param.Scan.AngioMode,'No') && strcmp(Param.Scan.FlowComp,'no') ...
            && Param.Scan.SENSEFactor(1,2)==1 && length(Param.Labels.CoilNrsPerStack)==2) %|| (isfield(Param.Labels,'RefScan') && strcmp(Param.Labels.RefScan,'SenseRefFast'))
            %Have introduced the check on the COILNRSPERSTACK, which may be problematic in other scans, but was necessary for 2018_07_17/LE_28230
            Par=lab2Mat(Param,2,modal,curStudMod);%SENSE            
        elseif strcmp(Param.Scan.Technique,'T1FFE') && strcmp(Param.Scan.FastImgMode,'None') && strcmp(Param.Scan.ScanMode,'M2D') && strcmp(Param.Scan.PartialEcho,'no') && Param.Labels.ZReconLength==2
            Par=lab2Mat(Param,3,modal,curStudMod);%B0
        elseif (strcmp(Param.Scan.Technique,'T1TFE') || strcmp(Param.Scan.Technique,'T1FFE') || strcmp(Param.Scan.Technique,'TFE')) && (strcmp(Param.Scan.FastImgMode,'TFE') || strcmp(Param.Scan.FastImgMode,'None')) && Param.Encoding.NrMixes==2 && ~strcmp(Param.Scan.ScanMode,'2D')
            Par=lab2Mat(Param,4,modal,curStudMod);%B1
        elseif strcmp(Param.Scan.Technique,'TSE') && strcmp(Param.Scan.FastImgMode,'TSE') && (strcmp(Param.Scan.ScanMode,'MS') || strcmp(Param.Scan.ScanMode,'M2D')) && size(Param.Scan.FOV,1)==1%To prevent surveys
            Par=lab2Mat(Param,5,modal,curStudMod);%MST2            
        elseif (strcmp(Param.Scan.Technique,'TIR') || strcmp(Param.Scan.Technique,'T1TFE') || strcmp(Param.Scan.Technique,'B-FFE') || strcmp(Param.Scan.Technique,'B-TFE')) && (strcmp(Param.Scan.FastImgMode,'TSE') || strcmp(Param.Scan.FastImgMode,'TFE') || strcmp(Param.Scan.FastImgMode,'None')) && (strcmp(Param.Scan.ScanMode,'MS') || strcmp(Param.Scan.ScanMode,'M2D')) && Param.Encoding.NrMixes==1
            Par=lab2Mat(Param,6,modal,curStudMod);%MST1     
        elseif (strcmp(Param.Scan.FastImgMode,'None') || strcmp(Param.Scan.FastImgMode,'TFE') || strcmp(Param.Scan.FastImgMode,'TSE')) && (strcmp(Param.Scan.Technique,'T1FFE') || strcmp(Param.Scan.Technique,'TSE') || strcmp(Param.Scan.Technique,'T1TFE') || strcmp(Param.Scan.Technique,'BALANCEDTFE') || strcmp(Param.Scan.Technique,'B-TFE') || strcmp(Param.Scan.Technique,'B-FFE') || strcmp(Param.Scan.Technique,'TIR')) && strcmp(Param.Scan.ScanMode,'3D')
            Par=lab2Mat(Param,7,modal,curStudMod);%3D-Technique refers to TSE-MPRAGE/SPGR-SSFP-TIR (MPRAGE and SPGR are not distinguishable at this point, it is done later based on the TFE factor)
        %Modality number 8 has been removed. It may be discretionally assigned in the future
        elseif strcmp(Param.Scan.Technique,'FFEEPI') || strcmp(Param.Scan.Technique,'FEEPI') || ((strcmp(Param.Scan.Technique,'DIFFSE') || strcmp(Param.Scan.Technique,'DwiSE') || strcmp(Param.Scan.Technique,'SEEPI')) && length(Param.Parameter2Read.extr1)==1 && length(Param.Parameter2Read.extr2)==1)
            Par=lab2Mat(Param,9,modal,curStudMod);%fMRI            
        elseif (strcmp(Param.Scan.Technique,'DIFFSE') || strcmp(Param.Scan.Technique,'DwiSE')) && (~strcmp(curStud.IssuId,'DiscardDWI') || f>curStud.IssuPa)
            %if strcmp(curStud.IssuId,'DiscardDWI') && f<=curStud.IssuPa;curStudMod.Stud='None';else curStudMod.Stud=curStud.Stud;end
            Par=lab2Mat(Param,10,modal,curStudMod);%dMRI
        end        
        if ~isempty(Par)
            Par=computeRAF(MR,Par);           
            if ~isempty(Par) && isfield(Par.Mine,'APhiAcq') && ~isempty(Par.Mine.APhiAcq) && (alreadyref || Par.Mine.Modal==2)
                alreadyref=1;
                targetFileCell{1}=targetFile;
                prot.B_FileName=[prot.B_FileName;targetFileCell];
                %This would probably be better from MatlabR2016b
                %prot.B_FileName=vertcat(prot.B_FileName,string(targetFile));
                prot.A_Modal=vertcat(prot.A_Modal,uint8(Par.Mine.Modal));
                save(destinFile,'Par',versSave); 
            end
        end
        delete(MR);
    end
end

%WRITTING THE PROTOCOL FILE AND SOME STATS
tdfwrite(protFile,prot);
sizeLab=0;sizeMat=0;
for n=1:length(prot.B_FileName)
    fileAux=dir(fullfile(headFolder,sprintf('%s.mat',prot.B_FileName{n})));    
    %fileAux=dir(sprintf('%s/%s.mat',headFolder,prot.B_FileName(n)));%When using strings
    sizeMat=sizeMat+fileAux.bytes;
    fileAux=dir(fullfile(pathIn,sprintf('%s.lab',prot.B_FileName{n})));
    %fileAux=dir(sprintf('%s/%s.lab',pathIn,prot.B_FileName(n)));%When using strings
    sizeLab=sizeLab+fileAux.bytes;
end
%fprintf('Size of .mat files: %s. Ratio of compression over .lab: %.2f %%\n',byteSize(sizeMat),100*sizeMat/sizeLab);
tend=toc(tsta);fprintf('Time parsing raw data folder: %.3f s\n\n',tend);
