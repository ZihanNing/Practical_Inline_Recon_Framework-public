function fileNameConversion(path,pathid,modal,rootOu,rootIn,rootSh,series,cpFlag)

%FILENAMECONVERSION   Performs filenaming conversion to the dHCP format
%(based on BIDS)
%   FREC=FILENAMECONVERSION(PATH,{MODAL},{ROOTOU},{ROOTIN})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {PATHID} contains the subject and session ids
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%   * {ROOTIN} is the folder where the raw data has been mounted. It
%   defaults to Data/rawSource (relative to the user's home folder)
%   * {ROOTSH} is the folder where the sharing data has to be placed. It
%   defaults to ROOTOU
%   * {SERIES} restricts the reconstructions to a specific set of series
%   * {CPFLAG} indicates whether to copy (1) or to move (0) the data 
%   (with overwritting)
%

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('\nChanging file names of study %s into sub-%s/ses-%s\n',path,pathid{1},pathid{2});
%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;

if nargin<8 || isempty(cpFlag);cpFlag=0;end

if nargin<3;modal=[];end
if nargin<4 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
if nargin<6 || isempty(rootSh);rootSh=rootOu;end
pathOu=fullfile(rootOu,path);
%DETECTING THE RAW FOLDER
if nargin<5 || isempty(rootIn);rootIn=fullfile(filesep,'home',user,'Data','rawSource');end
pathIn=rawFolderDetection(path,rootIn);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end
%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');parseLabFolder(path,[],rootOu,rootIn);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin<7 && ~isempty(series);nV=nV(ismember(nV,series));end
interveneNext=0;

protOu.A_Modal=[];protOu.B_Recon=[];protOu.C_FileName=[];
for n=nV        
    if interveneNext~=2        
        rec.Rec=1;
        rec.Names.Name=strtrim(prot.B_FileName(n,:));
        matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));rawFile=fullfile(pathIn,sprintf('%s.raw',rec.Names.Name));
        rec.Names.rawFile=rawFile;rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.pathIn=pathIn;rec.Names.path=path;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;
        if exist(rawFile,'file')
            %FILENAMING CONVERSION
            load(matFile);
            rec.Par=Par;Par=[];            
            a=strsplit(rec.Names.Name,'_');
            acqLabel=a{end};acqLabel=acqLabel(1:end-2);
            serLabel=sprintf('%02d',str2double(a{end-2}));
            niiPref=fullfile(pathOu,numbe2Modal(rec.Par.Mine.Modal),sprintf(rec.Names.Name));
            niiFile=sprintf('%s_Aq.json',niiPref);
            recData=0;
            if exist(niiFile,'file')
                recData=1;
                structFile=dir(sprintf('%s_*',niiPref));
                structFile={structFile.name};
                for l=1:length(structFile)
                    a=strsplit(structFile{l},'_');   
                    recLabel=a{end};                    
                    outFile=sprintf('%s/sub-%s/ses-%s/%s/sub-%s_ses-%s_%s_%s-%s_%s',rootSh,pathid{1},pathid{2},numbe2Modal(rec.Par.Mine.Modal),pathid{1},pathid{2},numbe2Modal(rec.Par.Mine.Modal),serLabel,acqLabel,recLabel);
                    inpFile=fullfile(pathOu,numbe2Modal(rec.Par.Mine.Modal),structFile{l});
                    outFold=fileparts(outFile);
                    if ~exist(outFold,'dir');mkdir(outFold);end
                    if ~cpFlag;movefile(inpFile,outFile);else copyfile(inpFile,outFile);end
                end
            end          
            %PHYSLOG ANONIMIZATION
            logFileIn=sprintf('%s/SCANPHYSLOG_%s.log',pathIn,sprintf(rec.Names.Name));  
            logFileOu=sprintf('%s/sub-%s/ses-%s/%s/sub-%s_ses-%s_%s_%s-%s.log',rootSh,pathid{1},pathid{2},numbe2Modal(rec.Par.Mine.Modal),pathid{1},pathid{2},numbe2Modal(rec.Par.Mine.Modal),serLabel,acqLabel);                
            if exist(logFileIn,'file')         
                %READ
                fid = fopen(logFileIn,'r');            
                A=cell(1,3);
                for l=1:2;A{l}=fgetl(fid);end;A{3}=fread(fid);      
                
                %ANONIMIZE
                anon='St. Thomas Hospital';
                ind=regexp(A{1},anon);
                A{1}(ind:ind+length(anon)-1)='X';
                A{2}(4:end)='X';              
                
                %WRITE
                outFold=fileparts(logFileOu);
                if ~exist(outFold,'dir');mkdir(outFold);end
                fid=fopen(logFileOu, 'w');
                for l=1:2;fprintf(fid,'%s\n',A{l});end;fwrite(fid,A{3});           
                fclose(fid);
            end
            
            %PROTOCOL INFORMATION
            fileOu=sprintf('sub-%s_ses-%s_%s_%s-%s',pathid{1},pathid{2},numbe2Modal(rec.Par.Mine.Modal),serLabel,acqLabel);                        
            protOu.A_Modal=[protOu.A_Modal;{numbe2Modal(rec.Par.Mine.Modal)}];   
            protOu.B_Recon=[protOu.B_Recon;uint8(recData)];
            protOu.C_FileName=[protOu.C_FileName;{fileOu}];
        end
    end
    rec=[];
end

%WRITTING THE PROTOCOL FILE
protFile=sprintf('%s/sub-%s/ses-%s/sub-%s_ses-%s_prot.txt',rootSh,pathid{1},pathid{2},pathid{1},pathid{2});
outFold=fileparts(protFile);
if ~exist(outFold,'dir');mkdir(outFold);end
tdfwrite(protFile,protOu);

%DELETING FOLDER
if ~cpFlag
    rmdir(pathOu,'s');
    day=strsplit(path,'/');day=day{1};day=fullfile(rootOu,day);
    if isempty(ls(day));rmdir(day);end
end
