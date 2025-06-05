function sTInfo=scanTimes(sTInfo,path,pathid,modal,rootOu,rootIn,rootSh,series,cpFlag)

%SCANTIMES   Obtains scanning times from raw data
%   STINFO=FILENAMECONVERSION(STINFO,PATH,{MODAL},{ROOTOU},{ROOTIN})
%   * STINFO is the scan time info to be updated
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
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
%   * {PATHID} contains the subject and session ids
%   * {CPFLAG} indicates whether to copy (1) or to move (0) the data. It
%   defaults to 0
%

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('\nChanging file names of study %s\n',path);
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

sTInfoC{1}=pathid{1};%Participant ID
sTInfoC{2}=pathid{2};%Scan ID
pathSpl=strsplit(path,'/');
sTInfoC{3}=pathSpl{1};%Date
pathSpl=strsplit(pathSpl{2},'_');
sTInfoC{4}=pathSpl{1};
pathSpl=strsplit(pathIn,'/');
sTInfoC{5}=pathSpl{6};
sTInfoC{6}=[];
sTInfoC{7}=[];


structFile=dir(sprintf('%s/*.raw',pathIn));
structFile={structFile.name};
for l=1:length(structFile)
    structFile=structFile([1 end]);
    for n=1:2
        pathSpl=strsplit(structFile{n},'_');
        pathSpl=pathSpl{3};
        sTInfoC{5+n}=[pathSpl(1:2) '.' pathSpl(3:4) '.' pathSpl(5:end)];        
    end
end
sTInfo=[sTInfo;sTInfoC];