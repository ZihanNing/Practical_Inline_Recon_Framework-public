clearvars
addpath(genpath(fileparts(mfilename('fullpath'))));
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
[path,pathid]=studiesDetection(1);
permMount='/home/lcg13/Data/pnrawDe';%Permanent mount
[~,hostname]=system('hostname');hostname=strtrim(lower(hostname));%Computer name
load('neonatalPermutation.mat');

inpRootPath=[];%Modify to specify the input (raw) path
outRootPath=fullfile(permMount,'ReconstructionsRelease07','dHCPNeonatal');%Modify to specify the output (reconstructions) path
modality=[];%Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%Modify to specify a particular series (with indexes as row position in prot.txt)

if strcmp(hostname,'perinatal130-pc');n0=1;
elseif strcmp(hostname,'gpubeastie02-pc');n0=2;
elseif strcmp(hostname,'gpubeastie03-pc');n0=3;
end
ne=length(path);%Final case
nd=3;%Jump

%neonatalPermutation=1:894;
for l=n0:nd:ne
    lock_check;
    n=neonatalPermutation(l);
    if ~exist(outRootPath,'dir');error('Output path %s not defined: check mounts',outRootPath);end
    pathOuA=fullfile(outRootPath,path{n},'Re-Se');
    pathOuB=fullfile(outRootPath,sprintf('sub-%s/ses-%s',pathid{n}{1},pathid{n}{2}));
    if ~exist(pathOuA,'dir') && ~exist(pathOuB,'dir')
        fprintf('\nRunning %d of %d (id %d)\n',l,length(path),n);
        protoPipeline(path{n},modality,outRootPath,inpRootPath,series);%RECONSTRUCTING
        procePipeline(path{n},10,outRootPath,series);%PREPROCESSING
        procePipelineAssembleSVR(path{n},5:6,outRootPath,series);%SVR       
        fileNameConversion(path{n},pathid{n},modality,outRootPath,inpRootPath,outRootPath,series,0);
        pathPerm=sprintf('sub-%s/ses-%s',pathid{n}{1},pathid{n}{2});
        changePermissions(pathPerm,outRootPath);
    end
end
