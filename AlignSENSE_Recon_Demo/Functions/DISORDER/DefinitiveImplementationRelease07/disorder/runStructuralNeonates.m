clearvars
addpath(genpath(fileparts(mfilename('fullpath'))));
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
path=studiesDetection(1);
permMount='/home/lcg13/Data/pnrawDe';%Permanent mount
[~,hostname]=system('hostname');hostname=strtrim(lower(hostname));%Computer name
load('neonatalPermutation.mat');

inpRootPath=[];%Modify to specify the input (raw) path
outRootPath=fullfile(permMount,'ReconstructionsRelease07','dHCPNeonatal');%Modify to specify the output (reconstructions) path
modality=[5 6 7];%Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%Modify to specify a particular series (with indexes as row position in prot.txt)

if strcmp(hostname,'perinatal130-pc');n0=1;
elseif strcmp(hostname,'gpubeastie02-pc');n0=2;
elseif strcmp(hostname,'gpubeastie03-pc');n0=3;
end
ne=length(path);%Final case
nd=3;%Jump

%neonatalPermutation=1:894;
for l=n0:nd:ne
    n=neonatalPermutation(l);
    if ~exist(outRootPath,'dir');error('Output path %s not defined: check mounts',outRootPath);end
    pathOu=fullfile(outRootPath,path{n},'An-S2');
    if ~exist(pathOu,'dir')
        fprintf('\nStructural %d of %d (id %d)\n',l,length(path),n);
        protoPipeline(path{n},modality,outRootPath,inpRootPath,series);%RECONSTRUCTING
        procePipelineAssembleSVR(path{n},5:6,outRootPath,series);%SVR
        changePermissions(path{n},outRootPath);
    end
end
