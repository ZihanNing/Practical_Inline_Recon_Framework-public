clearvars
addpath(genpath(fileparts(mfilename('fullpath'))));
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
[path,pathid]=studiesDetection(1);
permMount='/home/lcg13/Data/pnrawDe';%Permanent mount
load('neonatalPermutation.mat');

inpRootPath=[];%Modify to specify the input (raw) path
outRootPath=fullfile(permMount,'ReconstructionsRelease07','dHCPNeonatal');%Modify to specify the output (reconstructions) path
modality=[];%Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%Modify to specify a particular series (with indexes as row position in prot.txt)

%neonatalPermutation=1:894;
for l=1:length(path)
    n=neonatalPermutation(l);
    if ~exist(outRootPath,'dir');error('Output path %s not defined: check mounts',outRootPath);end
    pathOuA=fullfile(outRootPath,path{n});
    pathOuB=fullfile(outRootPath,sprintf('sub-%s/ses-%s',pathid{n}{1},pathid{n}{2}));
    if ~exist(pathOuA,'dir') && ~exist(pathOuB,'dir')
        fprintf('\nParsing %d of %d (id %d)\n',l,length(path),n);
        parseLabFolder(path{n},modality,outRootPath,inpRootPath);%PARSING RAW REQUIRES RECONFRAME LICENSE---WORKS IN perinatal130-pc AND gpubeastie03-pc
        changePermissions(path{n},outRootPath);
    end
end
