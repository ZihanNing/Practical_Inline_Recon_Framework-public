clearvars
addpath(genpath(fileparts(mfilename('fullpath'))));
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
[path,pathid]=studiesDetection(1);
permMount='/home/lcg13/Data/pnrawDe';%Permanent mount
load('neonatalPermutation.mat');

inpRootPath=[];%Modify to specify the input (raw) path
outRootPath=fullfile(permMount,'ReconstructionsRelease07','dHCPNeonatal','INFOPriv');%Modify to specify the output (reconstructions) path
modality=[];%Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%Modify to specify a particular series (with indexes as row position in prot.txt)

n0=1;
ne=length(path);%Final case
nd=1;%Jump

%neonatalPermutation=1:894;
fileID=fopen(fullfile(outRootPath,'dHCPRelease07RandomSpreadsheet.txt'),'w');
for l=n0:nd:ne
    n=neonatalPermutation(l);    
    %fprintf(fileID,'%s %s\n',pathid{n}{1},pathid{n}{2});
    fprintf(fileID,'%d %s %s %s\n',n,pathid{n}{1},pathid{n}{2},path{n});    
end
fclose(fileID);

fileID=fopen(fullfile(outRootPath,'dHCPRelease07RandomList.txt'),'w');
for l=n0:nd:ne
    n=neonatalPermutation(l);    
    %fprintf(fileID,'%s %s\n',pathid{n}{1},pathid{n}{2});
    fprintf(fileID,'%s %s\n',pathid{n}{1},pathid{n}{2});
end
fclose(fileID);

