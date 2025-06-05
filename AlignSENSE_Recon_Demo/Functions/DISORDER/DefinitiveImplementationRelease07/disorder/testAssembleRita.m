
addpath(genpath(fileparts(mfilename('fullpath'))));

outRootPath='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/RitaMotionCorrection';
folde{1}='RitaMotionCorrection/Train';
folde{2}='RitaMotionCorrection/Test';

for n=1:length(folde)
    pathFold{n}=strcat(outRootPath,filesep,folde{n});
    if ~exist(pathFold{n},'dir');mkdir(pathFold{n});end
end

neonatalBrainData;
for p=1:length(path)
    pathIn=strcat(outRootPath,filesep,path{p},filesep,'An-Ve',filesep,'*mat');
    structFil=dir(pathIn);%Returns all the files and folders in the directory
    
    for n=1:length(structFil)
        pathIn=strcat(structFil(n).folder,filesep,structFil(n).name);
        pathOu=strcat(pathFold{floor((p-1)/5)+1},filesep,sprintf('rec%03d-%02d.mat',p,n));
        copyfile(pathIn,pathOu);
    end
end
