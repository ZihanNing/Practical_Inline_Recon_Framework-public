addpath(genpath('.'))

if useGPU;gpuDevice(2);end

idData=1;%[2 5];%[1 3 4];%[];%1:3;%[];%Reconstruction to run, if empty all

%FOLDER AND FILES TO RUN
siemensStudies;
if ~isempty(idData)
    dataIn=dataIn(idData);
    surfIn=surfIn(idData);
end

%dataIn=[];

%RECONSTRUCTION PIPELINE
for n=1:length(dataIn)
    rec=[];
    rec=reconPipeline(caseIn,bodyIn,surfIn{n},dataIn{n});
end

return

