clearvars
addpath(genpath(fileparts(mfilename('fullpath'))));
%return
%GENERIC VALUES (BEASTIE01)
%inpRootPath='/pnraw01';
%outRootPath='/projects/perinatal/peridata/TestReconCode';
inpRootPath=[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
outRootPath=[];%'/home/lcg13/Data/pnrawDe/ReconstructionsDebug07/SlavaNewReconstruction';%'/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/fetalSAFETest';%[];%'/home/lcg13/Data/18ReplicaRelease06';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
%outRootPath=[];%'/home/lcg13/Data/pnrawDe/EpilepsyRelease07';
modality=[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%11:12;%15;%[];%15;%[];%[];%[];%[21 25];%[];%[9 8 7];%15:16;%[];%[];%[];%6;%8;%6;%9:100;%[];%[];%16:100;%[];
%series=13:14;%Second DWI neonatal%30:31;%16:17;%[];%9:40;%[];%3:100;%[];%[16 17];%[];%14:15;%[];%16:17;%[];%[];%[];%11:30;%[];%[];%3:6;%11;%[];%14:100;%[];%17:19;%17;%12;%[];%9-SB/6-MB4%[];%[10 12];%[2 6 7];%[10 12];%[];%7:20;%[];%3:6;%7:10;%3:6;%7:10;%8:9;%9;%[];%9;%[];%9%3:6;%7:10;%3:6;%7:10;%[];%7:10;%3:6;%[];%10%8;%[4 8];%[];%6:20;%[];%:20;%[];%[];%2;%[];%:100;%[];%:100;%16:100;%[];%[];%18:19;%[];%7:20;%2;%11:20;%:20;%:20;%[];%[];%17;%12;%22:100;%17:100;%8:20;%[];%10;%[];%[];%4:6;%10:12;%[];%12:100;%[];%:100;%[];%[];%6;%[];%5;%1:4;%[];%3;%4;%3:4;%[];%4;%[];%3:4;%:4;%[];%[];%[];%[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)
%series=10:11;%Second fMRI Neonatal
%series=31;%:31;%30:31;%DWI fetal
%series=16:17;%16:17;%16:17;%fMRI fetal
%series=14:15;%Third DWI neonatal
%series=11:12;%Third fMRI neonatal
%series=15;%14:15;%14:15;%14:15;%14:15;%14:15;%14:15;%14:15;%14:15;%First DWI neonatal
%series=10:11;%10:11;%10:11;%10:11;%First fMRI neonatal
specific='';%'Spikes';%'Detuning';

repeated=0;
caseFetal=0;
shaRootPath=outRootPath;
%shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Neonatal';
%if caseFetal;shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Fetal';end

%Lists of cases processed (from them, look for uncommented file in Studies to see the study that will be processed)
%abdominalData
%adultBrainData
%infantBrainData
%neonatalBrainData
%phantomData
%fetalBrainData
%cardiacData;

[path,pathid]=studiesDetection(1+caseFetal,[],outRootPath,path);

%path=path(slavaCases);%path=path(10:end);
%pathid=pathid(slavaCases);%pathid=pathid(10:end)

path=path(~cellfun('isempty',path));%To remove empty values from the path files

%%%NUMERO 178---DWI---2017_08_17/MA_245100

%for n=1:length(path);parseLabFolder(path{n},[],outRootPath,inpRootPath);end%PARSING RAW REQUIRES RECONFRAME LICENSE---WORKS IN perinatal130-pc AND gpubeastie03-pc
%for n=1:length(path);fRec=protoPipeline(path{n},modality,outRootPath,inpRootPath,series,specific,repeated);end%RECONSTRUCTING


%for n=1:length(path);fPro=procePipeline(path{n},modality,outRootPath,series);end%PREPROCESSING
%for l=0:2
%    for n=1:length(path);fPro=procePipelineTransform(path{n},modality,outRootPath,series,'Aq',l);end%PREPROCESSING
%end
%for n=1:length(path);fPro=procePipelineTransform(path{n},modality,outRootPath,series,'Aq');end

%for n=1:length(path);[fPro,~]=procePipelineSVR(path{n},5,outRootPath,series);end
%for n=1:length(path);reorientBrain(path{n},5,outRootPath);end

%for n=1:length(path);[fPro,~]=procePipelineAssembleSVR(path{n},modality,outRootPath,series);end


%for n=1:length(path);fileNameConversion(path{n},pathid{n},modality,outRootPath,inpRootPath,shaRootPath,series,0);end
%for n=1:length(path);fileNameConversion(path{n},pathid{n},modality,outRootPath,inpRootPath,shaRootPath,series,1);end
%sTInfo=[];
%for n=1:length(path);sTInfo=scanTimes(sTInfo,path{n},pathid{n},modality,outRootPath,inpRootPath,shaRootPath,series,1);end
%sTInfo=sTInfo(:,[2 3 6 7 5 4 1]);
%sTInfoTable=table(sTInfo);
%writetable(sTInfoTable,'dHCPScanningInfoNeonates.xls');

%for n=1:length(path);fAss=assessmentPipeline(path{n},modality,outRootPath,series);end
