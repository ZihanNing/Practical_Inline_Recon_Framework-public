
proces=1;

addpath(genpath(fileparts(mfilename('fullpath'))));

[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

caseFetal=0;

%IF RUNNING THIS WE SHOULD UNCOMMENT THE LINE THAT SAVES THE ANOMALIES

%False positives, THEY ARE ARGUABLE
%%%%%HEREHEREHERE
%1) Case 633/1-Risk 9.56 Perhaps /home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_08_18-SQ_245400-sq_18082017_1402015_2_2_dhcp8refneoheadV4.png (usual problem with reference scans where central lines show vertical halos, solved by discarding checks on these lines)
%2) Case 608/8-Risk 12.42/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_30-JE_223701-je_30052017_1416561_8_2_dhcp8sbfmriclearV4.png (probably due to background problems)
%3) Case 608/7-Risk 10.21/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_30-JE_223701-je_30052017_1433165_10_2_dhcp8sbfmrirepclearV4.png (same as before, probably due to background problems)
%4) Case 608/6-Risk INF/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_30-JE_223701-je_30052017_1414103_7_2_dhcp8sefmriclearV4.png (same as before, probably due to background problems)
%5) Case 605/7-Risk 10.69/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_23-MI_222000-mi_23052017_1028067_9_2_dhcp8sbfmriclearV4.png (same as before, probably due to background problems)
%6) Case 600/7-Risk INF/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_12-DE_219701-de_12052017_1456375_8_2_dhcp8t2w3ddrive32chshcsenseV4.png (probably due to quantifization problems)
%7) Case 600/6-Risk INF/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_12-DE_219701-de_12052017_1451393_7_2_dhcp8t2w3ddrive32chshcsenseV4.png (same as 6)
%8) Case 600/5-Risk INF/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Anomalies/Snapshots/2017_05_12-DE_219701-de_12052017_1444545_6_2_dhcp8t2w3ddrive32chshcsenseV4.png (same as 7)

[paths,pathid]=studiesDetection(1+caseFetal);

inpRootPath=[];%'/home/lcg13/Data/18ExperimentsTracking';%[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Neonatal';
outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/Reference';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
if caseFetal
    shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Fetal';
    outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/ReferenceFetal';
end
modality=[2:4 9];%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%[];%[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)

NCori=1;
NCend=length(paths);

%for n=NCori:NCend
%    fprintf('\n\nCase %d of %d: %s\n',n-NCori+1,NCend,paths{n});
%    parseLabFolder(paths{n},modality,outRootPath,inpRootPath);%PARSING RAW REQUIRES RECONFRAME LICENSE
%end

%return

for n=NCori:NCend%468%[1 425 NCend]%:1%length(paths)
    modality=2:4;%2:4;%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
    tsta=tic;
    fprintf('\n\nCase %d of %d: %s',n-NCori+1,NCend,paths{n});
    fRec=protoPipeline(paths{n},modality,outRootPath,inpRootPath,series);

    modality=[];  
    fileNameConversion(paths{n},pathid{n},modality,outRootPath,inpRootPath,shaRootPath,series);
    tend=toc(tsta);fprintf('Time: %.3f s\n',tend);
end