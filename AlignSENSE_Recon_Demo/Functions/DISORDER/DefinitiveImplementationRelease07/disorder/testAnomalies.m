
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

paths=studiesDetection(1+caseFetal);

inpRootPath=[];%'/home/lcg13/Data/18ExperimentsTracking';%[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/Anomalies';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
if caseFetal;outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/AnomaliesFetal';end
modality=[];%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%[];%[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)

%Case 15: 2015_03_10/GA_15104 anomalies not detected in series 17 (T1MS axial)---Series 14 really
%-Appear as lines in the readout direction, not clear why not picked up, perhaps because they grow linearly
%Case 8: 2015_02_24/RO_11002 anomalies not detected in series 15 (MB diffusion)---Series 12 really
%-Only 16 volumes detected, but there are 300?-Second MB not detected-check file with directions.
%-Regarding spikes whole spectrum damaged which may explain problems in detecting, perhaps with old method...
%Case 7: 2015_02_23/PA_10700 anomalies not detected in series 13 (MB diffusion)---Series 12 really
%-Regarding spikes, same as in case 8
%Case 1: 2015_02_06/HO_7201 anomalies not detected in series 15 (T1MS axial)---Series 13 really
%-Regarding spikes, same as in case 15, blocky appearance, massive

anomalies=0.5*ones(1,length(paths));
anomalyFile=sprintf('%s/anomalies.mat',outRootPath);
if exist(anomalyFile,'file');load(anomalyFile);end
if proces==1
    %for n=721:-1:1%Neonates%Checked from 781 to 722
    %for n=110:-1:1%Fetuses%Checked from 116 to 111
    %for n=length(paths):-1:1
    %%%Case 743-SB-Case 744-MB
    %for n=[625 550 512 466 448 432 426]% 781:-1:1%70%781:-1:761%52:-1:1%76:-1:1%75%781:-1:775%81:-1:770%1%[1 7 8 15]%781%768
    for n=781:-1:1%81%781:-1:1%762%762 vs 781%:-1:1%781:-1:1
        fprintf('\n\nCase %d: %s',n,paths{n});
        fRec=protoPipeline(paths{n},modality,outRootPath,inpRootPath,series);
        anomalies(n)=0;
        for l=1:length(fRec);anomalies(n)=anomalies(n)+fRec{l}.AnomalyDetected(1);end
        if anomalies(n)==0
            for l=1:length(fRec)
                if isfield(fRec{l},'maxProfRisk') && fRec{l}.maxProfRisk>6;anomalies(n)=0.5;break;end;
            end
        end
        %save(sprintf('%s/anomalies.mat',outRootPath),'anomalies');
    end    
else
    for s=1:length(paths)
        dat(s)=datetime(paths{s}(1:10),'InputFormat','yyyy_MM_dd');
        if s==1;dat=repmat(dat,[1 length(paths)]);end
    end
    %anomalies([772 699 686 674 659 645 626])=anomalies([772 699 686 674 659 645 626])-1;
    
    dat=dat(562:end);
    anomalies=anomalies(562:end);
    
    figure(1)
    plot(dat,anomalies,'*')
    xlabel('date')
    ylabel('Number of outliered series')    
    hold on
    for n=1:length(anomalies)
        if anomalies(n)>0
            plot([dat(n) dat(n)],[0 anomalies(n)],'b');
        end
    end
    grid on
%     figure(2)
%     plot(anomalies)
%     xlabel('dHCP Case')
%     ylabel('Number of outliered series')
%     grid on    
%     percProb=sum(anomalies>0)./(sum(anomalies==0)+sum(anomalies>0));
    fprintf('Corrupted cases: %.2f%%\n',100*percProb);
end