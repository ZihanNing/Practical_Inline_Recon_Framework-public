
proces=0;

pcur=fileparts(mfilename('fullpath'));
addpath(genpath(sprintf('%s/..',pcur)));

[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

caseFetal=0;

paths=studiesDetection(1+caseFetal);

inpRootPath=[];%'/home/lcg13/Data/18ExperimentsTracking';%[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/Detuning';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
if caseFetal;outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/DetuningFetal';end
modality=10;%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%[];%[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)

NP=length(paths);
percReads=50*ones(2,NP);
percSlices=50*ones(2,NP);
percVolumes=50*ones(2,NP);
coils=zeros(32,NP);
casesEnd=0;
detuningFile=sprintf('%s/detuning.mat',outRootPath);
if exist(detuningFile,'file')
    load(detuningFile);
    NC=size(percReads,2);
    if NC<NP;casesEnd=1;end
    percReads(:,end+1:NP)=50;
    percSlices(:,end+1:NP)=50;
    percVolumes(:,end+1:NP)=50;
    coils(:,end+1:NP)=0;
end
nOr=find(percReads(2,:)==50 & percSlices(2,:)==50 & percVolumes(2,:)==50,1,'last');
if isempty(nOr);nOr=NP;end
if proces==1 && nOr>0
    if casesEnd;nEn=NC+1;else nEn=1;end   
    for n=nOr:-1:nEn%Neonates%Checked from 781 to 722
    %for n=829
        fprintf('\n\nCase %d: %s',n,paths{n});
        fRec=protoPipeline(paths{n},modality,outRootPath,inpRootPath,series,'Detuning');  
        %return
        percReads(:,n)=0;
        percSlices(:,n)=0;
        percVolumes(:,n)=0;
        for l=1:length(fRec)
            if contains(fRec{l}.Names.Name,'mbdti2')
                percReads(:,n)=percReads(:,n)+fRec{l}.Anomaly.PercReads(:,2);
                percSlices(:,n)=percSlices(:,n)+fRec{l}.Anomaly.PercSlices(:,2);
                percVolumes(:,n)=percVolumes(:,n)+fRec{l}.Anomaly.PercVolumes(:,2);
                coils(fRec{l}.Anomaly.Coils,n)=1;
            end
        end
        if percReads(2,n)==50 && percSlices(2,n)==50 && percVolumes(2,n)==50;percReads(:,n)=0;percSlices(:,n)=0;percVolumes(:,n)=0;end
        
        save(detuningFile,'percReads','percSlices','percVolumes','coils');
    end
else
    for s=1:length(paths)
        dat(s)=datetime(paths{s}(1:10),'InputFormat','yyyy_MM_dd');
        if s==1;dat=repmat(dat,[1 length(paths)]);end
    end
    nOr=find(percReads(2,:)==50 & percSlices(2,:)==50 & percVolumes(2,:)==50,1,'last');
    nOr=nOr+1;
    if nOr<830
        typ={'reads','slices','volumes'};
        for n=1:3
            datRed=dat(nOr:end);
            coilsRed=coils(:,nOr:end);
            if n==1;percRed=percReads(1,nOr:end)./percReads(2,nOr:end);    
            elseif n==2;percRed=percSlices(1,nOr:end)./percSlices(2,nOr:end);    
            elseif n==3;percRed=percVolumes(1,nOr:end)./percVolumes(2,nOr:end);
            end
            datRed(isnan(percRed))=[];
            coilsRed(:,isnan(percRed))=[];
            percRed(isnan(percRed))=[];
            coilsRed(:,datRed<datetime('2018_09_01','InputFormat','yyyy_MM_dd'))=[];
            percRed(datRed<datetime('2018_09_01','InputFormat','yyyy_MM_dd'))=[];
            datRed(datRed<datetime('2018_09_01','InputFormat','yyyy_MM_dd'))=[];
            
            percRed=percRed*100;

            figure(n)
            plot(datRed,percRed,'*')
            xlabel('date')
            ylabel(sprintf('Percentage of damaged %s',typ{n}))    
            hold on
            for n=1:length(percRed)
                if percRed(n)>0
                    plot([datRed(n) datRed(n)],[0 percRed(n)],'b');
                end
            end
            %ylim([0 100])
            grid on
        end 
        coilsRedSum=sum(coilsRed,1);
        datRed=repelem(datRed,coilsRedSum);
        coilsRedIn=find(coilsRed==1);
        coilsRedIn=ind2subV(size(coilsRed),coilsRedIn);        
        figure(4)
        plot(datRed,coilsRedIn(:,1),'*')
        xlabel('date')
        ylabel('Damaged coils')
        hold on
        grid on
        
        %coils=sum(coils,2);
        %fprintf('Damaged coils%s\n',sprintf('-%d',find(coils~=0)'));
    end
end