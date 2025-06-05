
proces=0;

addpath(genpath(fileparts(mfilename('fullpath'))));

[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

paths=studiesDetection(2);

    inpRootPath=[];%'/home/lcg13/Data/18ExperimentsTracking';%[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
    outRootPath=[];%'/home/lcg13/Data/18ReplicaRelease06';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
    modality=5;%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
    series=[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)

    if proces==1    
        for n=1:length(paths)
            fRec=protoPipeline(paths{n},modality,outRootPath,inpRootPath,series);
            NV=length(fRec);     
            c=1;
            for v=1:NV
                if isfield(fRec{v},'Max')
                    Max{n}(c,:)=fRec{v}.Max;
                    c=c+1;
                end
            end    
        end    
        save Max Max
    else
        load Max
        for n=1:length(Max)
            if ~isempty(Max{n})
                figure(1)
                plot(n*ones(1,size(Max{n},1)),Max{n}(:,1),'b*')
                hold on
                figure(2)
                plot(n*ones(1,size(Max{n},1)),Max{n}(:,2),'r*')
                hold on  
                figure(3)
                plot(n*ones(1,size(Max{n},1)),Max{n}(:,2)./Max{n}(:,1),'g*')
                hold on
            end            
        end
        figure(1)
        xlabel('fHCP Case')
        ylabel('Max signal before PDA correction')
        grid on
    
        figure(2)
        xlabel('fHCP Case')
        ylabel('Max signal after PDA correction')
        grid on
        
        figure(3)
        xlabel('fHCP Case')
        ylabel('Ratio after/before PDA')
        grid on
    end
    
    




