%RUNS THE RECONSTRUCTION AFTER ALIGNMENT AT HIGH RESOLUTION (THAT IS, ONLY FINAL RECONSTRUCTIONS)

addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'));
clearvars
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVRTesting';
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR';
pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug07';
fetalBrainDataAlt;
pathf3=pathf(~cellfun('isempty',pathf));%To remove empty values from the path files

%pathf3=[];
%pathf3{1}='2018_11_15/DO_62030';
pathf3{1}='2019_06_19/HO_107030';

for p=1:length(pathf3)
    %pathSurro=sprintf('%s/%s/An-S2/SecondHalfFW2-20x',pathOu,pathf3{p});
    pathSurro=sprintf('%s/%s/An-S2/FirstHalf',pathOu,pathf3{p});
    filO=sprintf('%s/%s/An-S2/*_ExAAA0.80.nii',pathOu,pathf3{p});
    %if isempty(dir(filO))
        fil=sprintf('%s/svr3.mat',pathSurro);
        if ~isempty(dir(fil));load(fil);else continue;end

        svr.ParSVR.removeData=0;    
        if svr.ParSVR.removeData
            warning('off','MATLAB:DELETE:FileNotFound');
            ste=svr.ParSVR.Step+1;
            while 1
                fil=sprintf('%s/*Step%02d*.nii',svr.ParSVR.pathSurro,ste);
                if ~isempty(dir(fil));delete(fil);ste=ste+1;else break;end
            end
            warning('on','all');
        end               

        svr.ParSVR.Step=svr.ParSVR.Step+5000;%Steps of the algorithm to write information
        if svr.ParSVR.Nl>=3
            svrPrev.TP=svr.TP;svrPrev.TV=svr.TV;svrPrev.TE=svr.TE;svrPrev.Mx=svr.Mx;svrPrev.M=svr.M;svrPrev.ParSVR.Step=svr.ParSVR.Step;svrPrev.EnViews=svr.EnViews;
            svrPrev.x=svr.x;svrPrev.W=svr.W;svrPrev.We=svr.We;            
            svrPrev.ParSVR.tiSh=svr.ParSVR.tiSh;%SHEARLET REGULARIZATION
            %svrPrev.ParSVR.ti(2)=50;
            svr=svrPrev;

            svr.ParSVR.FOVSize=mean(svr.MS(:))*median(max(svr.Elp(1,4:6,:),[],2),3)*svr.ParSVR.multFOV;
            svr.ParSVR.FOVSize=16*svr.ParSVR.Resol*ceil(svr.ParSVR.FOVSize/(16*svr.ParSVR.Resol)); 
            svr.ParSVR.MS=1;
            if modal==5;svr.ParSVR.MS=1;else svr.ParSVR.MS=1.25;end
            %svr.ParSVR.Resol=0.8;

            svr.ParSVR.MS=svr.ParSVR.MS*svr.ParSVR.Resol;

            fprintf('Setting up SVR %s\n',path);tsta=tic;   
            svr=svrSetUp(svr);
            tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

            %ACCURATE RECONSTRUCTION
            nItFinal=7;
            fprintf('Final reconstruction %s\n',path);tsta=tic;    
            svr=svrCG(svr,nItFinal,[],[],2*nItFinal+1);
            tend=toc(tsta);fprintf('Time final reconstruction: %.3f s\n\n',tend);

            %WE WRITE THE DATA            
            svrWriteData(svr,'Ma',svr.Mx,[],0.8,[],[],1);
            svrWriteData(svr,'Ex',svr.x,[],0.8,[],[],[],0.2);%WE WRITE THE DATA CROPPING TO POSITIVE AND AT 0.8
        end        
    %end
end
