%RUNS THE RECONSTRUCTION AFTER ROTATIONAL AND TRANSLATIONAL TRACKING

addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'));
clearvars
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVRTesting';
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR';
pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug07';
fetalBrainDataAlt;
pathf1=pathf(~cellfun('isempty',pathf));%To remove empty values from the path files

%pathf1=[];
%pathf1{1}='2019_05_21/TA_100930';
%pathf1{1}='2019_05_24/FO_101830';
pathf1{1}='2019_06_19/HO_107030';

for p=1:length(pathf1)
    pathSurro=sprintf('%s/%s/An-S2/SurrVoluRecoFirstHalf',pathOu,pathf1{p});
    filO=sprintf('%s/svr2.mat',pathSurro);
    %if isempty(dir(filO))
        fil=sprintf('%s/svr1.mat',pathSurro);  
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
                
        %svr.EnViews=logical([1 0 0 0 0 0 0 0 0]);
        svr.ParSVR.Step=svr.ParSVR.Step+3000;%Steps of the algorithm to write information
        %svr.ParSVR.convL=100;%1%10%100%1000;%TO ACCELERATE CONVERGENCE    
                
        if svr.ParSVR.Nl>=1
            fprintf('Solve SVR %s\n',path);tsta=tic;     
            svr=svrAlternateMinimization(svr,svr.ParSVR.maxNIt(1));
            tend=toc(tsta);fprintf('Time solving SVR: %.3f s\n\n',tend);    
        end

        %tic;save(sprintf('%s/svr2.mat',svr.ParSVR.pathSurro),'-v7.3');toc
    %end
end
