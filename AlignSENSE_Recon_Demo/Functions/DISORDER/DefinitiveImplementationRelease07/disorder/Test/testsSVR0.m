%RUNS THE TRANSLATIONAL AND ROTATIONAL TRACKING

addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'));
clearvars
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVRTesting';
%pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR';
pathOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug07';
fetalBrainDataAlt;
pathf0=pathf(~cellfun('isempty',pathf));%To remove empty values from the path files

%pathf0=[];
%pathf0{1}='2018_11_15/DO_62030';
pathf0{1}='2019_06_19/HO_107030';


for p=1:length(pathf0)
    pathSurro=sprintf('%s/%s/An-S2/SurrVoluReco',pathOu,pathf0{p});    
    filO=sprintf('%s/svr1.mat',pathSurro);   
    if isempty(dir(filO))
        fil=sprintf('%s/svr0.mat',pathSurro);
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

        svr.ParSVR.Step=0;%Steps of the algorithm to write information

        svr=svrExcitationStructures(svr);%Slice orders, number of packages, number of slices per package and indexes of slices in package
        svr=svrRearrangeAxes(svr);
        
        %svr.ParSVR.FWHM=2;%FULL WIDTH HALF MAXIMUM (RATIO OF THE SLICE THICKNESS VS SLICE SEPARATION)

        svr.ParSVR.Prerun=1;
        svrPrev=svr;%To store default values
        svr.ParSVR.FOVSize=320;%FOV of the reconstruction
        if modal==5;svr.ParSVR.MS=2;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
        fprintf('Setting up SVR %s\n',path);tsta=tic;     
        svr=svrSetUp(svr);
        tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

        svrPrev.ParSVR.Step=svr.ParSVR.Step;svrPrev.ParSVR.Prerun=0;
        svr=svrPrev;
        svr.ParSVR.FOVSize=320;%FOV of the reconstruction
        if modal==5;svr.ParSVR.MS=2;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
        fprintf('Setting up SVR %s\n',path);tsta=tic;     
        svr=svrSetUp(svr);
        tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

        %TRACKING
        fprintf('Track SVR %s\n',path);tsta=tic;     
        svr=svrTracking(svr,0);
        tend=toc(tsta);fprintf('Time tracking SVR: %.3f s\n\n',tend);

        svrPrev.El=svr.El;svrPrev.Elp=svr.Elp;svrPrev.ParSVR.Step=svr.ParSVR.Step;%Update relevant fields
        svrPrev.ParSVR.OverEncode=1;%Modify reset structure
        if svr.ParSVR.Nl>=1
            svr=svrPrev;%Reset svr structure      
            svr.ParSVR.FOVSize=mean(svr.MS(:))*median(max(svr.Elp(1,4:6,:),[],2),3)*svr.ParSVR.multFOV;%It was 3 before
            svr.ParSVR.FOVSize=16*svr.ParSVR.Resol*ceil(svr.ParSVR.FOVSize/(16*svr.ParSVR.Resol)); 
            fprintf('Automatic FOV size: %.2f\n',svr.ParSVR.FOVSize); 
            %svr.ParSVR.FOVSize=144;%FOV of the reconstructiono

            %%%MODIFIED
            svr.ParSVR.fracOrd=0;%0.25;

            
            if modal==5;svr.ParSVR.MS=2;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
            fprintf('Setting up SVR %s\n',path);tsta=tic;
            svr=svrSetUp(svr);
            tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

            fprintf('Track SVR %s\n',path);tsta=tic;     
            svr=svrTracking(svr,1);
            tend=toc(tsta);fprintf('Time tracking SVR: %.3f s\n\n',tend);                    
        end

        tic;save(sprintf('%s/svr1.mat',svr.ParSVR.pathSurro),'-v7.3');toc
        
        if svr.ParSVR.Nl>=1
            fprintf('Solve SVR %s\n',path);tsta=tic;     
            svr=svrAlternateMinimization(svr,svr.ParSVR.maxNIt(1));
            tend=toc(tsta);fprintf('Time solving SVR: %.3f s\n\n',tend);    
        end

        tic;save(sprintf('%s/svr2.mat',svr.ParSVR.pathSurro),'-v7.3');toc
        
    end
end
