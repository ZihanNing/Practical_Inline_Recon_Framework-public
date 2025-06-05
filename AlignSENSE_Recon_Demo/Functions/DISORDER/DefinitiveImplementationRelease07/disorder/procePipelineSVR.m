function [fPro,svr]=procePipelineSVR(path,modal,rootOu,series,specific)

%PROCEPIPELINE   Runs a preprocessing pipeline for a given modality
%   [FPRO,SVR]=PROCEPIPELINE(PATH,{MODAL},{ROOTOU},{SERIES},{SPECIFIC})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%   * {SERIES} restricts the reconstructions to a specific set of series
%   * {SPECIFIC} indicates to use a specific configuration of parameters as
%   stated in reconSpecific.m
%   ** FPRO returns the rec structures where the method failed
%   ** SVR returns the SVR structure
%

addpath(genpath(fileparts(mfilename('fullpath'))));

%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode,~,~,~,~,pathNt]=versRecCode;

if nargin<2;modal=5;end%WE USE A SINGLE MODALITY, THE T1 CODE IS DEEMED AS EXPERIMENTAL
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end

%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');error('Protocol file %s not found',protFile);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

fPro=[];
svr=[];
csvr=1;
cinh=1;
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin>=4 && ~isempty(series);nV=nV(ismember(nV,series));end
if nargin<5;specific=[];end

for n=nV
    rec.Names.Name=strtrim(prot.B_FileName(n,:));
    matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));
    rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;rec.Names.Specific=specific;
    if exist(matFile,'file')
        load(matFile);
        rec.Par=Par;Par=[];rec.Par.Mine.Proce=1;
        rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);
        if ~rec.Fail
            rec.Dyn.Typ2Wri(:)=0;modal=rec.Par.Mine.Modal;
            if modal==5 || modal==6          
                fileRECON=fullfile(rec.Names.pathOu,numbe2Modal(modal),rec.Names.Name,'');
                if modal==5
                    if ~isempty(regexp(fileRECON,'_sf')) || ~isempty(regexp(fileRECON,'_fd')) || ~isempty(regexp(fileRECON,'_f\d'));express='t2mbz';noexpress='';%dHCP project
                    %elseif ~isempty(regexp(fileRECON,'_dev')) || ~isempty(regexp(fileRECON,'_pip5'));express='pack';noexpress='';%PIP project or dev patch
                    elseif ~isempty(regexp(fileRECON,'_dev')) || ~isempty(regexp(fileRECON,'_pip5'));express='t2';noexpress='pack';%PIP project or dev patch, to reconstruct standard scans
                    else express='.';noexpress='';
                    end
                elseif modal==6;express='t1';noexpress='';
                end
                if ~isempty(regexp(fileRECON,express)) && isempty(regexp(fileRECON,noexpress))                     
                    suff=[];
                    suff{1}='Aq';
                    cont=2;
                    existFile=1;
                    for e=1:length(suff)
                        suff{e}=strcat(suff{e},rec.Plan.Suff);
                        if ~exist(strcat(fileRECON,'_',suff{e},'.nii'),'file');existFile=0;break;end
                    end                                        
                    if existFile    
                        %if csvr==1
                            for l=n-1:-1:1
                                recB1.Names.Name=strtrim(prot.B_FileName(l,:));
                                matFileB1=fullfile(headFolder,sprintf('%s.mat',recB1.Names.Name));
                                recB1.Names.matFile=matFileB1;recB1.Names.pathOu=pathOu;recB1.Names.prot=prot;recB1.Names.ind=n;recB1.Names.headFolder=headFolder;recB1.Names.versCode=versCode;                                
                                if exist(matFile,'file')
                                    load(matFileB1);
                                    recB1.Par=Par;Par=[];recB1.Par.Mine.Proce=1;
                                    recB1.Fail=0;recB1=reconPlanner(recB1);recB1=reconAlgorithm(recB1);
                                    if ~recB1.Fail
                                        recB1.Dyn.Typ2Wri(:)=0;modal=recB1.Par.Mine.Modal;
                                        if modal==4
                                            fileRECONB1=fullfile(recB1.Names.pathOu,numbe2Modal(modal),recB1.Names.Name,'');
                                            existFileB1=1;
                                            for e=1:length(suff)
                                                suff{e}=strcat(suff{e},recB1.Plan.Suff);
                                                if ~exist(strcat(fileRECONB1,'_',suff{e},'.nii'),'file');existFileB1=0;break;end
                                            end
                                            if existFileB1
                                                fprintf('Reading B1 %s\n',recB1.Names.matFile);tsta=tic;
                                                [yB1,MSB1,MTB1]=readNII(fileRECONB1,suff,0);
                                                svr.yB1{cinh}=yB1{1};
                                                svr.MSB1{cinh}=MSB1{1};
                                                svr.MTB1{cinh}=MTB1{1};
                                                svr.recB1{cinh}=recB1;
                                                cinh=cinh+1;
                                                tend=toc(tsta);fprintf('Time reading B1: %.3f s\n\n',tend);yB1=[];recB1=[];
                                                break;
                                            else
                                                fprintf('Reconstruction data file %s not found\n',strcat(fileRECONB1,'_',suff{1},'.nii'));
                                            end
                                        end
                                    end
                                end
                            end                            
                        %end
                        
                        fprintf('Reading %s\n',rec.Names.matFile);tsta=tic;
                        %NOISE
                        suff{2}='No';
                        suff{2}=strcat(suff{2},rec.Plan.Suff);
                        if ~exist(strcat(fileRECON,'_',suff{2},'.nii'),'file');suff=suff(1);end                        
                        [y,MS,MT]=readNII(fileRECON,suff,rec.Dyn.GPU);
                        N=size(y{1});N(end+1:4)=1;
                        fprintf('Series: %s\n',rec.Names.Name);
                        fprintf('Resolution:%s\n',sprintf(' %.2f',MS{1}));
                        fprintf('Acquired resolution:%s\n',sprintf(' %.2f',rec.Par.Scan.AcqVoxelSize));
                        fprintf('Slice gap:%s\n',sprintf(' %.2f',rec.Par.Scan.SliceGap(1)));
                        fprintf('Number of dynamics: %d\n',N(4));
                        fprintf('Scan duration:%s\n',sprintf(' %.2f',rec.Par.Labels.ScanDuration(1)));
                        for d=1:N(4)%Extract the averages
                            svr.y{csvr}=dynInd(y{1},d,4);
                            if length(y)==2;svr.no{csvr}=dynInd(y{2},2,4);end
                            svr.MS{csvr}=MS{1};
                            svr.MT{csvr}=MT{1};
                            svr.rec{csvr}=rec; 
                            if N(4)>1;svr.rec{csvr}.Names.Name=sprintf('%s_%d',svr.rec{csvr}.Names.Name,d);end
                            svr.B1Ind(csvr)=cinh-1;
                            svr.dynamic{csvr}=d;
                            csvr=csvr+1;
                        end       
                        
                        tend=toc(tsta);fprintf('Time reading: %.3f s\n\n',tend);y=[];
                    else
                        fprintf('Reconstruction data file %s not found\n',strcat(fileRECON,'_',suff{1},'.nii'));
                    end
                end
                %%if modal==10;break;end%The heuristic says that only the first one can be processed, as Anthony acquired the doubly reversed thing
            end    
        end
    end
    rec=[];
end

MSRes=cat(3,MS{:});
MSRes=min(MSRes,[],3);
%svr.ParSVR.FWHM=4;%*max(MSRes)/min(MSRes);%4;%FULL WIDTH HALF MAXIMUM (RATIO OF THE SLICE THICKNESS VS THE SLICE SEPARATION)
if strcmp(express,'t2mbz');svr.ParSVR.FWHM=4;else svr.ParSVR.FWHM=2;end%FULL WIDTH HALF MAXIMUM (RATIO OF THE SLICE THICKNESS VS THE SLICE SEPARATION)

%%%SLIGHTLY CHANGED
svr.ParSVR.Resol=min(MSRes);
%%%SLIGHTLY CHANGED
fprintf('Target resolution: %.2f\n',svr.ParSVR.Resol);

csvr=csvr-1;
if csvr==0;fprintf('NO EXISTING VOLUMES\n\n');return;end

modal=svr.rec{1}.Par.Mine.Modal;
%PARAMETERS
svr.ParSVR.UseApo=0.1;

svr.ParSVR.EstT=3;%3;

svr.ParSVR.convL=1;%1%10%100%1000;%TO ACCELERATE CONVERGENCE    
if modal==5;svr.ParSVR.fracOrd=0.125;else svr.ParSVR.fracOrd=1;end
%if modal==5;svr.ParSVR.fracOrd=0.25;else svr.ParSVR.fracOrd=1;end
%if modal==5;svr.ParSVR.fracOrd=0.5;else svr.ParSVR.fracOrd=1;end%0.25 shown to be more stable for small fetuses than 0.5

svr.ParSVR.targetResFact=1.1;%TARGET RESOLUTION AS A FACTOR OF THE GRID RESOLUTION, THE BIGGER THE LOWER RESOLUTION
%assert(svr.ParSVR.targetResFact>1,'The target resolution factor should be strictly higher than 1 and it is %.3f',svr.ParSVR.targetResFact);
%ordFracDer=cos(pi/(2*svr.ParSVR.targetResFact)).^(-2);
%fprintf('Fractional derivative order for target resolution %.2f: %.2f\n',svr.ParSVR.targetResFact,ordFracDer);

%svr.ParSVR.ti=[1 50];%[1 1 1];%TIKHONOV REGULARIZER
svr.ParSVR.path=path;
svr.ParSVR.OverEncode=1;
svr.ParSVR.Debug=0;
%svr.ParSVR.ti=[1 5];

%NEW---WE ONLY RELY ON THE SHEARLET TRANSFORM
svr.ParSVR.ti=[0 0];
%NEW-SHEARLET TRANSFORM-IT HAS TO BE SET TO SOMETHING DIFFERENT FROM 0 LATER ON
svr.ParSVR.tiSh=0;

svr.ParSVR.spti=[1 1/svr.ParSVR.targetResFact];%[1 1 1];%SPACING TIKHONOV REGULARIZER
svr.ParSVR.regFracOrd=[0 100];%[2 1 0];%REGULARIZATION FRACTIONAL ORDER
svr.ParSVR.Lp=0;%NORM OF THE RECONSTRUCTION
svr.ParSVR.Nl=3;%LEVELS OF CORRECTION
svr.ParSVR.maxNIt=[600 16];%NUMBER OF ITERATIONS OF CORRECTION
svr.ParSVR.GibbsRingi=0;%TO APPLY GIBBS RINGING BEFORE PROCESSING---BETTER TO INCORPORATE TO THE REGULARIZER
svr.ParSVR.CorrectInhomPowe=0.75;%1 before
if modal==5;svr.ParSVR.winic=1e-2;else svr.ParSVR.winic=1;end
svr.ParSVR.SlBefore=0;
svr.ParSVR.pathNt=pathNt;%Path to pretrained neural network
svr.ParSVR.Visualize=0;%To visualize intermediate results of SVR
svr.ParSVR.Step=0;%Steps of the algorithm to write information
svr.ParSVR.quickRecon=0;%For quick reconstructions
svr.ParSVR.usePrecond=0;
svr.ParSVR.pathModal=fullfile(svr.rec{1}.Names.pathOu,numbe2Modal(svr.rec{1}.Par.Mine.Modal));
suffAux='SurrVoluReco';%FWHM4';
svr.ParSVR.pathSurro=fullfile(svr.ParSVR.pathModal,suffAux);
svr.ParSVR.removeData=0;
svr.ParSVR.multFOV=3;%To avoid boundary issues without incurring in too much computational penalty
%if strcmp(express,'t2mbz') || strcmp(express,'pack');svr.ParSVR.maxNoViews=6;%Maximum number of views used for actual reconstruction
if strcmp(express,'t2mbz');svr.ParSVR.maxNoViews=6;%Maximum number of views used for actual reconstruction
else svr.ParSVR.maxNoViews=12;
end

if svr.ParSVR.removeData==1
    warning('off','MATLAB:DELETE:FileNotFound');
    delete(sprintf('%s/*',svr.ParSVR.pathSurro));
    delete(sprintf('%s/*_El.nii',svr.ParSVR.pathModal));
    delete(sprintf('%s/*_Re*nii',svr.ParSVR.pathModal));
    delete(sprintf('%s/*_Ex*.nii',svr.ParSVR.pathModal));
    delete(sprintf('%s/*_Ma*.nii',svr.ParSVR.pathModal));
    warning('on','all');
    svr.ParSVR.removeData=0;
elseif svr.ParSVR.removeData==2
    warning('off','MATLAB:DELETE:FileNotFound');
    delete(sprintf('%s/*.nii',svr.ParSVR.pathSurro));
    for n=2:3;delete(sprintf('%s/svr%d.mat',svr.ParSVR.pathSurro,n));end  
    delete(sprintf('%s/*_Ex*.nii',svr.ParSVR.pathModal));
    delete(sprintf('%s/*_Ma*.nii',svr.ParSVR.pathModal));
    warning('on','all');
    svr.ParSVR.removeData=0;
end

if svr.ParSVR.Visualize;svrVisualization(svr.y,[],'Original data');end

if svr.ParSVR.CorrectInhomPowe && isfield(svr,'yB1');svr=svrCorrectInhom(svr);end

svr=svrCNNBrainDetection(svr);

svr=svrFiltering(svr,2);%Filtering stacks

if ~exist(svr.ParSVR.pathSurro,'dir');mkdir(svr.ParSVR.pathSurro);end
tic;save(sprintf('%s/svr0.mat',svr.ParSVR.pathSurro),'-v7.3');toc

svr=svrExcitationStructures(svr);%Slice orders, number of packages, number of slices per package and indexes of slices in package
svr=svrRearrangeAxes(svr);

svr.ParSVR.Prerun=1;
svrPrev=svr;%To store default values
svr.ParSVR.FOVSize=320;%FOV of the reconstruction
if modal==5;svr.ParSVR.MS=2;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
%svr.ParSVR.MS=svr.ParSVR.MS*svr.ParSVR.Resol;
fprintf('Setting up SVR %s\n',path);tsta=tic;     
svr=svrSetUp(svr);
tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

svrPrev.ParSVR.Step=svr.ParSVR.Step;svrPrev.ParSVR.Prerun=0;
svr=svrPrev;
svr.ParSVR.FOVSize=320;%FOV of the reconstruction
if modal==5;svr.ParSVR.MS=2;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
%svr.ParSVR.MS=svr.ParSVR.MS*svr.ParSVR.Resol;
fprintf('Setting up SVR %s\n',path);tsta=tic;     
svr=svrSetUp(svr);
tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

%TRANSLATIONAL TRACKING
fprintf('Track SVR %s\n',path);tsta=tic;     
svr=svrTracking(svr,0);
tend=toc(tsta);fprintf('Time tracking SVR: %.3f s\n\n',tend);

svrPrev.El=svr.El;svrPrev.Elp=svr.Elp;svrPrev.ParSVR.Step=svr.ParSVR.Step;%Update relevant fields
svrPrev.ParSVR.OverEncode=1;%Modify reset structure
if svr.ParSVR.Nl>=1
    svr=svrPrev;%Reset svr structure      
    svr.ParSVR.tiSh=1;%SHEARLET REGULARIZATION
    svr.ParSVR.FOVSize=mean(svr.MS(:))*median(max(svr.Elp(1,4:6,:),[],2),3)*svr.ParSVR.multFOV;
    svr.ParSVR.FOVSize=16*svr.ParSVR.Resol*ceil(svr.ParSVR.FOVSize/(16*svr.ParSVR.Resol)); 
    fprintf('Automatic FOV size: %.2f\n',svr.ParSVR.FOVSize); 
    if modal==5;svr.ParSVR.MS=1.5;else svr.ParSVR.MS=2.5;end%Resolution of the reconstruction
    svr.ParSVR.MS=2;%svr.ParSVR.MS*svr.ParSVR.Resol;
    fprintf('Setting up SVR %s\n',path);tsta=tic;     
    svr=svrSetUp(svr);
    tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);
    
    %ROTATIONAL TRACKING
    fprintf('Track SVR %s\n',path);tsta=tic;     
    svr=svrTracking(svr,1);
    tend=toc(tsta);fprintf('Time tracking SVR: %.3f s\n\n',tend);    
    
    tic;save(sprintf('%s/svr1.mat',svr.ParSVR.pathSurro),'-v7.3');toc    

    fprintf('Solve SVR %s\n',path);tsta=tic;     
    svr=svrAlternateMinimization(svr,svr.ParSVR.maxNIt(1));
    tend=toc(tsta);fprintf('Time solving SVR: %.3f s\n\n',tend);    
    
    tic;save(sprintf('%s/svr2.mat',svr.ParSVR.pathSurro),'-v7.3');toc
end


if svr.ParSVR.Nl>=2    
    svrPrev.TP=svr.TP;svrPrev.TV=svr.TV;svrPrev.TE=svr.TE;svrPrev.ParSVR.Step=svr.ParSVR.Step;svrPrev.EnViews=svr.EnViews;%Update relevant fields
    svrPrev.x=svr.x;svrPrev.W=svr.W;
    svr=svrPrev;
    svr.ParSVR.tiSh=1;%SHEARLET REGULARIZATION
    svr.ParSVR.FOVSize=mean(svr.MS(:))*median(max(svr.Elp(1,4:6,:),[],2),3)*svr.ParSVR.multFOV;
    svr.ParSVR.FOVSize=16*svr.ParSVR.Resol*ceil(svr.ParSVR.FOVSize/(16*svr.ParSVR.Resol)); 
    if modal==5;svr.ParSVR.MS=1;else svr.ParSVR.MS=1.25;end
    %svr.ParSVR.Resol=0.8;
    
    svr.ParSVR.MS=svr.ParSVR.MS*svr.ParSVR.Resol;
    fprintf('Setting up SVR %s\n',path);tsta=tic;     
    svr=svrSetUp(svr);
    tend=toc(tsta);fprintf('Time setting up SVR: %.3f s\n\n',tend);

    fprintf('Solve SVR %s\n',path);tsta=tic;     
    svr=svrAlternateMinimization(svr,svr.ParSVR.maxNIt(2));
    tend=toc(tsta);fprintf('Time solving SVR: %.3f s\n\n',tend);
    
    tic;save(sprintf('%s/svr3.mat',svr.ParSVR.pathSurro),'-v7.3');toc
end

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
    svrWriteData(svr,'MaFW4',svr.Mx,[],0.8,[],[],1);
    svrWriteData(svr,'ExFW4',svr.x,[],0.8,[],[],[],0.2);%WE WRITE THE DATA CROPPING TO POSITIVE AND AT 0.8
end


% if csvr>=2%ROUNDING ISSUES IN DEFINITION OF FOV FOR A SINGLE VOLUME. A TESTING CASE TO IMPROVE ON THIS COULD BE 2018_03_28/VA_732
% else fprintf('Not enough %s studies have been detected for case %s\n',numbe2Modal(modal),path);
% end