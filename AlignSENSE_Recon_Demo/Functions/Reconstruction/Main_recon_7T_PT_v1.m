
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));%cd /home/ybr19/Projects/Reconstruction/

%% SET PATHS
% Remove directories from path to remove unwanted repos loaded into MATLAB path
warning('off');
rmpath(genpath('/home/ybr19/Data')); 
rmpath(genpath('/home/ybr19/Software')); 
rmpath(genpath('/home/ybr19/Projects')); 
warning('on');

% Add data and repositories to MATLAB path
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
addpath(genpath('/home/ybr19/Software/Utilities'))
addpath(genpath('/home/ybr19/Projects/PTSignal'))
addpath(genpath('/home/ybr19/Projects/Reconstruction/'))

%% SELECT STUDIES
% List with all the names of the studies
if ~exist('study','var'); study = 17;end
if study == 1
    studies_20210302_GeomTesting_Phantom; addpath(genpath('/home/ybr19/Data/2021_03_02_GeomTesting_Phantom'))
    idPath = 1; idFile = [1]; 
    idRef = 1*ones(size(idFile)) ; 
elseif study == 3
    studies_20201002_DISORDER_MotionFree; addpath(genpath('/home/ybr19/Data/2020_10_02_DISORDER_MotionFree'))
    idPath= 1;  idFile= {1}; %9 is the hybrid k-space
    idRef = 4*ones(size(idFile)); % Take sensitivities from 2nd low-res scan
elseif study == 4
    studies_20201116_DISORDER_LowMotion;addpath(genpath('/home/ybr19/Data/2020_11_16_DISORDER_LowMotion'))
    idPath= 1;    idFile=[1,2,3];
    idRef = 1*ones(1,length(idFile)); % Take sensitivities from 1st low-res scan
elseif study == 5
    studies_20211012_PT_Phantom; addpath(genpath('/home/ybr19/Data/2021_10_12_PT_Phantom'))
    idPath= 1; idFile = 1:6;
    idRef = 1*ones(1,length(idFile)); % Take sensitivities from 1st low-res scan
elseif study==6
    studies_20210507_DISORDER_Motion_Jonny; addpath(genpath('/home/ybr19/Data/2021_05_07_DISORDER_Motion_Jonny'))
    idPath= 1; idFile = [5];
    idRef = 1*ones(1,length(idFile)); 
elseif study==7
    studies_20210811_Bandwidth_MoCo; addpath(genpath('/home/ybr19/Data/2021_08_11_Bandwidth_MoCo'))
    idPath= 1; idFile = {[10]};%[1,2,4,5,6,7,10];
    idRef = 4*ones(1,length(idFile)); 
elseif study==8
    studies_20210923_DISORDER_MotionFree_B0B1; addpath(genpath('/home/ybr19/Data/2021_08_11_DISORDER_MotionFree_B0B1'))
    idPath= 1; idFile = 1:9;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 9
    studies_20211025_MoCo_HighRes_PT; addpath(genpath('/home/ybr19/Data/2021_10_25_MoCo_HighRes_PT'))
    idPath=1; idFile = 2;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 10
    studies_20211029_MF_MoCo_PT; addpath(genpath('/home/ybr19/Data/2021_10_29_MF_MoCo_PT'))
    idPath= 1; idFile={1,6};%[6:9 11:12];
    idRef = 12*ones(1,length(idFile)); %This is the reference scan without PT
elseif study== 11
    studies_20211104_MoFreeHiResT1w; addpath(genpath('/home/ybr19/Data/2021_11_04_MoFreeHiResT1w'))
    idPath= 1; idFile = {[1]};
    idRef = 1*ones(1,length(idFile)); 
elseif study== 12
    studies_20211104_MoFreeHiResT1w_SJM; addpath(genpath('/home/ybr19/Data/2021_11_04_MoFreeHiResT1w_SJM'))
    idPath= 1; idFile = 2;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 13
    studies_20220329_DISORDER_Vida_Test; addpath(genpath('/home/ybr19/Data/2022_03_29_DISORDER_Vida_Test'))
    studies_20220405_VidaTestBody;addpath(genpath('/home/ybr19/Data/2022_04_05_VidaTestBody'))
    idPath= 1; idFile = 1;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 14
    studies_20220407_Vida_DISORDER_PT;addpath(genpath('/home/ybr19/Data/2022_04_07_Vida_DISORDER_PT'))
    idPath= 1; idFile = {4,9,10,11,12,[9 10]};
    idRef = 1*ones(1,length(idFile)); 
elseif study== 15
    studies_20220411_HiRes_PT_ShotExp;addpath(genpath('/home/ybr19/Data/2022_04_11_HiRes_PT_ShotExp/'))
    idPath= 1; idFile = 1:10;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 16
    studies_20220419_PT_Test_Osc;addpath(genpath('/home/ybr19/Data/2022_04_19_PT_Test_Osc/'))
    idPath= 1; idFile = 1:22;
    idRef = 1*ones(1,length(idFile)); 
elseif study== 17
    studies_20211014_PT_Phantom;addpath(genpath('/home/ybr19/Data/2021_10_14_PT_Phantom/'))
    idPath= 1; idFile = 1;
    idRef = 1*ones(1,length(idFile)); 
end
% Only consider files selected to reconstruct
%addpath(genpath(pathData{1}));
[pathIn, fileIn, refIn, refBIn, B0In, B1In, fileUnique, refUnique, refBUnique, B0Unique, B1Unique, isPTUnique] = extractStudies (pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath,idFile, idRef,[],[],isPT);

%% SET GLOBAL PARAMETERS
noiseFile=[];
%noiseFile = strcat(pathIn{1},filesep,refIn{1}{1})
supportReadout = [];
supportReadoutRecon = [];
resRec = [];
resRecRecon=[];

%%% Set parameters
estSensExplicit = 0;
estB1Explicit=0;
estB0Explicit=0;
invDataExplicit = 0;
runReconExplicit = 1;

%% PROCESS DATA 
for p=1:length(pathIn)
    
    %BUILD THE ACQUISITION DATA - start with this if an acquisition file is also used as a reference, you don't read it ouy twice (dat2Se checks if rec structure exist whereas dat2Rec reads data out by default)
    for f=1:length(fileUnique{p})
        fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
        if (~exist(strcat(fileName,'.mat'),'file') || invDataExplicit) && ~strcmp(fileName,'')
            dat2Rec(fileName,supportReadout,[],1,[],isPTUnique{p}{f},[],noiseFile, resRec);
        end                
    end
    
    %BUILD THE REFERENCE DATA - SENSITIVITY MAPS
    for f=1:length(refUnique{p})
        fileRef=strcat(pathIn{p},filesep,refUnique{p}{f});
        fileRefB=strcat(pathIn{p},filesep,refBUnique{p}{f});
        if ( ~exist(strcat(fileRef,'.mat'), 'file') || ~ismember('recS',who('-file', strcat(fileRef,'.mat')) ) ) || estSensExplicit && ~strcmp(fileRef,'')
            dat2Se(fileRef,fileRefB,[],1,1,[],study==50);
        end
    end
    
    %BUILD THE REFERENCE DATA - B1 MAPS
    for f=1:length(B1Unique{p})
        fileB1=strcat(pathIn{p},filesep,B1Unique{p}{f});
        fileRef=strcat(pathIn{p},filesep,refIn{p}{f});%Not unique!Assume reference in that scan is appropriate for B1 mapping
        if (~exist(strcat(fileB1,'.mat'), 'file')  || estB1Explicit) && ~strcmp(B1In{p}{f},'')
            dat2B1(fileRef,fileB1, [], 2, []);
        end
    end
    
    %BUILD THE REFERENCE DATA - B0 MAPS
    for f=1:length(B0Unique{p})
        fileB0=strcat(pathIn{p},filesep,B0Unique{p}{f});
        fileRef=strcat(pathIn{p},filesep,refIn{p}{f});%Assume reference in that scan is appropriate for B1 mapping
        if (~exist(strcat(fileB0,'.mat'), 'file') || estB0Explicit) && ~strcmp(B0In{p}{f},'')
            dat2B0(fileRef,fileB0, [], 2, []);
        end
    end
end


%% RUN RECONSTRUCTIONS
for p=1:length(pathIn)
    assert( length(idFile) == length(fileIn{p}),'File handling not correct.');
    
    for f=1:length(idFile)
        %% PREPARE RECONSTRUCTION STRUCTURE
        %%% ADD PATH TO ALL THE FILENAMES
        if iscell(fileIn{p}{f})
            for n=1:length(fileIn{p}{f});fileIn{p}{f}{n}=strcat(pathIn{p},filesep,fileIn{p}{f}{n});end
        else
            fileIn{p}{f} = strcat(pathIn{p},filesep,fileIn{p}{f});
        end
        
        fileRef=strcat(pathIn{p},filesep,refIn{p}{f});
        fileB1=strcat(pathIn{p},filesep,B1In{p}{f});
        fileB0=strcat(pathIn{p},filesep,B0In{p}{f});
        fileAcq=fileIn{p}{f};%Is allowed to be a cell array
        
        %%% MERGE DATA FROM ALL STRUCTURES
        writeSWCCFlag=0;
        rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon, supportReadoutRecon);
        
        %%% FILTER IF SLAB DETECTION
        [filterSize , isSlab] = filterForSlab(rec.Enc.FOVmm);
        if isSlab
            NY = size(rec.y);
            rec.y=bsxfun(@times,rec.y,ifftshift(buildFilter(NY(1:3),'tukey',filterSize,[],.4))); %Apodize y in image domain
            warning('Slab filtering activated. Make sure this is only done when estimating motion using this slab and not applying motion parameters.')
        end

        %%% GENERAL RECONSTRUCTION PARAMETERS
        rec.Alg.WriteSnapshots=0;       
        rec.Alg.UseSoftMasking= 0;

        Batch = 1; %YB: see solveXTB
        rec.Par.Mine.Modal = 7; %YB: Anatomical Volumetric Encoding

        rec.Alg.parXT.writeInter=0;

        rec.Dyn.GPUbS = [6 7];%[2 4]
        rec.Dyn.MaxMem = [6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
        
        %% MOTION CORRECTION
        rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
        %Disabling B0 estimation
        rec.Alg.parXB.useSH = 0;
        rec.Alg.parXB.useTaylor = 0;
        rec.Alg.parXB.useSusc = 0;
        
        %Motion estimation parameters
        rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
            %rec.Alg.parXT.resolMax = 8;
        rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
        rec.Alg.parXT.redFOV=0.3;
        rec.Alg.parXT.disableGrouping = 1;
        rec.Alg.parXT.convTransformJoint = 1;
        rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
        rec.Alg.parXT.meanT=0;
        
        %Optimisation scheme
        rec.Alg.resPyr =  [ 1   1];
        rec.Alg.nExtern = [ 15   1];
        rec.Alg.parXT.estT = [ 1  0 ];
        if ~rec.Enc.DISORDER %Need to specify the shot definition as it will not be detected automatically.
            shotDuration = 2;%In seconds
            numShots = [];
            if ~isempty(shotDuration);rec.Alg.parXT.sampleToGroup = shotDuration/(1e-3*rec.Par.Labels.RepetitionTime);end
            if ~isempty(numShots);rec.Alg.parXT.sampleToGroup = round(length(rec.Assign.z{2})/numShots);end
        end
        
        %redFac = 0;
        %[recT] = removeKSpaceSamples(rec, redFac);
        %recT.Alg.parXT.UseGiRingi=0;
        %solveXTB_PT(recT);  
        
        %Run        
        if ~exist(sprintf('%s%03d.mat',strcat(generateNIIFileName(rec),rec.Plan.Suff),Batch),'file') || runReconExplicit  
            solveXTB_PT(rec);         
        end
        
        %% MOTION + B0 CORRECTION (BASIS FUNCTIONS)
        rec.Plan.Suff='_MotB0Corr';

        rec.Alg.parXB.useSH = 0;%0 if disabled, 1 if in head frame, 2 if in scanner frame
        rec.Alg.parXB.BasisType = 'SH';% 'SHOW' / 'BSplines' - not yet implemented ar
        rec.Alg.parXB.SHorder = 2; % disable by setting < 0
        rec.Plan.SuffOu= sprintf('_SH%d', rec.Alg.parXB.SHorder);
        rec.Alg.parXB.nGD_B = 3;
        rec.Alg.parXB.deTaylor = [2 inf];%0 if de-activated / 1 if move cr to D / 2 to make D estimate with SH and discard cr / 3 if estimated 1 SH term and no motion dependence
        rec.Alg.parXB.Optim.basisDelay=0;  
        
        %Run
        if ~exist(sprintf('%s%03d.mat',strcat(generateNIIFileName(rec),rec.Plan.Suff),Batch),'file') || runReconExplicit  
            %solveXTB_PT(rec);         
        end
        
        %% MOTION + B0 CORRECTION (TAYLOR MODEL)
        rec.Alg.parXB.useTaylor = 1;
        rec.Plan.Suff='_MotB0Corr';rec.Plan.SuffOu= sprintf('_Taylor1');
        
        %General parameters
        rec.Alg.parXB.nGD_Taylor = 2;
        rec.Alg.parXB.weightedEst = 1;
        rec.Alg.parXB.modelShim = 0;
        
        %Optimisation
        rec.Alg.parXB.Optim.taylorDelay=6;
        rec.Alg.parXB.Optim.alphaList=10.^(-1* [-12:0.5:6]);%10.^(-1* [-8:0.5:6])
        
        %Contraining
        rec.Alg.parXB.C.filterTaylor.sp = 2.5;%In mm
        rec.Alg.parXB.C.filterTaylor.gibbsRinging = 0.1;
        rec.Alg.parXB.C.unwrapTaylor=1;
        
        %Run
        if ~exist(sprintf('%s%03d.mat',strcat(generateNIIFileName(rec),rec.Plan.Suff),Batch),'file') || runReconExplicit  
            %solveXTB_PT(rec);         
        end
        
        %% MOTION CORRECTION (PILOT TONE)
        if study==10; rec.Alg.parXB.useTaylor=0; else;rec.Alg.parXB.useTaylor = 0; end

        % ------ TEMP
        rec.Alg.WriteSnapshots=1;
        rec.Alg.parXT.groupSweeps = 1;
        rec.Alg.AlignedRec=1;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
        % ------ TEMP
        
        %Optimisation scheme
        rec.Alg.resPyr               = [   .5   1 ];
        rec.Alg.parXT.PT.usePT       = [     0   2 ]; %Don't use PT (0), use PT for alternating estimation (1), use PT to estimate motion and run one image est(2)
        rec.Alg.parXT.PT.subDiv      = [     1   0 ];
        rec.Alg.parXT.estT           = [     1   0 ];
        rec.Alg.nExtern              = [    15   2 ];
        rec.Alg.parXT.PT.nItActivate = [     0   0 ];

        %PT filtering
        rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthms = 0;%In ms, e.g. for TR=10, 40 samples
        rec.Alg.parXT.PT.PTFilter = []; 

        %PT signal to use
        rec.Alg.parXT.PT.signalUsage.useRealImag=1;
        rec.Alg.parXT.PT.signalUsage.whiten=1;
        rec.Alg.parXT.PT.signalUsage.whitenPhase=1;
        rec.Alg.parXT.PT.signalUsage.normaliseCoils=1;
        rec.Alg.parXT.PT.signalUsage.orderPreProcessing = [3 1];%in which order to apply the before. [1 3] good if useRealImag and [2 3] for useComplex
        rec.Alg.parXT.PT.signalUsage.eigTh = 0.03;%Threshold for eigenvalue thresholding
        rec.Alg.parXT.PT.signalUsage.useRealImagSVD=1;%Whether to apply SVD on real/imag parts instead of complex signal

        %PT signal usage for weighted recon
        rec.Alg.parXT.PT.weightedRecon.Flag=0;%Whether to use weights from PT instead of residual-based weights in last image reconstruction (1) or in every iteration(2)
        rec.Alg.parXT.PT.weightedRecon.Method='outlWePT';% =====  NOT GOOD DETECTION YET ====

        %PT calibration
        rec.Alg.parXT.PT.Calibration.implicitCalibration=1;
        rec.Alg.parXT.PT.Calibration.implicitCalibrationDelay=0;%Use LMsolver during this delay and do explicit calibration
        rec.Alg.parXT.PT.Calibration.implicitCalibrationNIt=1;
        rec.Alg.parXT.PT.Calibration.calibrationModel='backward';
        rec.Alg.parXT.PT.Calibration.calibrationOffset=0;
        rec.Alg.parXT.PT.Calibration.useRASMotion=0;%Whether to use the RAS motion parameters (1) and possibly the logarithm of the RAS rotation matrix (2) as this might linearise better
        rec.Alg.parXT.PT.Calibration.weightedCalibrationFlag=0;%To have weighted LS for PT calibration
        rec.Alg.parXT.PT.Calibration.weightedCalibrationMethod='outlWePT';%whether to use soft weights from PT trace (WePT), hard weights (outlWePT) or discrepancy between PT and DISORDER (WeDisPT --> To be implemented!!)

        %PT-based clustering
        rec.Alg.parXT.PT.profileClustering.FlagFinal=0;
        rec.Alg.parXT.PT.profileClustering.FlagInitialise=0;%Whether to use PT clustering to infer initial clustering
        rec.Alg.parXT.PT.profileClustering.clusterMotion=1;
        rec.Alg.parXT.PT.profileClustering.numClusters=300;
        rec.Alg.parXT.PT.profileClustering.traLimX=[];%Threshold in within-cluster-distance for determining numClusters (over-writes the parameters numClusters)

        %Run
        if ~exist(sprintf('%s%03d.mat',strcat(generateNIIFileName(rec),rec.Plan.Suff),Batch),'file') || runReconExplicit  
            %solveXTB_PT(rec);         
        end
        
        rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu= sprintf('_imp_eigTh10_NoClust');
        rec.Alg.parXT.PT.Calibration.implicitCalibration=1;
        rec.Alg.parXT.PT.signalUsage.eigTh = -10;%Threshold 
        rec.Alg.parXT.PT.profileClustering.FlagFinal=0;
        solveXTB_PT(rec);  
        
        rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu= sprintf('_imp_eigTh10_Clust');
        rec.Alg.parXT.PT.profileClustering.FlagFinal=1;
        solveXTB_PT(rec);  
        
    end
end
