
clc; close all;
clear all
cd /home/ybr19/Projects/Reconstruction/
addpath(genpath('/home/ybr19/Projects/Reconstruction/'))

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
addpath(genpath('/home/ybr19/Software/Utilities'))
addpath(genpath('/home/ybr19/Projects/PTSignal'))

warning('off');rmpath(genpath('/home/ybr19/Data')); warning('on');

% List with all the names of the studies
study = 10;
if study ==1
    studies_20210302_GeomTesting_Phantom; addpath(genpath('/home/ybr19/Data/2021_03_02_GeomTesting_Phantom'))
    idPath = 1; idFile = [1]; 
    idRef = 1*ones(size(idFile)) ; 
elseif study ==2
    studies_20191023_DISORDER_Lucilio; addpath(genpath('/home/ybr19/Data/2019_10_23_DISORDER_Lucilio'))
    idPath = 5; idFile = [1]; % 4 is the alternating continuous motion (fig 5 in DISORDER7T.pdf)
    idRef = 1*ones(size(idFile)) ; 
elseif study ==3
    studies_20201002_DISORDER_MotionFree; addpath(genpath('/home/ybr19/Data/2020_10_02_DISORDER_MotionFree'))
    idPath= 1;  idFile= {1}; % 9 is the hybrid k-space
    idRef = 4*ones(size(idFile)); % Take sensitivities from 2nd low-res scan
elseif study ==4
    studies_20201116_DISORDER_LowMotion;addpath(genpath('/home/ybr19/Data/2020_11_16_DISORDER_LowMotion'))
    idPath= 1;    idFile=[1,2,3];
    idRef = 1*ones(1,length(idFile)); % Take sensitivities from 1st low-res scan
elseif study==5
    studies_20201214_PilotTone; addpath(genpath('/home/ybr19/Data/2020_12_14_PilotTone'))
    idPath= 1; idFile = [1];
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
elseif study==9
    studies_20211025_MoCo_HighRes_PT; addpath(genpath('/home/ybr19/Data/2021_10_25_MoCo_HighRes_PT'))
    idPath= 1; idFile = 2;
    idRef = 1*ones(1,length(idFile)); 
elseif study==10
    studies_20211029_MF_MoCo_PT; addpath(genpath('/home/ybr19/Data/2021_10_29_MF_MoCo_PT'))
    idPath= 1; idFile = {6};%[6:9 11:12];
    idRef = 12*ones(1,length(idFile)); %This is the reference scan without PT
elseif study==11
   studies_20211104_MoFreeHiResT1w; addpath(genpath('/home/ybr19/Data/2021_11_04_MoFreeHiResT1w'))
    idPath= 1; idFile = {[1]};
    idRef = 1*ones(1,length(idFile)); 
elseif study==12
   studies_20211104_MoFreeHiResT1w_SJM; addpath(genpath('/home/ybr19/Data/2021_11_04_MoFreeHiResT1w_SJM'))
    idPath= 1; idFile = 2;
    idRef = 1*ones(1,length(idFile)); 
end

idB1=1*ones(1,length(idFile));
idB0=1*ones(1,length(idFile));
[pathIn, fileIn, refIn, B0In, B1In, fileUnique, refUnique, B0Unique, B1Unique] = extractStudies(pathIn, fileIn, refIn, B0In, B1In, idPath,idFile, idRef,idB0,idB1);

%fileIn{1}{1} = strcat(fileIn{1}{1},'_Slab')

noiseFile = strcat(pathIn{1},filesep,refIn{1}{1});  
noiseFile=[]
supportReadout = [];%[120 160]/248
supportReadoutRecon = [];
resRec = [];
resRecRecon=[];
isPTList = [5 9:12 ];

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
            dat2Rec(fileName,supportReadout,[],1,[],ismember(study,isPTList),[],noiseFile, resRec);
        end                
    end
    
    %BUILD THE REFERENCE DATA - SENSITIVITY MAPS
    for f=1:length(refUnique{p})
        fileRef=strcat(pathIn{p},filesep,refUnique{p}{f});
        if ( ~exist(strcat(fileRef,'.mat'), 'file') || ~ismember('recS',who('-file', strcat(fileRef,'.mat')) ) ) || estSensExplicit && ~strcmp(fileRef,'')
            dat2Se(fileRef,[],1,1,[],study==50);
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
        
        %%% MERGE DATA FRO ALL STRUCTURES
        writeSWCCFlag=0;
        rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon, supportReadoutRecon);
        
        %%% RECONSTRUCT
        rec.Alg.parXT.accel=[1 0];%[1 0]%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
        rec.Alg.WriteSnapshots=0;
        
        %rec.Alg.parXT.fractionOrder=0.5;
        rec.Alg.parXT.perc=[0.95 0.95 0.95];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
        
        Batch = 1; %YB: see solveXTB
        rec.Par.Mine.Modal = 7; %YB: Anatomical Volumetric Encoding
                
        rec.Plan.SuffOu='DephCorrSH4';
        rec.Plan.Suff='DephCorrSH4';
        fileName=generateNIIFileName(rec);
        fileName=strcat(fileName,rec.Plan.Suff);           

        if ~exist(sprintf('%s%03d.mat',fileName,Batch),'file') || runReconExplicit  % The file would contain raw reconstruction data of last levels.
            if isfield(rec , 'Plan.Types'); rec.Plan=rmfield(rec.Plan,'Types');end
            if isfield(rec , 'Plan.TypeNames'); rec.Plan=rmfield(rec.Plan,'TypeNames');end          
            if isfield(rec , 'Dyn.Typ2Wri'); rec.Dyn.Typ2Wri(end+1)=0;end

            rec.Alg.parXT.computeCSRecon = 0; %YB: Note that if this is activated, the typ2Rec(18) will automatically be set to 1
            %parXT.corrFact = 0.6;

            rec.Alg.parXT.exploreMemory = 0;
            rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
            rec.Alg.parXT.writeInter=0;
            rec.Alg.parXT.meanT=0;
            rec.Alg.UseSoftMasking= 0;

            rec.Dyn.GPUbS = [6 7];%[2 4]
            rec.Dyn.MaxMem=[6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
            rec.Alg.parXT.maximumDynamics = 6;
            
            %%% MotCorr parameters
            rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
            rec.Alg.parXT.resolMax = 8;
            rec.Alg.parXT.saveFinal=0;
            rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
            rec.Alg.parXT.NWend = 1; 
            rec.Alg.parXT.groupSweeps=2;% Factor to group sweeps
            rec.Alg.parXT.redFOV=0.3;
            rec.Alg.parXB.useSH = 0;
            rec.Alg.parXB.useTaylor = 0;
            rec.Alg.parXB.useSusc = 0;
rec.Alg.parXT.disableGrouping = 1;
rec.Alg.resPyr = [0.25 0.5 1];
rec.Alg.nExtern = [40 30 15];
rec.Alg.parXT.convTransformJoint = 1;
            rec.Alg.parXT.percRobustShot=[0.25 0.75];%Percentiles for robust computation of expected inter-shot dispersion
            rec.Alg.parXT.enerRobustShot=0.9;%0.95;%Ratio of the error for acceptance of shots
            rec.Alg.parXT.corrFact = 1;
            
            %%% Basis functions parameters - DISABLED
            rec.Plan.Suff='_DephCorr';

            rec.Alg.parXB.useSH = 0;%0 if disabled, 1 if in head frame, 2 if in scanner frame
            rec.Alg.parXB.BasisType = 'SH';% 'SHOW' / 'BSplines' - not yet implemented ar
            rec.Alg.parXB.SHorder = 2; % disable by setting < 0
            rec.Plan.SuffOu= sprintf('_SH%d', rec.Alg.parXB.SHorder);
            rec.Alg.parXB.nGD_B = 3;
            rec.Alg.parXB.deTaylor = [2 inf];%0 if de-activated / 1 if move cr to D / 2 to make D estimate with SH and discard cr / 3 if estimated 1 SH term and no motion dependence
            rec.Alg.parXB.Optim.basisDelay=0;  
            %solveXTB_ext(rec);  
            
            %% With Taylor model
            rec.Alg.parXB.useTaylor = 1;
            rec.Alg.parXB.orderTaylor = 1;
            rec.Alg.parXB.useSusc = 0;
            if rec.Alg.parXB.useSH == 1
                rec.Plan.Suff='_DephCorr';rec.Plan.SuffOu= sprintf('_SH%dTaylor%d', rec.Alg.parXB.SHorder,rec.Alg.parXB.orderTaylor);
            elseif rec.Alg.parXB.useSH == 0
                rec.Plan.Suff='_DephCorr';rec.Plan.SuffOu= sprintf('_Taylor%d',rec.Alg.parXB.orderTaylor);
            end
            
            rec.Alg.parXB.Reg.lambda_sparse = 0e+0; %1e+5; % Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_smooth = 0e+2; % 3e+5 Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_susc = 0e+3;
            rec.Alg.parXB.nGD_Taylor = 2;
            rec.Alg.parXB.filterD.sp = 2.5;%In mm
            rec.Alg.parXB.filterD.gibbsRinging = 0.6;
            rec.Alg.parXB = rmfield(rec.Alg.parXB,'filterD');
            rec.Alg.parXB.deSH = 0;
            rec.Alg.parXB.Optim.FISTA=1;
            rec.Alg.parXB.Optim.RMSProp=0; rec.Alg.parXB.Optim.beta=0.9;
            rec.Alg.parXB.Optim.DupLim=inf;
            rec.Alg.parXB.Optim.taylorDelay=6;
            rec.Alg.parXB.Optim.alphaList=10.^(-1* [-8:0.5:6]);
            rec.Alg.parXB.unwrapTaylor=1;
            rec.Alg.parXB.denoiseTaylor=0;
            
            rec.Alg.parXB.weightedEst = 1;
            rec.Alg.parXB.intraVoxDeph = 0;
            rec.Alg.parXB.modelShim = 0;
            
            %solveXTB(rec);
            close all
            
            %% Pilot Tone
            if iscell(fileAcq) && length(fileAcq)>1; rec.Names.Name = 'Combined'; end
            
            if study==10; rec.Alg.parXB.useTaylor = 0;else;rec.Alg.parXB.useTaylor = 0; end
            
            rec.Alg.parXB.weightedEst = 0;
            rec.Alg.parXT.redFOV=0.35

if study==9 && contains(fileAcq,'Slab') 
    rec.Alg.parXT.redFOV=0
    NY=size(rec.y);
    rec.y=bsxfun(@times,rec.y,ifftshift(buildFilter(NY(1:3),'tukey',[1 20 20],[],.4))); %Apodize % YB: y in image domain
end
            rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
            rec.Alg.WriteSnapshots=1;
            
            rec.Alg.parXT.groupSweeps = 1;
            rec.Alg.resPyr          = [.5  1];
            rec.Alg.parXT.PT.usePT  = [0   1]; % Don't use PT (0), use PT for alternating estimation (1), use PT at end to predict next step (2)
            rec.Alg.parXT.PT.subDiv = [0   0 ];
            rec.Alg.parXT.PT.estT   = [1   0];
            rec.Alg.nExtern = [ 10 10  10  10];
            rec.Alg.parXT.PT.nItActivate = [5 0 0 0 ];
                        
            rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthms = 0;%In ms, so for TR=10, 40 samples
            rec.Alg.parXT.PT.PTFilter = []; 
            
            rec.Alg.parXT.PT.signalUsage.useRealImag=1;
            rec.Alg.parXT.PT.signalUsage.whiten=1;
            rec.Alg.parXT.PT.signalUsage.whitenPhase=1;
            rec.Alg.parXT.PT.signalUsage.normaliseCoils=1;
            rec.Alg.parXT.PT.signalUsage.orderPreProcessing = [1 3];%in which order to apply the before. [1 3] good if useRealImag and [2 3] for useComplex
            
            rec.Alg.parXT.PT.weightedRecon.Flag=0;%Whether to use weights from PT instead of residual-based weights in last image reconstruction (1) or in every iteration(2)
            rec.Alg.parXT.PT.weightedRecon.Method='outlWePT';% =====  NOT GOOD DETECTION YET ====
            
            rec.Alg.parXT.PT.Calibration.implicitCalibration=1;
            rec.Alg.parXT.PT.Calibration.calibrationModel='backward';
            rec.Alg.parXT.PT.Calibration.calibrationOffset=0;
            rec.Alg.parXT.PT.Calibration.useRASMotion=0;%Whether to use the RAS motion parameters (1) and possibly the logarithm of the RAS rotation matrix (2) as this might linearise better
            rec.Alg.parXT.PT.Calibration.weightedCalibrationFlag=0;%To have weighted LS for PT calibration
            rec.Alg.parXT.PT.Calibration.weightedCalibrationMethod='outlWePT';%whether to use soft weights from PT trace (WePT), hard weights (outlWePT) or discrepancy between PT and DISORDER (WeDisPT --> To be implemented!!)
            
            rec.Alg.parXT.PT.profileClustering.Flag=1;
            rec.Alg.parXT.PT.profileClustering.clusterMotion=0;
            rec.Alg.parXT.PT.profileClustering.numClusters=200;
            rec.Alg.parXT.PT.profileClustering.traLimX=0.01;%Threshold in within-cluster-distance for determining numClusters (over-writes the parameters numClusters)
            
            rec.Plan.Suff='_PT_1mm_1tr';rec.Plan.SuffOu= sprintf('');
            solveXTB_PT(rec);
            
%             rec.Alg.parXT.PT.subDiv = [1   1  1  0];
%             rec.Plan.Suff='_PT_1mm_p5tr';rec.Plan.SuffOu= sprintf('');
%             solveXTB_PT(rec);
%              
%             rec.Alg.parXT.PT.subDiv = [1   1  0  0] %No point of subdividing if not using PT
%             rec.Alg.parXT.PT.usePT  = [0  0  0 0];
%             rec.Plan.Suff='_DIS_1mm_1tr';
%             solveXTB_PT(rec);
            
        end
    end
end



