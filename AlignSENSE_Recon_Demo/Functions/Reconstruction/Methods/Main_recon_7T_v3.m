
clc; close all;
clear all
cd /home/ybr19/Reconstruction/
addpath(genpath('.'))

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
addpath(genpath('/home/ybr19/dB0_analysis/PhaseUnwrapping'))
rmpath(genpath('/home/ybr19/Data'))

% List with all the names of the studies
study = 4;

if study ==2
    studies_YB_testrecon;
    addpath(genpath('/home/ybr19/Data/YB_testrecon'))
    % Select studies to reconstruction
    id_path = 5;
    id_file = [1:8]; % 4 is the alternating continuous motion (fig 5 in DISORDER7T.pdf)
    id_file = {[2,3],[1,2,8]};
    id_ref = 1*ones(size(id_file)) ; 
elseif study ==3
    studies_20201002_DISORDER_MotionFree;
    addpath(genpath('/home/ybr19/Data/2020-10-02_DISORDER_MotionFree'))
    % Select studies to reconstruction
    id_path= 1;  id_file= [1]; % 9 is the hybrid k-space
    %id_file = {[2:6 ],[1,6,8],[3,5,7],[2,3,5]};
    id_ref = 1*ones(size(id_file)); % Take sensitivities from 2nd low-res scan
    disp('Taking the sensitivity estimation of the second scan - assuming coil loading negligible.')
elseif study ==4
    studies_20201116_DISORDER_LowMotion;
    addpath(genpath('/home/ybr19/Data/2020-11-16_DISORDER_LowMotion'))
    id_path= 1;    
    id_file=[ 3 ];
    %id_file= {[4],[5],[4,5]};  
    id_ref = 2*ones(1,length(id_file)); % Take sensitivities from 1st low-res scan
elseif study==5
    studies_20201214_PilotTone;
    addpath(genpath('/home/ybr19/Data/2020-12-14_PilotTone'))
    id_path= 1; id_file= {[1,2]};   
    id_file = [1];
    id_ref = 1*ones(1,length(id_file)); % Take sensitivities from 1st low-res scan
end

[pathIn, refIn, fileIn] = extract_studies(pathIn, refIn, fileIn, id_path, id_ref, id_file);

%%% Set parameters
estSens_explicit = 0;
invData_explicit = 0;
runRecon_explicit = 1;

gibbsRinging=1;
corrMotion=1;
corrPhase=1;

supportReadS=[];supportRead=[]; % change depending on scan - especially if not acquired yourself

for p=1:length(pathIn)
    
    %BUILD THE REFERENCE DATA - COIL ESTIMATIONS
    for f=1:length(refIn{p})
        fileName=strcat(pathIn{p},filesep,refIn{p}{f});
        if ~exist(strcat(fileName,'.mat'), 'file') || estSens_explicit 
            TWS=mapVBVD(strcat(fileName,'.dat'));            
            if study==5; recS=invert7T_PilotTone(TWS,supportReadS);else recS=invert7T(TWS,supportReadS);end
            
            rmpath(genpath('/home/ybr19/dB0_analysis/PhaseUnwrapping'))
            recS=solveSensit7T(recS);
            addpath(genpath('/home/ybr19/dB0_analysis/PhaseUnwrapping'))

            fileNameS = strcat(pathIn{p}, 'Re-Se/',refIn{p}{f});if ~exist( strcat(pathIn{p}, 'Re-Se/'),'dir');mkdir(strcat(pathIn{p}, 'Re-Se/'));end
            xW=[];xW{1} = recS.S; MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
            writeNII(fileNameS, {'Se'},xW, MSW, MTW);
            
            save(strcat(fileName,'.mat'),'recS','-v7.3');    
%%%TO DO: estimateSensitivities (fileName, pathOu, writeNII);
%%%TO DO: dat2recSiemens(fileName, pathOu);
%%%TO DO: assassignSensToRec(fileName, fileref);
        end
    end
   
    %RECONSTRUCT
    for f=1:length(id_file)
        rec=[];fileName=[];
        if iscell(id_file)
            for n=1:length(id_file{f});fileName{n}=strcat(pathIn{p},filesep,fileIn{p}{id_file{f}(n)});end
        else
            fileName=strcat(pathIn{p},filesep,fileIn{p}{f});
        end
        if ~iscell(fileName);fileNametemp = fileName;else fileNametemp = fileName{1};end%latter in case of multiple files
        if ~exist(strcat(fileNametemp,'.mat'),'file') || invData_explicit                      
            refName=strcat(pathIn{p},filesep,refIn{p}{f});           
            TW=mapVBVD(strcat(fileNametemp,'.dat'));
            if study==5; rec=invert7T_PilotTone(TW,supportRead);else rec=invert7T(TW,supportRead);end
            load(strcat(refName,'.mat'));
            NS=size(recS.S);NS=NS(1:3);
            NY=size(rec.y);NY=NY(1:3);
            NSY=max(round(NS.*recS.Enc.AcqVoxelSize./rec.Enc.AcqVoxelSize),NY(1:3));

            if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(2*NS,'tukeyIso',0.5,0,gibbsRinging,1);elseif gibbsRinging==-1;HY=buildFilter(2*NS,'CubicBSpline',[],0,1);else HY=[];end
            if ~isempty(HY)
                rec.S=filtering(recS.S,HY,1);
                rec.W=abs(filtering(recS.W,HY,1));
                HY=buildFilter(NS,'tukeyIso',1,0,gibbsRinging);
                rec.M=abs(filtering(abs(recS.x),HY)/20);           
            else            
                rec.S=recS.S;rec.W=recS.W;rec.M=recS.x;
            end
            rec.S=resampling(rec.S,NSY);rec.W=resampling(rec.W,NSY);rec.M=resampling(rec.M,NSY); 
            rec.S=resampling(rec.S,NY,2);rec.W=resampling(rec.W,NY,2);rec.M=resampling(rec.M,NY,2);
            if isfield(rec,'N');rec.S=standardizeCoils(rec.S,rec.N);end
            rec=rmfield(rec,'W');
            
            NY=size(rec.y);NY(end+1:4)=1;
            
            offset = -1;
            Lines = mod( TW.image.Lin -1 +offset, NY(2) ) +1;
            Partitions = mod( TW.image.Par -1 +offset, NY(3) ) +1;
            
            rec.Assign.z{2}= (NY(2) - (Lines-1))-floor((NY(2))/2) -1; % 2nd PE  %Minus comes from convention in invert7T.m
            rec.Assign.z{3}= (NY(3) - (Partitions-1))-floor((NY(3))/2)-1; % 3rd PE = slices
            
% rec.Assign.z{2}=TW.image.Lin-floor(NY(2)/2)-1; %YB: 2nd PE
% rec.Assign.z{3}=TW.image.Par-floor(NY(3)/2)-1; %YB: 3rd PE - slices

            [~,rec.Names.Name] = fileparts(fileNametemp);
            rec.Names.pathOu=pathIn{p};

            rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
            rec.Par.Labels.ZReconLength=1;
            
            rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
            rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
            rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
            rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;
                                   
            save(strcat(fileNametemp,'.mat'),'rec','-v7.3');    
            fprintf('File saved:\n   %s', fileNametemp);
        else           
            if ~iscell(fileName)
                load(strcat(fileName,'.mat'));              
            else
                for n=1:length(fileName)
                    load(strcat(fileName{n},'.mat'));                    
                    if n==1
                       y=rec.y;
                       z{2}=rec.Assign.z{2};
                       z{3}=rec.Assign.z{3};
                       Name=rec.Names.Name;

                    else
                       y=cat(5,y,rec.y);
                       z{2}=cat(2,z{2},rec.Assign.z{2});
                       z{3}=cat(3,z{3},rec.Assign.z{3});
                       Name=strcat(Name,rec.Names.Name);
                    end
                end
                rec.y=y;y=[];
                rec.Assign.z{2}=z{2};z{2}=[];
                rec.Assign.z{3}=z{3};z{3}=[];
                rec.Names.Name=Name;
            end
        end 

        %rec.Alg.parXT.discardHighRes=0.4;%Those volumetric scans whose resolution is below this value (in mm) are not reconstructed
        rec.Alg.parXT.accel=[1 0];%[1 0]%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
        rec.Alg.WriteSnapshots=0;
        
        %rec.Alg.parXT.fractionOrder=0.5;
        rec.Alg.parXT.perc=[0.95 0.95 0.95];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
        rec.Dyn.Log=1;
        
        Batch = 1; %YB: see solveXTB
        rec.Par.Mine.Modal = 7; %YB: Anatomical Volumetric Encoding
                
        rec.Plan.SuffOu='DephCorrSH4';
        rec.Plan.Suff='DephCorrSH4';
        fileName=generateNIIFileName(rec);
        fileName=strcat(fileName,rec.Plan.Suff);           

        if ~exist(sprintf('%s%03d.mat',fileName,Batch),'file') || runRecon_explicit  % The file would contain raw reconstruction data of last levels.
            if isfield(rec , 'Plan.Types'); rec.Plan=rmfield(rec.Plan,'Types');end
            if isfield(rec , 'Plan.TypeNames'); rec.Plan=rmfield(rec.Plan,'TypeNames');end          
            if isfield(rec , 'Dyn.Typ2Wri'); rec.Dyn.Typ2Wri(end+1)=0;end

            rec.Alg.parXT.computeCSRecon = 0; %YB: Note that if this is activated, the typ2Rec(18) will automatically be set to 1
            parXT.corrFact = 0.6;

            rec.Alg.parXT.exploreMemory = 0;
            rec.Alg.parXT.traLimXT=[0.009 0.007];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
            rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
            rec.Alg.parXT.writeInter=0;
            rec.Alg.parXT.meanT=0;
            rec.Alg.UseSoftMasking= 0;

            rec.Dyn.GPUbS = [6 7];%[2 4]
            rec.Dyn.MaxMem=[6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
            rec.Alg.parXT.maximumDynamics = 6;
            close all
            
            %% MotCorr
            rec.Plan.Suff='_MotCorr';
            rec.Alg.parXT.resolMax = 5;
            rec.Alg.parXT.saveFinal=0;
            rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
            rec.Alg.parXT.NWend = 1; 
            rec.Alg.parXT.groupSweeps=2;% Factor to group sweeps
            rec.Alg.parXT.redFOV=0.28;
            rec.Plan.SuffOu='';
            rec.Alg.parXB.dephaseCorrection=0;
            rec.Alg.parXB.useSH = 0;
            rec.Alg.parXB.useTaylor = 0;
            rec.Alg.parXB.useSusc = 0;
rec.Alg.parXT.disableGrouping = 1;
rec.Alg.nExtern = 20;
rec.Alg.parXT.convTransformJoint = 1;
            rec.Alg.parXT.percRobustShot=[0.25 0.75];%Percentiles for robust computation of expected inter-shot dispersion
            rec.Alg.parXT.enerRobustShot=0.9;%0.95;%Ratio of the error for acceptance of shots
            rec.Alg.parXT.corrFact = 1;

            %recD=solveXTB_ext(rec);  
            
            %% Basis functions
            rec.Plan.Suff='_DephCorr';

            rec.Alg.parXB.useSH = 2;%1 if in head frame, 2 if in scanner frame
            rec.Alg.parXB.BasisType = 'SH';% 'SHOW' / 'BSplines' - not yet implemented ar
            rec.Alg.parXB.SHorder = 2; % disable by setting < 0
            rec.Plan.SuffOu= sprintf('_SH%d', rec.Alg.parXB.SHorder);
            rec.Alg.parXB.nGD_B = 3;
            rec.Alg.parXB.deTaylor = [2 inf];
            rec.Alg.parXB.Optim.basisDelay=0;  
            %solveXTB_ext(rec);  
            
            rec.Alg.parXB.useB1 = 0;%%%B1 correction
            rec.Alg.parXB.B1order = 3; 
            rec.Alg.parXB.Optim.B1Delay=15;
            
            %% With Taylor model
            close all
            rec.Alg.parXB.useTaylor = 1;
            rec.Alg.parXB.orderTaylor = 1;
            rec.Alg.parXB.useSusc = 0;
            rec.Plan.Suff='_DephCorr';
            rec.Plan.SuffOu= sprintf('_SH%dTaylor%d', rec.Alg.parXB.SHorder,rec.Alg.parXB.orderTaylor);
          
            rec.Alg.parXB.Reg.lambda_sparse = 0e+0; %1e+5; % Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_smooth = 0e+2; % 3e+5 Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_susc = 0e+3;
            rec.Alg.parXB.nGD_Taylor = 2;
            rec.Alg.parXB.filterD.sp = 2.5;%In mm
            rec.Alg.parXB.filterD.gibbsRinging = 0.6;
            rec.Alg.parXB = rmfield(rec.Alg.parXB,'filterD')
            rec.Alg.parXB.deSH = 0;
            rec.Alg.parXB.Optim.FISTA=1;
            rec.Alg.parXB.Optim.RMSProp=0; rec.Alg.parXB.Optim.beta=0.9;
            rec.Alg.parXB.Optim.DupLim=inf;
            rec.Alg.parXB.Optim.taylorDelay=4;
            rec.Alg.parXB.Optim.alphaList=10.^(-1* [-3:0.7:6]);
            rec.Alg.parXB.unwrapTaylor=1;
            rec.Alg.parXB.denoiseTaylor=1;
            
            rec.Alg.parXB.weightedEst = 1;
            rec.Alg.parXB.redFOV = 1/(2.7);
            rec.Alg.parXB.intraVoxDeph = 0;
           
            solveXTB(rec);
        end

    end
end



