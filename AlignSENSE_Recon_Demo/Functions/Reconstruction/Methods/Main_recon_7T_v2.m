
clc; close all;
clear all
cd /home/ybr19/Reconstruction/
addpath(genpath('.'))

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
rmpath(genpath('/home/ybr19/Data'))

% List with all the names of the studies
study = 4;
if study ==1
    studies_29092020_DISORDER_Agarose;
    addpath(genpath('/home/ybr19/Data/29-09-2020_DISORDER_Agarose'))
    % Select studies to reconstruction
    id_path= 1;  id_file= [1];  id_ref = [1] ; 
elseif study ==2
    studies_YB_testrecon;
    addpath(genpath('/home/ybr19/Data/YB_testrecon'))
    % Select studies to reconstruction
    id_path = 5;
    id_file = [5]; % 4 is the alternating continuous motion (fig 5 in DISORDER7T.pdf)
    id_ref = [1] ; 
elseif study ==3
    studies_20201002_DISORDER_MotionFree;
    addpath(genpath('/home/ybr19/Data/2020-10-02_DISORDER_MotionFree'))
    % Select studies to reconstruction
    id_path= 1;  id_file= [9]; % 9 is the hybrid k-space
    id_file = {[1,2,3,7]};
    id_ref = 2*ones(size(id_file)); % Take sensitivities from 2nd low-res scan
    disp('Taking the sensitivity estimation of the second scan - assuming coil loading negligible.')
elseif study ==4
    studies_20201116_DISORDER_LowMotion;
    addpath(genpath('/home/ybr19/Data/2020-11-16_DISORDER_LowMotion'))
    id_path= 1;    
    id_file= [3];  
    %id_file = {[1:3]};
    id_ref = 1*ones(1,length(id_file)); % Take sensitivities from 1st low-res scan
end

[pathIn, refIn, fileIn] = extract_studies(pathIn, refIn, fileIn, id_path, id_ref, id_file);

%%% Set parameters
gpu=(gpuDeviceCount>0 && ~blockGPU);
ND=12; % YB:Dimensions you use to store all relevant data (first 3 are for spatial/spectral coordinates)
estSens_explicit = 1;
invData_explicit = 1;
runRecon_explicit = 1;

gibbsRinging=1;
corrMotion=1;
corrPhase=1;

supportReadS=[];supportRead=[]; % change depending on scan - especially if not acquired yourself

for p=1:length(pathIn)
    
    %BUILD THE REFERENCE DATA - COIL ESTIMATIONS
    for f=1:length(refIn{p})
        fileName=strcat(pathIn{p},filesep,refIn{p}{f});
        if ~exist(strcat(fileName,'recS.mat'), 'file') || estSens_explicit 
            TWS=mapVBVD(fileName);            
            recS=invert7T_old(TWS,supportReadS);
            recS=solveSensit7T(recS);
            save(strcat(fileName,'recS.mat'),'recS','-v7.3');    
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
        if ~iscell(fileName);fileNametemp = fileName;else fileNametemp = fileName{1};end
        if ~exist(strcat(fileNametemp,'rec.mat'),'file') || invData_explicit                      
            refName=strcat(pathIn{p},filesep,refIn{p}{f});           
            TW=mapVBVD(fileName);
            rec=invert7T_old(TW,supportRead);
            load(strcat(refName,'recS.mat'));
            NS=size(recS.S);NS=NS(1:3);
            NY=size(rec.y);NY=NY(1:3);
            NSY=max(round(NS.*recS.Enc.AcqVoxelSize./rec.Enc.AcqVoxelSize),NY(1:3));

            % Should be as below probably! since in DCT domain size==2*NS
            if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(2*NS,'tukeyIso',0.5,0,gibbsRinging,1);elseif gibbsRinging==-1;HY=buildFilter(2*NS,'CubicBSpline',[],0,1);else HY=[];end
            %if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(NS,'tukeyIso',0.5,0,gibbsRinging);elseif gibbsRinging==-1;HY=buildFilter(NS,'CubicBSpline',[],0);else HY=[];end
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
            rec.Assign.z{2}=TW.image.Lin-floor(NY(2)/2)-1; %YB: 2nd PE
            rec.Assign.z{3}=TW.image.Par-floor(NY(3)/2)-1; %YB: 3rd PE - slices
            
            if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))>sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med'))
                perm = [1 2 3 4]; % YB: second PE is the fastest - leave as it is
            else
                perm = [1 3 2 4]; % YB: put fastest PE in 2nd dimension
                rec.Assign.z{3}=TW.image.Lin-floor(NY(2)/2)-1;  
                rec.Assign.z{2}=TW.image.Par-floor(NY(3)/2)-1; 
                tt = rec.Par.Scan.MPS;
                tt(4:5) = rec.Par.Scan.MPS(7:8); tt(7:8) = rec.Par.Scan.MPS(4:5);
                rec.Par.Scan.MPS = tt;
                fprintf('Re-ordering dimenions results in: %s\n',rec.Par.Scan.MPS);
                fprintf('Fastest PE direction: %s\n',rec.Par.Scan.MPS(4:5));
            end
            
            rec.S=permute(rec.S,perm);rec.y=permute(rec.y,perm);rec.M=permute(rec.M,perm);
            rec.Enc.AcqVoxelSize=rec.Enc.AcqVoxelSize(perm(1:3));
            rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);
            rec.Par.Mine.Arot=rec.Par.Mine.Arot(:,perm);   
            rec.Par.Mine.Atra(:,4)=rec.Par.Mine.Atra(perm,4);  
            %rec.Par.Mine.MTT=rec.Par.Mine.MTT(:,perm); 
            %rec.Par.Mine.MTT=inv(rec.Par.Mine.MTT);
            rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;%YB: need to update this one as well (see invert7T.m)
            
            rec.Names.Name=fileIn{p}{f};
            rec.Names.pathOu=pathIn{p};

            rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
            rec.Par.Labels.ZReconLength=1;
            
            rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
            rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
            rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
            rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;
                                   
            %save(strcat(fileName,'rec.mat'),'rec','-v7.3');          
        else           
            if ~iscell(fileName)
                load(strcat(fileName,'rec.mat'));              
            else
                for n=1:length(fileName)
                    load(strcat(fileName{n},'rec.mat'));                    
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
        rec.Alg.parXT.perc=[0.9 0.9 0.9];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
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
            rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
            rec.Alg.parXT.NWend = 1; 
            rec.Alg.parXT.groupSweeps=2;% Factor to group sweeps
            rec.Alg.parXT.redFOV=0.4;
            rec.Alg.parXT.exploreMemory = 0;
            rec.Alg.parXT.traLimXT=[0.009 0.007];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
            rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
            rec.Par.Labels.FatShiftDir= 'H';% probably because of error in MPS
            rec.Alg.parXT.writeInter=0;
            rec.Alg.parXT.meanT=1;
            rec.Alg.UseSoftMasking= 0;

            rec.Dyn.GPUbS = [6 7];%[2 4]
            rec.Dyn.MaxMem=[6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

            close all
            
            %% MotCorr
            rec.Plan.Suff='_MotCorr';
            rec.Alg.parXT.resolMax = 5;
            rec.Alg.parXT.saveFinal=0;
            rec.Plan.SuffOu='';
            rec.Alg.parXB.dephaseCorrection=0;
            rec.Alg.parXB.useSH = 0;
            rec.Alg.parXB.useTaylor = 0;
            rec.Alg.parXB.useSusc = 0;
rec.Alg.parXT.disableGrouping = 1;
rec.Alg.nExtern = 70;
rec.Alg.parXT.convTransformJoint = 1;
            %recD=solveXTB_ext(rec);  

            %visSegment(recD.d,[],0)
            %visMotion(recD,[],[],0)

            %% Basis functions
            rec.Plan.Suff='_DephCorr';
            rec.Alg.parXT.saveFinal=0;

            rec.Alg.parXB.useSH = 1;
            rec.Alg.parXB.BasisType = 'SH';% 'SHOW' / 'BSplines' - not yet implemented ar
            rec.Alg.parXB.SHorder = 0; % disable by setting < 0
            rec.Plan.SuffOu= sprintf('_SH%d', rec.Alg.parXB.SHorder);
            rec.Alg.parXB.nGD_B = 1;
            rec.Alg.parXB.deTaylor = 0;
rec.Alg.parXT.disableGrouping = 1;
            %solveXTB_ext(rec);  
            %visSegment(recD.d,[],0)
            %visMotion(recD,[],[],0)

            %% With Taylor model
            close all
            rec.Alg.parXB.useTaylor = 1;
            rec.Alg.parXB.orderTaylor = 1;
            rec.Plan.Suff='_DephCorr';
            rec.Plan.SuffOu= sprintf('_SH%dTaylor%d', rec.Alg.parXB.SHorder,rec.Alg.parXB.orderTaylor);
            rec.Alg.parXT.saveFinal=0;
          
            rec.Alg.parXB.Reg.lambda_sparse = 0e+0; %1e+5; % Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_smooth = 0e+2; % 3e+5 Ad hoc - need to refine
            rec.Alg.parXB.Reg.lambda_susc = 0e+3;
            rec.Alg.parXB.nGD_Taylor = 2;
            rec.Alg.parXB.filterD.sp = 2;%In mm
            rec.Alg.parXB.filterD.gibbsRinging = 0.4;
            rec.Alg.parXB.deSH = 1;
            rec.Alg.parXB.Optim.FISTA=1;
            rec.Alg.parXB.Optim.RMSProp=0; rec.Alg.parXB.Optim.beta=0.9;
            rec.Alg.parXB.Optim.DupLim=inf;
            rec.Alg.parXB.Optim.taylorDelay=10;
            rec.Alg.parXB.Optim.alphaList=10.^(-1* [-3:0.7:6]);
            rec.Alg.parXB.unwrapTaylor=1;
            rec.Alg.nExtern = 30;
 
            rec.Alg.parXB.weightedEst = 0;
            rec.Alg.parXB.redFOV = 0.3;
            rec.Alg.parXB.intraVoxDeph = 0;
            rec.Alg.parXB.weightedEst = 1;
            rec.Alg.parXT.percRobustShot=[0.25 0.75];%Percentiles for robust computation of expected inter-shot dispersion
            rec.Alg.parXT.enerRobustShot=0.9;%0.95;%Ratio of the error for acceptance of shots
rec.Alg.parXT.corrFact = 1;
rec.Alg.parXT.disableGrouping = 1;
rec.Alg.parXB.useSusc = 0;
rec.Alg.parXT.convTransformJoint = 1;%all motion parameters estimated at every iteration

            rec.Alg.parXB.SHorder = 1; % disable by setting < 0
            rec.Plan.SuffOu= sprintf('_SH%dTaylor%d', rec.Alg.parXB.SHorder,rec.Alg.parXB.orderTaylor);
            recD=solveXTB_ext(rec);
                        
                %rec.Alg.parXB.useLin = 0;
                %rec.Alg.parXB.useSusc = 1;
                %rec.Alg.parXB = rmfield(rec.Alg.parXB,'Reg');
                %rec.Alg.parXB.Reg.lambda_smooth = 1e+8; % Ad hoc - need to refined
                                
                %solveXTB_ext(rec);
                %visSegment(recD.d,[],0)
                %visMotion(recD,[],[],0)
                %visSegment(dynInd(recD.D,1,6),[],0)
                %visSegment(dynInd(recD.D,2,6),[],0)
                %visBasisCoef(recD.E.Db, [], recD.timeState)
                
        end
    end

end



