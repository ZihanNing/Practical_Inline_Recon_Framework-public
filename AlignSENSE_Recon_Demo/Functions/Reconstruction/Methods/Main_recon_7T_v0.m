
clc; close all;
clear all
cd /home/ybr19/Reconstruction//
addpath(genpath('.'))

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))

% List with all the names of the studies
studies_20201002_DISORDER_MotionFree;
addpath(genpath('/home/ybr19/Data/2020-10-02_DISORDER_MotionFree'))

% Select studies to reconstruction
id_path=1;
id_file=[2 ];
id_ref =[2 ] ; % All the same

[pathIn, refIn, fileIn] = extract_studies(pathIn, refIn, fileIn, id_path, id_ref, id_file);

% Set parameters
gpu=(gpuDeviceCount>0 && ~blockGPU);
ND=12; % YB:Dimensions you use to store all relevant data (first 3 are for spatial/spectral coordinates)
estSens_explicit = 0;
invData_explicit = 0;
runRecon_explicit = 1;

gibbsRinging=1;
corrMotion=1;
corrPhase=1;

%port=[];%[0.4 0.75]; 
if id_path==7;supportReadS=[0.4 0.85];supportRead=[0.45 0.8];
elseif ismember(id_path,8:9);supportReadS=[0.275 0.725];supportRead=[0.325 0.675];%supportRead=[0.275 0.725];
else supportReadS=[];supportRead=[];
end


for p=1:length(pathIn)
    
    %BUILD THE REFERENCE DATA - COIL ESTIMATIONS
    for f=1:length(refIn{p})
        fileName=strcat(pathIn{p},filesep,refIn{p}{f});
        if ~exist(strcat(fileName,'recS.mat'), 'file') || estSens_explicit 
            TWS=mapVBVD(fileName);            
            recS=invert7T(TWS,supportReadS);
            recS=solveSensit7T(recS);
            save(strcat(fileName,'recS.mat'),'recS','-v7.3');    
        end
    end
    
    %RECONSTRUCT
    for f=1:length(id_file)
        rec=[];
        if iscell(id_file)
            for n=1:length(id_file{f});fileName{n}=strcat(pathIn{p},filesep,fileIn{p}{id_file{f}(n)});end
        else
            fileName=strcat(pathIn{p},filesep,fileIn{p}{f});
        end
        if ~exist(strcat(fileName,'rec.mat'),'file') || invData_explicit                      
            refName=strcat(pathIn{p},filesep,refIn{p}{f});           
            TW=mapVBVD(fileName);
            rec=invert7T(TW,supportRead);
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
            rec.S=permute(rec.S,[1 3 2 4]);rec.y=permute(rec.y,[1 3 2 4]);rec.M=permute(rec.M,[1 3 2 4]);% YB: Now: HF-AP-LR
            rec.Enc.AcqVoxelSize=rec.Enc.AcqVoxelSize([1 3 2]);
            rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);
            rec.Par.Mine.Arot=rec.Par.Mine.Arot(:,[1 3 2 4]);                           

            rec.Names.Name=fileIn{p}{f};
            rec.Names.pathOu=pathIn{p};

            rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
            rec.Par.Labels.ZReconLength=1;
            
            NY=size(rec.y);NY(end+1:4)=1;
            rec.Assign.z{2}=TW.image.Par-floor(NY(2)/2)-1; %YB: Refers to AP direction
            rec.Assign.z{3}=permute(TW.image.Lin-floor(NY(3)/2)-1,[1 3 2]); %YB: Refers to LR direction %YB: added floor. Not sure if right
            
            rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
            rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
            rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
            rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;
                                   
            save(strcat(fileName,'rec.mat'),'rec','-v7.3');          
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

        %rec.Alg.parXT.resolMax=1;%2;%4;%Coarser resolution on which to compute motion (in mm)
        %rec.Alg.parXT.discardHighRes=0.4;%Those volumetric scans whose resolution is below this value (in mm) are not reconstructed
        rec.Alg.parXT.groupSweeps=16;%Factor to group sweeps
        rec.Alg.AlignedRec=2;
        rec.Alg.parXT.NWend=1;%No intra-shot correction
        rec.Alg.parXT.accel=[1 0];%[1 0]%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
        rec.Alg.parXT.maximumDynamics=6;%Maximum number of dynamics allowed
        rec.Alg.WriteSnapshots=0;
        
        %rec.Alg.parXT.fractionOrder=0.5;
        %rec.Alg.parXT.redFOV=0;%Factor to reduce the FOV in the inferior direction for motion estimation
        rec.Alg.parXT.perc=[0.9 0.9 0.9];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
        rec.Alg.parXT.saveFinal=1;%Save final results for further inspection
        rec.Dyn.Log=1;
        
        Batch = 1; %YB: see solveXTB
        rec.Par.Mine.Modal = 7; %YB: Anatomical Volumetric Encoding
        fileName=generateNIIFileName(rec);
        %fileName=strcat(fileName,rec.Plan.Suff);

        if ~exist(sprintf('%s%03d.mat',fileName,Batch),'file') || runRecon_explicit  % The file would contain raw reconstruction data of last levels.
            % this should have names with suffou! Change
            if ~corrMotion 
                rec.Plan.SuffOu='NoMot';%Suffix to add to output data
                rec.Plan.Suff='NoMot';%Suffix to add to images
                rec.Dyn.Typ2Rec=12;
                %rec.x=bsxfun(@times,sum(bsxfun(@times,conj(rec.S),rec.y),4),1./(normm(rec.S,[],4)+1e-6));
                %writeData(rec);

            elseif ~corrPhase
                rec.Plan.SuffOu='MotCorr';
                rec.Plan.Suff='MotCorr';
                                
                rec.Alg.parXT.groupSweeps=16;% Factor to group sweeps
                rec.Alg.parXT.accel=[0 0];%[1 0]%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
                rec.Alg.parXT.redFOV=0;%Only for phantom
                rec = solveXT(rec);
                visSegment(rec.xRes)

            else
                rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
                if isfield(rec , 'Plan.Types'); rec.Plan=rmfield(rec.Plan,'Types');end
                if isfield(rec , 'Plan.TypeNames'); rec.Plan=rmfield(rec.Plan,'TypeNames');end          
                if isfield(rec , 'Dyn.Typ2Wri'); rec.Dyn.Typ2Wri(end+1)=0;end
                rec.Plan.SuffOu='DephCorr32Group1';
                rec.Plan.Suff='DephCorr32Group1';
                rec.Alg.parXB.dephaseCorrection=32;
                rec.Alg.parXB.DCTdims=[1 1 1];
                rec.Alg.parXT.computeCSRecon = 0; %YB: Note that if this is activated, the typ2Rec(18) will automatically be set to 1
                rec.Alg.parXT.NWend = 1; % YB: No intra-shot correction
                rec.Alg.parXT.exploreMemory = 0;
                
                rec.Plan.Suff='_MotCorr';
                rec.Alg.parXT.saveFinal=0;
                rec.Plan.SuffOu='';
                rec.Alg.parXB.dephaseCorrection=0;
                rec.Alg.parXB.useSH = 0;
                rec.Alg.parXB.useLin = 0;
                rec.Alg.parXB.useSusc = 0;
                recD=solveXTB(rec); 
                
            end

        end
    end
end



