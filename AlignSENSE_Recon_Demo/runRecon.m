
function [] = runRecon(idFile, idRef, idRefB, seqType)

%% 0. Initialization
% List with all the names of the studies
studyNames;
B0In = fillCell(fileIn, '');
B1In = fillCell(fileIn, '');
isPT = fillCell(fileIn, 0);
supportReadout = fillCell(fileIn,[]);
resRec = fillCell(fileIn,[]);
if idRefB == 0
    refBIn{1}{1} = '';
    refBIn = fillCell(fileIn, '');
end
noiseFile = fillCell(fileIn,[]);

idPath= 1; 

supportReadoutRecon = fillCell(cell(1,length(idFile)),[]);
resRecRecon = fillCell(cell(1,length(idFile)),[]);

% Extract studies
[pathIn, fileIn, refIn, refBIn, B0In, B1In,...
fileUnique, refUnique, refBUnique, B0Unique, B1Unique, ...
isPTUnique,noiseFileUnique, supportReadoutUnique, resRecUnique, RDesiredUnique] ...
    =  extractStudies (pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath,idFile, idRef,[],[], isPT, noiseFile, supportReadout, resRec, RDesired);

% DATA CONVERSION PARAMETERS
estSensExplicit=0;
invDataExplicit = 0;

%% 1. Generate the input of reconstruction (REC based on the twix-like structure
%%% CONVERT TO REC STRUCTURE
for p=1:length(pathIn)  
    %BUILD THE ACQUISITION DATA - start with this if an acquisition file is also used as a reference, you don't read it ouy twice (dat2Se checks if rec structure exist whereas dat2Rec reads data out by default)
    for f=1:length(fileUnique{p})
        fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
        if (~exist(strcat(fileName,'.mat'),'file') || invDataExplicit) && ~strcmp(fileName,'')
            [rec, TW] = dat2Rec(fileName,supportReadoutUnique{p}{f},[],1,[],isPTUnique{p}{f},[],noiseFileUnique{p}{f},resRecUnique{p}{f},[],RDesiredUnique{p}{f});
        end                
    end
end
TW = [];

%% 2. Estimate coil sensitivity map 
% ZN: estimation method
solve_espirit = 1;
isFailed = 0;

fprintf('=========== Estimating coil sensitivity map using ACS line.\n'); 
sens_name = ['REF_DATA_ACS_',fileName,'.mat']; % ZN: the ref file will contain seq name in the file name
nameRef = fullfile( rec.Names.pathOu, sens_name);

recS = rec; % ZN: borrow the rec structure
recS.y = rec.ACS;
recS.Enc.AcqVoxelSize = rec.Enc.UnderSampling.ACSVoxelSize;
recS.Par.Mine.APhiRec = rec.Par.Mine.APhiACS;
recS.Plan.Suff=''; recS.Plan.SuffOu='';
recS=solveSensit7T(recS); 
save(nameRef, 'recS','-v7.3');

%%
    %BUILD THE REFERENCE DATA - SENSITIVITY MAPS
    useACS_flag = 1; % ZN: 1- use ACS line for coil sensitivity map estimation; 0 - use external ref scan for coil sensitivity map estimation
    if useACS_flag == 0 % ZN: use external ref
        for f=1:length(refUnique{p})
            fileRef=strcat(pathIn{p},filesep,refUnique{p}{f});
            fileRefB=strcat(pathIn{p},filesep,refBUnique{p}{f});
            if ( ~exist(strcat(fileRef,'.mat'), 'file') || ~ismember('recS',who('-file', strcat(fileRef,'.mat')) ) ) || estSensExplicit && ~strcmp(fileRef,'')
                dat2Se(fileRef,fileRefB,[],1,1,[]);
            end
        end
    else % ZN: use ACSline
        fileAcq=strcat(pathIn{p},filesep,fileIn{p}{f});%Is allowed to be a cell array
        if ( ~exist(strcat(fileName,'_ACS.mat'), 'file') || ~ismember('recS',who('-file', strcat(fileName,'_ACS.mat')) ) ) || estSensExplicit && ~strcmp(fileName,'_ACS')
            ACS2rec(fileAcq,1,1);
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
        if useACS_flag == 0 % ZN: use external ref for coil sensitivity map estimation
            rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f},useACS_flag);
        else % ZN: use ACS line for coil sensitivity map estimation
            rec=prepareRec(fileAcq, [], fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f},useACS_flag);
        end

        
        %%% ADD SPATIAL INFROMATION ABOUT COILS
%         rec.Par.preProcessing.coilGeom.RAS = coilCentroids(rec.S,rec.Par.Mine.APhiRec);
%         [~,rec.Par.preProcessing.coilGeom.idxIS] = sort(rec.Par.preProcessing.coilGeom.RAS(3,:));
%         ccm = eye(size(rec.S,4));
%         ccm = ccm(flip(rec.Par.preProcessing.coilGeom.idxIS),:);
%         rec.Par.preProcessing.coilGeom.ccmIS = ccm;ccm=[];

        
        %%% COMPRESS COILS BEFORE CALLING MOCO
        percCoil = 17; % when you use -23, it means use 23 coils % used to be 0.9
        if ~isempty(percCoil)
            rec.preProcessing.compressCoils.NChaOrig = size(rec.y,4);
            rec.preProcessing.compressCoils.percCoil = percCoil;
            [rec.S,rec.y] = compressCoils(rec.S,percCoil,rec.y);
            rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction
            fprintf('Coil compression: Number of coils compressed from %d to %d elements (%d%%).\n',rec.preProcessing.compressCoils.NChaOrig,size(rec.y,4),percCoil*100);
        end

        %%% GENERAL RECONSTRUCTION PARAMETERS
        rec.Alg.WriteSnapshots = 1;       
        rec.Alg.UseSoftMasking = 2; %%% ZN modified, 0 for hard masking
        % 1 for soft masking, using SOS of the coil sensitivity map
        % 2 for don't mask
        % 3 for ZN defined mask
        % for accleration by ZN
        rec.Alg.disabledisplay = 1; % 1 - not display the image during the recon % ZN
        rec.Alg.computeEnergy = 0; % 0 - do not compute EnBefore & EnAfter during recon (to save some time) % ZN
        rec.Alg.setResAniManually = 0; % 1 - set ResAni manually, instead of using the result by pyramidplan
        rec.Alg.set_resAni = [0.25 0.25 0.25];
        rec.Alg.nIt = [300 1];  % for CG; defualt [300 1]
        
        
        

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
        rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
        if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
        rec.Alg.parXT.disableGrouping = 0; % 0 - enable group combined (check motion between shots and shots)
        rec.Alg.parXT.convTransformJoint = 1;
        rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
        rec.Alg.parXT.meanT=0;

        %Optimisation scheme
        rec.Alg.resPyr =  [ .5   1];
        rec.Alg.nExtern = [ 15   1]; % number of iteration % defualt [25 1]
        rec.Alg.parXT.estT = [ 1  0 ];
        if ~rec.Enc.DISORDER %Need to specify the shot definition as it will not be detected automatically.
            shotDuration = 1;%In seconds
            numShots = [];
            if ~isempty(shotDuration);rec.Alg.parXT.sampleToGroup = shotDuration/(1e-3*rec.Par.Labels.RepetitionTime);end
            if ~isempty(numShots);rec.Alg.parXT.sampleToGroup = round(length(rec.Assign.z{2})/numShots);end
        end
        
        % Set EllipsMask % add by ZN 
%         rec.Alg.UseSoftMasking = 0; % ZN modified, used to be 0, for hard masking
%         rTh=1;
%         rec.M = getEllipsMask(rec.M,rTh)==1;
        

        %Run
        if seqType == 0 
            solveXTB_PT_me(rec,[]);
        elseif seqType == 1 
            if  ~rec.Enc.DISORDER % ZN: to deal with non-DISORDER trajectory case
                rec.Alg.parXT.sampleToGroup = round( rec.Par.Labels.TFEfactor/prod([1 1])); % ZN: for non-DISORDER one, just one shot -> no motion correction at all
            else
                rec.Alg.parXT.sampleToGroup = round( rec.Par.Labels.TFEfactor/prod(rec.Enc.DISORDERInfo.tileSize));
            end
            
            recTemp = rec;
            for e=1:size(rec.y,8)
                recTemp.Plan.Suff='_MotCorr';
                recTemp.Plan.SuffOu=sprintf('_echo%d',e);
                recTemp.y = dynInd(rec.y,e,8);
                if e==1
                   temp =  solveXTB_PT_me(recTemp,[]);
                   T = temp.T;temp = [];
                else
%                    solveXTB_PT_me(recTemp,T);
                   temp =  solveXTB_PT_me(recTemp,[]);
                   T = temp.T;temp = [];
                end
            end
        elseif seqType == 2
            recTemp = rec;
            for e=1:size(rec.y,10)
                recTemp.Plan.Suff='_MotCorr';
                recTemp.Plan.SuffOu=sprintf('_INV%d',e);
                recTemp.y = dynInd(rec.y,e,10);
                solveXTB_PT_me(recTemp,[]);
            end
            
        end

    end
end

