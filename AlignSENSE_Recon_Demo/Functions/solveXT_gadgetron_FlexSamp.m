

function [img_Di,img_Aq] = solveXT_gadgetron_FlexSamp(rec, flag_generateAq)
    % ZN: for CAIPI undersampling, the Aq image (image before moco) need to be generated
    % seperately, which takes extra time
    % flag_generateAq: 0-only generate Di image (image after moco);
    % 1-generate both Aq & Di images
    
    %%% CHECK SETTINGS
    if ~isfield(rec,'Alg')
        rec.Alg.CoilCompress_flag=1; rec.Alg.CoilCompressNum = 17; % in default compress coil to 17
        rec.Alg.Gibbsring_flag=0; % in default, do gibbs ringing reduction with MRTRIX3   
    else
        if isfield(rec.Alg,'CoilCompress_flag')
            if rec.Alg.CoilCompress_flag==1
                if ~isfield(rec.Alg,'CoilCompressNum');rec.Alg.CoilCompressNum = 17; end% ZN: default, compress coil to 17 channel
            end
        else
            rec.Alg.CoilCompress_flag=1;rec.Alg.CoilCompressNum = 17; % ZN: do coil compression by default
        end
        if isfield(rec.Alg,'Gibbsring_flag')
            if rec.Alg.Gibbsring_flag==1 % ZN: do gibbs ringing reduction with matlab
                if ~isfield(rec.Alg,'gibbsRinging');rec.Alg.gibbsRinging = [0.35 0.35]; end% ZN: default, gibbs ringing reduction factor
            end
        else
            rec.Alg.Gibbsring_flag=0; % ZN: do gibbs ringing reduction with mrtrix3 by default
        end
    end

    %%% GENERAL RECONSTRUCTION PARAMETERS
%     rec.Alg.WriteSnapshots=1;       
    rec.Alg.UseSoftMasking= 0; % 0-hardMask; 1-softMask; 2-NoMask; 3-ZNMask

    rec.Par.Mine.Modal = 7; %YB: Anatomical Volumetric Encoding

    rec.Alg.parXT.writeInter=0;

    rec.Dyn.GPUbS = [6 7];%[2 4]
    rec.Dyn.MaxMem = [6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

    
    %%% MOTION CORRECTION
%     rec.Plan.Suff='';rec.Plan.SuffOu='';
    
    %Disabling B0 estimation
    rec.Alg.parXB.useSH = 0;
    rec.Alg.parXB.useTaylor = 0;
    rec.Alg.parXB.useSusc = 0;

    %Motion estimation parameters
    rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
    %rec.Alg.parXT.resolMax = 8;
    rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
    if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
    rec.Alg.parXT.disableGrouping = 1; % ZN: 0 - enable group combined (check motion between shots and shots)
    rec.Alg.parXT.convTransformJoint = 1;
    rec.Alg.parXT.UseGiRingi=0; %ZN: added
    rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
%     rec.Alg.parXT.meanT=0;
    
%     rec.Alg.disabledisplay = 1; % ZN: 1 - not display the image during the recon
    rec.Alg.computeEnergy = 0; % 0 - do not compute EnBefore & EnAfter during recon (to save some time) % ZN
    
    %%% ADD SPATIAL INFROMATION ABOUT COILS
    rec.Par.preProcessing.coilGeom.RAS = coilCentroids(rec.S,rec.Par.Mine.APhiRec);
    [~,rec.Par.preProcessing.coilGeom.idxIS] = sort(rec.Par.preProcessing.coilGeom.RAS(3,:));
    ccm = eye(size(rec.S,4));
    ccm = ccm(flip(rec.Par.preProcessing.coilGeom.idxIS),:);
    rec.Par.preProcessing.coilGeom.ccmIS = ccm;ccm=[];
    
    %%% COMPRESS COILS BEFORE CALLING MOCO
    % ZN: readd to gadgetron code
    if rec.Alg.CoilCompress_flag
    percCoil = rec.Alg.CoilCompressNum; % when you use -23, it means use 23 coils % used to be 0.9; usually used 17
        if ~isempty(percCoil)
            rec.preProcessing.compressCoils.NChaOrig = size(rec.y,4);
            rec.preProcessing.compressCoils.percCoil = percCoil;
            [rec.S,rec.y] = compressCoils(rec.S,percCoil,rec.y);
            rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction
            fprintf('Coil compression: Number of coils compressed from %s to %s elements (%s%%).\n',...
                num2str(rec.preProcessing.compressCoils.NChaOrig,'%.1f'),sprintf('%.1f',size(rec.y,4)),...
                num2str(percCoil*100,'%.1f'));
        end
    else
        rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction
    end

    %Optimisation scheme
%     rec.Alg.resPyr =  [ .5   1];
%     rec.Alg.nExtern = [ 15   1];
%     rec.Alg.parXT.estT = [ 1  0 ];
    nIt = 7;
    rec.Alg.resPyr =  [ .5  .5  1];
    rec.Alg.nExtern = [ nIt   nIt   1];
    rec.Alg.parXT.estT = [ 1 1 0 ];
    rec.Alg.parXT.PT.subDiv = [1 0 0];
    
    % Run
    if mod(length(rec.Assign.z{2})/round(rec.Par.Labels.TFEfactor),1)~=0%Inconsistent data
        rec.Alg.parXT.sampleToGroup=round(length(rec.Assign.z{2})/(rec.Par.Labels.NShots*size(rec.y,5)));
        rec.Par.Labels.TFEfactor=length(rec.Assign.z{2})/size(rec.y,5);%Fake steady state sequence 
    end
    
    % disable mask
    rec.M = ones(size(rec.M));
    
    rec.Alg.disabledisplay = 1; % 1 - not display the image during the recon % ZN
%     rec.Enc.UnderSampling.R = [2 2]; % ZN: manually set the acceleration factor temporally 
    recTemp = rec;
    dimToLoop = 8;suff = '';
    isMultiAve = length(size(rec.y))>4; % sometimes the average dimension will be merged into the 5th dimension
    if isMultiAve
        if length(size(rec.y))==5
            numAve = size(rec.y,5);
            ave_dim = 5;
        elseif length(size(rec.y))==6
            numAve = size(rec.y,6);
            ave_dim = 6;
        else
            disp('Wrong dimension configuration for multi-average sequences!!')
        end
    end
    if isMultiAve
        recTemp.y = cat(5,dynInd(rec.y,1,ave_dim),dynInd(rec.y,2,ave_dim)); 
        for i=2:3
            recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
        end
    else
        recTemp = rec;
        idToRun = 1:size(rec.y,dimToLoop);
        doEveryEchoIndividually = 0;
    end
    numShots = size(recTemp.y,5)*rec.Enc.DISORDERInfo.NShots;
%     numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
    recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);          
    temp = solveXTB_PT_FlexSamp_Gadget(recTemp);
    temp = gatherStruct(temp);
    img_Di = temp.d; % ZN: get the images after motion correction

    % Uncorrected
    if flag_generateAq
        recTemp.Alg.AlignedRec=1;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
        % ZN: for the caipi undersampling pattern, the Aq images need to be
        % generated seperately
        % while generating Aq image, do not reject the outliers
        recTemp.Plan.Suff='_UnCorr';recTemp.Plan.SuffOu='_outl';
        temp = solveXTB_PT_FlexSamp_Gadget(recTemp, zerosL(temp.T));
        temp = gatherStruct(temp);
        img_Aq = temp.d; % ZN: get the images before motion correction
    else
        fprintf('For this CAIPI undersampled kspace, Aq images will not be generated.\n');
        img_Aq = [];
    end
    
    