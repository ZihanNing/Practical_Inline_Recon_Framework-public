

function [rec] = solveXT_gadgetron(rec,T)
    
    if nargin < 2 || isempty(T); T=[];end  % ZN: T is the estimated motion trace, you can apply it to anther reconstruction

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
    rec.Alg.UseSoftMasking= 2; % 0-hardMask; 1-softMask; 2-NoMask; 3-ZNMask

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
    if isfield(rec.Par.Scan,'Mine');if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end;end
    rec.Alg.parXT.disableGrouping = 0; % ZN: 0 - enable group combined (check motion between shots and shots)
    rec.Alg.parXT.convTransformJoint = 1;
    rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
    rec.Alg.parXT.meanT=0;
    
%     rec.Alg.disabledisplay = 1; % ZN: 1 - not display the image during the recon
    rec.Alg.computeEnergy = 0; % 0 - do not compute EnBefore & EnAfter during recon (to save some time) % ZN
    
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
    rec.Alg.resPyr =  [ .5   1];
    rec.Alg.nExtern = [ 15   1];
    rec.Alg.parXT.estT = [ 1  0 ];
    
%     if ~rec.Enc.DISORDER %Need to specify the shot definition as it will not be detected automatically.
%         shotDuration = 2;%In seconds
%         numShots = [];
%         if ~isempty(shotDuration);rec.Alg.parXT.sampleToGroup = shotDuration/(1e-3*rec.Par.Labels.RepetitionTime);end
%         if ~isempty(numShots);rec.Alg.parXT.sampleToGroup = round(length(rec.Assign.z{2})/numShots);end
%     end

%     %Set mask
%     rTh=[1 1];
%     rec.M = getEllipsMask(rec.M,rTh)==1;
%     rec.Alg.UseSoftMasking= 0; % 0-hardMask; 1-softMask; 2-NoMask; 3-ZNMask

    rec = solveXTB_PT_zihan(rec,T);  
    
    