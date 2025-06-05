function rec = twixlike2rec(twix,recon_flexFOV,removeZeroSamples,rawFilter)
%%% This function is to extract information from the twix-like data
%%% (converted from ISMRMRD raw) to generate a rec structure as the input
%%% of DISORDER reconstruction
%%%
%%% Some general functions:
%%%
%%% by Zihan Ning @ King's
%%% 17-Mar-2025

    fprintf('----------- Generate REC based on converted twix-like data (from ISMRMRD raw).\n');
    
    %% PRE-DEFINED PARAMETERS
    [seq_type,~] = recogSeqType(twix.hdr.sequenceParameters.sequence_type,...
        twix.hdr.Meas.tScanningSequence,1);
    
    if nargin < 2 || isempty(recon_flexFOV); recon_flexFOV=[];end % ZN: reconstruct to recon matrix size defined in the sequence (usually cancelled the PE OS)
    if nargin < 3 || isempty(removeZeroSamples); removeZeroSamples=1;end 
    if nargin < 4 || isempty(rawFilter); rawFilter=0;end 
    
    if isequal(seq_type,'mege'); readoutMode = 'bipolar'; else; readoutMode = [];end % ZN: should be possible to get this in twix header and converted to ismrmrd header
    TI = [880000,1000000]; % /1000 in ms; should be possible to get in twix header but now hard coded -> only used for MPRAGE I think
    %         rec.Par.Labels.inversionTime = cat(2,TW.hdr.MeasYaps.alTI{1:length(TW.hdr.MeasYaps.alTI)})/1000;%In ms
    
    % parameters from the header
    seq_name = twix.hdr.measurementInformation.protocolName;
    isRef = contains(seq_name,'REF')&&~contains(seq_name,'exterREF'); % ZN: whether it is an external REF for coil sensitivity map estimation
    
    
    
    %% MATRIX SIZES
    rec=[];
    
    
    % ZN: You can find the matrix related parameters in twix.hdr
    % Acquired matrix size & fov_mm (twix.hdr.acqSpace)
    %     without zero-padding, with OS
    % Encoded matrix size (twix.hdr.encodedSpace)
    %     with OS, including zero-padded points due to partial fourier or asymmetric echo
    % Recon matrix size & fov_mm (twix.hdr.encodedSpace)
    %     without OS, including zero-padded points due to PF and asymmetric echo, considering reconstruct the image to a higher resolution
    % 
    % ZN: Double check before you use!!
    %     If any issue, try to correct them in the ParameterMap instead of
    %     in BucketToBuffer_matlab function or here, since the matrix size
    %     will also influence the receiving of raw & retrieving of image
    
    
    
    %% Accleration 
    % ZN: here I used summarized acceleration-related parameters for
    % calculation, but could also use twix parameters for calculation (ref:
    % TWIX2rec)
    if twix.hdr.acceleration.isAccel; rec.Enc.UnderSampling.isAccel = 1; else; rec.Enc.UnderSampling.isAccel = 0; end
    if isfield(rec.Enc, 'UnderSampling') && rec.Enc.UnderSampling.isAccel 
        accel_type = twix.hdr.acceleration.type; 
        if isequal(accel_type,'CAIPI') && (removeZeroSamples~=0);removeZeroSamples = 0;end % ZN: for CAIPI acceleration, we need to use another different reconstruction method, thus no removal of unsampled points
        if ~isRef
            fprintf('Accelerated scan detected.\n');
            if isfield(twix,'reference')
                rec.ACS = twix.reference.data; % ACS data: [Col Cha Lin Par]
            else
                fprintf('No calibration scans detected, an extra REF scan needed!\n');
            end
            rec.Enc.UnderSampling.R = twix.hdr.acceleration.R; % acceleration factor
        else % REF scan
            fprintf('Fetal Error: This scan is recognized as a reference scan, but it seems to be accelerated! \n');
        end
    end
    
    %% EXTRACT K-SPACE DATA
    %%% DATA
    rec.y = twix.data; % [Col, Cha, Lin, Par, Ave, Con], with OS, without zero-padding due to partial Fourier
    
    % ZN: add raw filter
    if rawFilter; rec.y = simple_raw_filter(rec.y,4,100,123,'exp',5,0); end% ZN: exp filter in Partition direction of the 1-2 row of tiles to get rid of the influence in the early echoes of the echo train (stripping artifact)
    % default: 112 123
    
    %%% ZERO-PADDING
    % ZN: taking care both partial Fourier and (reconstruction resolution < acquisition resolution)
    % ZN: only pad the Lin & Par, supposing that asymmetric echo is handled by AsymmetricEcho gadget
    matrix_size = twix.hdr.encodedSpace.matrixSize; % ZN: for the case recon_resol < acq_resol, to be tested
    Acq_matrix_size = twix.hdr.acqSpace.matrixSize;
    Lin_center = max(twix.image.centerLin);
    Par_center = max(twix.image.centerPar);
    Lin = twix.image.sort_Lin(:,1);
    Par = twix.image.sort_Par(:,1);
    [rec.y,kspace_Pad,Lin,Par] = zeroPadKSpace(rec.y,matrix_size,Acq_matrix_size,Lin_center,Par_center,Lin,Par,1);
    
    %%% REMOVE UNDERSAMPLING
    % ZN: you should remove undersampling after zero-padding > mind the data padding and Lin & Par
    if rec.Enc.UnderSampling.isAccel && removeZeroSamples
        rec.y = RmvUndersampling(rec.y,Lin,Par,rec.Enc.UnderSampling.R); % use updated LIn & Par after zero-padding here
    end
    
    %%% TRANSFER TO IMAGE DOMAIN
    % transfer data to image domain
    if exist('accel_type','var')
        if isequal(accel_type,'CAIPI');fftOrder = 'revert';else;fftOrder = 'ori';end % ZN: not sure whether we need to use revert fft due to flair or caipi; double check when you use
    else
        fftOrder = 'ori';
    end
    rec.y = kspace2image_gpu(rec.y,fftOrder); % [Col Lin Par Cha Ave Con]
    % transfer calibration to image domain
    if isfield(rec,'ACS'); rec.ACS = kspace2image_gpu(rec.ACS,fftOrder);end % [Col Lin Par Cha]
    
    % handle MEGE with bipolar readout (flip the even idx echo image)
    if isequal(readoutMode,'bipolar') && twix.image.NEco > 1 % multi-echo sequence
        Eco_flip = 2:2:size(rec.y,6);
        for i = 1:length(Eco_flip)
            rec.y(:,:,:,:,:,Eco_flip(i)) = flipdim(rec.y(:,:,:,:,:,Eco_flip(i)),1); 
        end
    end
    

    %% NOISE DECORRELATION
    %%% NOISE
    if isfield(twix,'noise') % read-in noise data
        if isRef % for external REF, we tend to do the noise decorrelation after body coil estimation & coil sensitivity map estimation (keepy consistency)
            rec.N = permute(twix.noise,[1 3 4 2]); % [Col, Cha, Rep] > [Col Rep 1 Cha] 
        else
            fprintf('De-correlating channels from noise samples.\n');
            rec.N = permute(twix.noise,[1 3 4 2]); % [Col, Cha, Rep] > [Col Rep 1 Cha]
            NCha = size(rec.N,4);
            rec.N = dynInd(rec.N,1:NCha,4);
            [rec.y, rec.Par.preProcessing.deCorrNoise.covMatrix,rec.Par.preProcessing.deCorrNoise.ccm] = standardizeCoils(rec.y,rec.N);
            rec.preProcessing.deCorrNoise.Flag=1;
        end
    end    
    
    %% MATRIX AND RECON INFO
    %%% DIMENSIONS IN IMAGE DOMAIN: [Col Lin Par Cha Ave Con]
    ND = 16; NY=size(rec.y);NY(end+1:ND)=1;
    % to recon to this size
    rec.Enc.FOVSize = twix.hdr.encodedSpace.matrixSize; rec.Enc.FOVSize(1) = NY(1);  % rmv RO OS, zero padded (contains PE OS)
    rec.Enc.FOVmm = twix.hdr.encodedSpace.FOVmm; rec.Enc.FOVmm(1) = twix.hdr.reconSpace.FOVmm(1); % rmv RO OS, zero padded (contains PE OS)
    rec.Enc.AcqVoxelSize=rec.Enc.FOVmm ./ rec.Enc.FOVSize;
    % backup the acquisition resolution here
    ROOverSamp = twix.hdr.Dicom.flReadoutOSFactor;
    PEOverSamp = [0 0];
    if ~isempty(twix.hdr.Meas.flPhaseOS); PEOverSamp(1) = twix.hdr.Meas.flPhaseOS;end
    if ~isempty(twix.hdr.Meas.flSliceOS); PEOverSamp(2) = twix.hdr.Meas.flSliceOS;end
    rec.Enc.AcqSize=[ 2^(ROOverSamp==2)*NY(1) NY(2:3)];
    
    
    %% GEOMETRY COMPUTATION
    % arrange input
    AcqVoxelSize = rec.Enc.AcqVoxelSize;
    matrixSize = rec.Enc.FOVSize;
    if isfield(rec,'ACS'); ACSSize = multDimSize(rec.ACS,1:3); end
    image_hdr = twix.image;
    hdr = twix.hdr;
    TablePosTra = twix.hdr.Dicom.lGlobalTablePosTra;
    % geom computation
    if isfield(rec,'ACS')
        geom = computeGeom(AcqVoxelSize,matrixSize,hdr,image_hdr,TablePosTra,1,ACSSize);
    else
        geom = computeGeom(AcqVoxelSize,matrixSize,hdr,image_hdr,TablePosTra,0); 
    end
    rec.Par.Mine = geom;
    rec.Par.Scan.MPS = geom.MPS;
    rec.Par.Scan.Mine.slicePlane = geom.slicePlane;
    rec.Par.Scan.Mine.RO = geom.RO;
    rec.Par.Scan.Mine.PE1 = geom.PE1;
    rec.Par.Scan.Mine.PE2 = geom.PE2;
    rec.Par.Labels.FoldOverDir = geom.FoldOverDir;
    rec.Par.Labels.FatShiftDir = geom.FatShiftDir;
    if isfield(geom,'ACSVoxelSize');rec.Enc.UnderSampling.ACSVoxelSize = geom.ACSVoxelSize;end
    

    %% SEQUENCE INFORMATION
    % General sequence parameters
    fprintf('Reading sequence parameters.\n');
    rec.Par.Labels.RepetitionTime = twix.hdr.sequenceParameters.TR;
    rec.Par.Labels.TE = twix.hdr.sequenceParameters.TE;
    rec.Par.Labels.FlipAngle = twix.hdr.sequenceParameters.flipAngle_deg;
%     rec.Par.Labels.RFSpoiling = 1;
    rec.Par.Labels.dwellTime = 2 * twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;
    rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz/FOV
    rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / rec.Enc.FOVSize(1); %In Hz/voxel in the readout direction
    rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm
    rec.Par.Labels.fieldStrength = twix.hdr.acquisitionSystemInformation.systemFieldStrength_T; %In T
    rec.Par.Labels.ScanDuration = twix.hdr.MeasYaps.lScanTimeSec; % ZN: this is essential parameter to shotBased sequence
    
    % tse related
    if isequal(twix.hdr.sequenceParameters.sequence_type,'TurboSpinEcho')
        rec.Par.Labels.turboFactor = twix.hdr.MeasYaps.sFastImaging.lTurboFactor;
        rec.Par.Labels.echoTrainDuration = twix.hdr.MeasYaps.sFastImaging.lEchoTrainDuration;
        rec.Par.Labels.echoSpacing = rec.Par.Labels.echoTrainDuration/rec.Par.Labels.turboFactor;
        rec.Par.Labels.RepetitionTime_Long = twix.hdr.sequenceParameters.TR;
        rec.Par.Labels.RepetitionTime = rec.Par.Labels.echoSpacing;%In ms;
        rec.Par.Labels.RFEchoTrainLength = twix.hdr.Meas.RFEchoTrainLength;
    end
    
    isShotBased = strcmp(twix.hdr.Meas.tScanningSequence,'GR\IR');  
    if isShotBased
        if twix.hdr.MeasYaps.sKSpace.unReordering == 1 %Linear --> Each shot covers all partition samples = conventional MP-RAGE
            rec.Par.Labels.shotLayout='Linear';
            rec.Par.Labels.TFEfactor = size(rec.y,3);
        elseif twix.hdr.MeasYaps.sKSpace.unReordering == 256  %Linear rotated --> Each shot covers all line samples
            rec.Par.Labels.shotLayout='Linear-Rotated';
            rec.Par.Labels.TFEfactor = size(rec.y,2);%TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor  
        end
        rec.Par.Labels.NShots = length(Par)/rec.Par.Labels.TFEfactor; %rec.Par.Labels.TFEfactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor
        fprintf('Shot-based sequence with %d shots (%d samples per shot).\n',rec.Par.Labels.NShots, rec.Par.Labels.TFEfactor)
        assert(mod(rec.Par.Labels.NShots,1)==0,'TWIX2Rec:: Number of samples per shot is not an integer. Probably a corrupted scan.');
        rec.Par.Labels.RepetitionTime_Long=rec.Par.Labels.RepetitionTime;
        rec.Par.Labels.RepetitionTime=rec.Par.Labels.TE*2;%Old TR was the TR per shot --> assumed TR is twice the echo time
        rec.Par.Labels.inversionTime = cat(2,TI(1),TI(2))/1000;%In ms
    else
        rec.Par.Labels.TFEfactor=length(Par); 
    end
    rec.Par.Labels.ZReconLength=1;
    
    
    %% PE TRAJECTORY
    fprintf('Reading sampling information.\n');
    NY=size(rec.y);NY=NY(1:3);

    %%% Take out under-sampling intervals
    if rec.Enc.UnderSampling.isAccel && removeZeroSamples
        Lines = ceil(Lin/rec.Enc.UnderSampling.R(1));
        Partitions = ceil(Par/rec.Enc.UnderSampling.R(2));
    else
        Lines = Lin;
        Partitions = Par;
    end

    %%% Flip sampling and offset (because Siemens fft vs. ifft convention)
    % ZN: it is a little bit tricky here, not sure whether mprage seems to
    % be an exception
    if isequal(seq_type,'mprage')
        offset =-1;
        Lines = mod( Lines-1 + offset, NY(2) ) +1;
        Partitions = mod( Partitions -1 + offset, NY(3) ) +1;
    else
        if mod(NY(2),2)==0;offset =-1;else;offset=0;end
        Lines = mod( Lines-1 + offset, NY(2) ) +1;
        if mod(NY(3),2)==0;offset =-1;else;offset=0;end
        Partitions = mod( Partitions -1 + offset, NY(3) ) +1;
    end

    %%% Shift to make compatible with DISORDER recon
    kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
    rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE 
    rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices
    %ATTENTION: this mod() operator gives big spikes in the sampling. This is just a processing implication.

    rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};


    %%%DISORDER related information
    if isRef
        rec.Enc.DISORDER=0;
    else
        rec.Enc.DISORDER=1;
        if ~isShotBased %DISORDER parameters not reliable when having a shot based sequence (To be fixed)
            rec.Enc.DISORDERInfo.tileSize = cat(2,twix.hdr.MeasYaps.sWipMemBlock.alFree{2},twix.hdr.MeasYaps.sWipMemBlock.alFree{3});
            rec.Enc.DISORDERInfo.NShots=prod(rec.Enc.DISORDERInfo.tileSize);
            rec.Enc.DISORDERInfo.segmentDuration = size(rec.Assign.z{2},2)/rec.Enc.DISORDERInfo.NShots*rec.Par.Labels.RepetitionTime/1000;%In s - Use rec.Assign.z since NY does not take into account shutter
            %if rec.Enc.DISORDER; assert(rec.Par.Labels.NShots==rec.Enc.DISORDERInfo.NShots,'DISORDER encoded shot-based sequence should have compatible shot/segment durations.');end
        end
    end

    %%% Shutter
    rec.Enc.ellipticalShutter = ~isempty(twix.hdr.Meas.ucEnableEllipticalScanning) && twix.hdr.Meas.ucEnableEllipticalScanning==1;

    
    %% MAKE FLEXIABLE RECONSTRUCTION FOV
    % ZN: handle the scenario for PE OS or you simply would like to
    % reconstruct to a flexible size of FOV
    if ~isRef % ZN: for REF, we fft it to its acquired matrix size should be okay
        reconMatrixSize = twix.hdr.reconSpace.matrixSize; % This could be not identical with the encoded matrixSize, such as if PE OS is not 0
        % ZN: now decide the actual recon FOV
        reconFOVSize = rec.Enc.FOVSize; 
        if ~isempty(recon_flexFOV) && ~isequal(recon_flexFOV,reconMatrixSize)
            reconFOVSize = floor(recon_flexFOV./rec.Enc.AcqVoxelSize); % ZN: the point is to not change the resolution, but to enlarge FOV by changing the matrixSize
        else
            reconFOVSize = reconMatrixSize;
        end
        % Calc RTarget based on PE OS
        rec.Enc.UnderSampling.RTarget = rec.Enc.UnderSampling.R./(rec.Enc.FOVSize(2:3)./reconFOVSize(2:3)); % PE OS
        % Change the geom & resol accordingly
        if ~isequal(reconFOVSize,rec.Enc.FOVSize) % ZN: change the reconstruction matrix size/FOV/VoxelSize/geom
            fprintf('Change the reconstruction matrix size from [%s %s %s] to [%s %s %s]; FOVmm, VoxelSize and geom changed accordingly!\n',...
                num2str(rec.Enc.FOVSize(1)),num2str(rec.Enc.FOVSize(2)),num2str(rec.Enc.FOVSize(3)),...
                num2str(reconFOVSize(1)),num2str(reconFOVSize(2)),num2str(reconFOVSize(3)));
            %Store original parameters if ever needed
            rec.Enc.UnderSampling.Acq.FOVSize = rec.Enc.FOVSize;
            rec.Enc.UnderSampling.Acq.AcqSize = rec.Enc.AcqSize;
            rec.Enc.UnderSampling.Acq.R = rec.Enc.UnderSampling.R;
            rec.Enc.UnderSampling.Acq.FOVmm = rec.Enc.FOVmm;
            rec.Enc.UnderSampling.Acq.PEOverSamp = PEOverSamp;
            %Assign new values
            if removeZeroSamples
                rec.Enc.FOVSize(2:3) = round(rec.Enc.UnderSampling.RTarget .* rec.Enc.AcqSize(2:3));
            else
                rec.Enc.FOVSize(2:3) = reconFOVSize(2:3);
            end
            if ~isequal(reconFOVSize,rec.Enc.FOVSize);warning('Might be error acceleration rate compuation!\n');end
            rec.Enc.FOVmm = rec.Enc.FOVSize .* rec.Enc.AcqVoxelSize;
            if removeZeroSamples
                rec.Enc.UnderSampling.R = rec.Enc.FOVSize(2:3)./rec.Enc.AcqSize(2:3);
            else
                rec.Enc.UnderSampling.R = rec.Enc.UnderSampling.RTarget;
            end
            % recompute geom
            [rec.Enc.AcqVoxelSize,rec.Par.Mine.APhiRec,rec.Par.Mine.APhiRecOrig]=alignResolGeom(rec.Enc.FOVSize,...
                rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.AcqVoxelSize,...
                rec.Par.Mine.APhiRec,rec.Par.Mine.APhiRecOrig,rec.Par.Mine.permuteHist);
            fprintf('Changing undersampling from R=%.1fx%.1f to R=%.2f-%.2f by changing the FOV size.\n',...
                rec.Enc.UnderSampling.Acq.R(1),rec.Enc.UnderSampling.Acq.R(2),...
                rec.Enc.UnderSampling.R(1),rec.Enc.UnderSampling.R(2));
        end
    end
    
    %% REPORT
    if isfield(twix.hdr,'fileName')
        rec.Names.Name = twix.hdr.fileName; % ZN: with time tag to prevent overwriting
    else
        rec.Names.Name = twix.hdr.measurementInformation.protocolName;
    end
    
    % added by ZN % identify the folder name for saving (patientID and
    % study ID)
    if isfield(twix.hdr,'subjectInformation') && isfield(twix.hdr.subjectInformation,'patientID')
        rec.Names.pathOu = dumpDir_ID(twix.hdr.subjectInformation.patientID);
    else % cannot find the information 
        rec.Names.pathOu = dumpDir();  % just save under DATA_DUMP folder
    end

    %%% WRITE METADATA TO JSON
    recJSON = rec;
    %Remove all memmory intensice arrays
    recJSON.y=[];
    recJSON.S=[];recJSON.N=[];recJSON.Assign=[];
    recJSON.x=[];recJSON.M=[];
    recJSON.Par.preProcessing = []; % ZN: specific for noise decorrelation
    recJSON.ACS=[];
    %Save JSON
    nameSave = sprintf('%s/%s.json',rec.Names.pathOu, rec.Names.Name);
    savejson('',recJSON,nameSave);%writeJSON has specific fields 
    
    %Plot JSON in the command view to see the parameters on the host
    fid = fopen(nameSave);
    raw = fread(fid,inf);
    str = char(raw')
    
    %%% SAVE NIFTI: added by Ning
    isRef = contains(rec.Names.Name, 'REF');
    writeNIIFlag = 1;
    if writeNIIFlag>0
        outDir = fullfile(rec.Names.pathOu,'Parsing'); if exist(outDir,'dir')~=2; mkdir(outDir);end
        xW=[];xW{1} = RSOS(rec.y); 
        MSW=[];MSW{1} = rec.Enc.AcqVoxelSize; 
        MTW =[];MTW{1} = rec.Par.Mine.APhiRec;
        if isRef
            writeNII( strcat(outDir,filesep,rec.Names.Name) , {'y_Ref'},xW, MSW, MTW);
        else
            writeNII( strcat(outDir,filesep,rec.Names.Name) , {'y_highres'},xW, MSW, MTW);
        end
        fprintf('NIFTI file with raw coil images saved:\n   %s\n', strcat(outDir,filesep,rec.Names.Name));
    end
    if writeNIIFlag>0 && isfield(rec,'ACS') && ~isRef
        outDir = fullfile(rec.Names.pathOu,'Parsing'); if exist(outDir,'dir')~=2; mkdir(outDir);end
        xW=[];xW{1} = RSOS(rec.ACS); 
        MSW=[];MSW{1} = rec.Enc.UnderSampling.ACSVoxelSize; 
        MTW =[];MTW{1} = rec.Par.Mine.APhiACS;
        writeNII( strcat(outDir,filesep,rec.Names.Name) , {'ACS'},xW, MSW, MTW);
        fprintf('NIFTI file with ACS images saved:\n   %s\n', strcat(outDir,filesep,rec.Names.Name));
    end
    
end