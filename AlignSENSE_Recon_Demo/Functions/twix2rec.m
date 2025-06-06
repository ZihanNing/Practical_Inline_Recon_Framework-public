function rec = twix2rec(twix,recon_flexFOV,removeZeroSamples,rawFilter)
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
    if nargin < 2 || isempty(recon_flexFOV); recon_flexFOV=[];end % ZN: reconstruct to recon matrix size defined in the sequence (usually cancelled the PE OS)
    if nargin < 3 || isempty(removeZeroSamples); removeZeroSamples=1;end 
    if nargin < 4 || isempty(rawFilter); rawFilter=0;end 
    gpu = 1 && (gpuDeviceCount>0 && ~blockGPU);
    
    
    %% MATRIX SIZES
    rec=[];
    %Parameters regarding TWIX object
    ND=length(twix.image.dataSize); %Number of dimensions that can be called in the TWIX object (see public method subsref)
    twix.image.flagIgnoreSeg=1;%See mapVBVD tutorial

    %Number of receiver channels
    NCha=twix.image.NCha;
    %NCha=min(NCha,2);

    %FOV array size 
    NImageCols = twix.hdr.Config.NImageCols;%TW.hdr.Config.BaseResolution
    NImageLins = twix.hdr.Config.NImageLins;
    NImagePars = twix.hdr.Config.NImagePar;

    %K-Space array size 
    if isfield(twix.hdr.Config,'NColMeas'); NKSpaceCols = twix.hdr.Config.NColMeas; else; NKSpaceCols=twix.image.NCol;end
    NKSpaceLins = twix.hdr.Config.NLinMeas;
    NKSpacePars = twix.hdr.Config.NParMeas;

    assert(twix.hdr.Dicom.flReadoutOSFactor*NImageCols==NKSpaceCols,'TWIX2rec:: Parameters involving over-sampling dealing not consistent.')
    
    %% Accleration 
    % ZN: here I used summarized acceleration-related parameters for
    % calculation, but could also use twix parameters for calculation (ref:
    % TWIX2rec)
    if twix.hdr.MeasYaps.sPat.ucPATMode>1;rec.Enc.UnderSampling.isAccel = 1;else; rec.Enc.UnderSampling.isAccel = 0;end %ucPATMode: No accel (1) GRAPPA (2) or mSENSE (3)
    if isfield(rec.Enc, 'UnderSampling') && rec.Enc.UnderSampling.isAccel 
        if  twix.hdr.MeasYaps.sPat.ucPATMode==2; accel_type = 'GRAPPA'; else ; accel_type = 'CAIPI'; end
        if isequal(accel_type,'CAIPI') && (removeZeroSamples~=0);removeZeroSamples = 0;end % ZN: for CAIPI acceleration, we need to use another different reconstruction method, thus no removal of unsampled points

        fprintf('Accelerated scan detected.\n');
        if isfield(twix,'refscan')
            rec.ACS = dynInd(twix.refscan,':',ND); % ACS data: [Col Cha Lin Par]
        else
            fprintf('No calibration scans detected, an extra REF scan needed!\n');
        end
        rec.Enc.UnderSampling.R = cat(2, twix.hdr.MeasYaps.sPat.lAccelFactPE , twix.hdr.MeasYaps.sPat.lAccelFact3D );
    end
    
    %% EXTRACT K-SPACE DATA
    %%% DATA & NOISE
    perm=1:ND;
    perm(2:4)=[3:4 2]; %Re-arrange so that 1st = Readout (RO), 2nd = First PE direction (PE1 or Lines), 3rd = Second PE direction (PE2 or Partitions) and 4th = Channel

    %%% Read data
    fprintf('Reading out channels: %d/%d\n',NCha,twix.image.NCha);
    rec.y=dynInd(twix.image,{1:NCha ':'},[2 ND]); %RO-Channel-PE1-PE2
       
    % ZN: add raw filter
    if rawFilter; rec.y = simple_raw_filter(rec.y,4,100,123,'exp',5,0); end% ZN: exp filter in Partition direction of the 1-2 row of tiles to get rid of the influence in the early echoes of the echo train (stripping artifact)
    % default: 112 123
    
    %%% ZERO-PADDING
    Lin = twix.image.Lin;
    Par = twix.image.Par;
    if any(ismember(fieldnames(twix.image),'IsRawDataCorrect'));idDataCorr = twix.image.IsRawDataCorrect;else;idDataCorr=[];end
    isShotBased = strcmp(twix.hdr.Meas.tScanningSequence,'GR\IR');% TW.hdr.MeasYaps.sFastImaging.lTurboFactor~=1;
    isMP2RAGE = isShotBased && twix.image.dataSize(10)>1;
    isMPRAGE = isShotBased && twix.image.dataSize(10)==1;
    isMeGRE = ~isShotBased && twix.image.dataSize(8)>1;
    isSPACE = contains(twix.hdr.Meas.tScanningSequence,'SE');
    if isShotBased && (~isMP2RAGE) % ZN; the TFE need to be recalculate for MP2RAGE
        %Shot ordering information
        if twix.hdr.MeasYaps.sKSpace.unReordering==1; TFE = twix.image.dataSize(4); %Linear
        else; TFE = twix.image.dataSize(3);%Linear-Rotated
        end
        numTI=twix.image.dataSize(10);
        %Make sure full inversions are extracted
        idFullInversion = 1: (floor(length(Lin)/(numTI*TFE))*(numTI*TFE));
        Lin = dynInd(Lin,idFullInversion,2); Par = dynInd(Par,idFullInversion,2); 
        if ~isempty(idDataCorr);idDataCorr = dynInd(idDataCorr,idFullInversion,2);end
        idFullInversion=[];
        %Remove duplicate sampling information
        if numTI>1
            id = makeBins(Lin, length(Lin)/TFE );
            id = mod(id,numTI)==1;
            Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); 
            if ~isempty(idDataCorr);idDataCorr=dynInd(idDataCorr,id,2);end
            id=[];
        end
    elseif isMP2RAGE % ZN: enable this for MP2RAGE
        if twix.hdr.MeasYaps.sKSpace.unReordering==1%Linear
            TFE = twix.image.dataSize(4); 
            TFE = findTFE(twix,3);
        else %Linear-Rotated
            TFE = twix.image.dataSize(3);
            TFE = findTFE(twix,2);
        end
        id = ceil(length(Lin)/TFE * ((1:length(Lin))-0.5)/length(Lin));
        id = mod(id,2)==1;
        Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); id=[];
    elseif isMeGRE
        %Remove duplicate sampling information
        numTE = twix.image.dataSize(8);
        id = 1:numTE:length(Lin);
        Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); 
        if ~isempty(idDataCorr);idDataCorr=dynInd(idDataCorr,id,2);end
        id=[];
    end
    
    %%% Pad for asymmetric echo
    NKSpaceColsGrid = twix.hdr.Meas.iRoFTLength;
    kSpacePad(1,1) = centerIdx(NKSpaceColsGrid) - twix.image.centerCol(1); % pre-pad
    kSpacePad(2,1) = NKSpaceColsGrid - twix.image.sqzSize(1) - kSpacePad(1,1); % post-pad
    if sum(kSpacePad); fprintf('Asymmetric echo detected: [Pre-pad, Post-pad] [%s, %s]\n',...
            num2str(kSpacePad(1,1)),num2str(kSpacePad(2,1))); end
    % Pad
    if multDimSum(kSpacePad (1,1))>0; rec.y = padArrayND(rec.y, [kSpacePad(1,1) 0 0 0],[],0,'pre'); end
    if multDimSum(kSpacePad (2,1))>0; rec.y = padArrayND(rec.y, [kSpacePad(2,1) 0 0 0],[],0,'post'); end

    
    %%% Pad for Lin and Par directions
    %Oversampling in both PE dimensions
    ROOverSamp = twix.hdr.Dicom.flReadoutOSFactor;
    PEOverSamp = [0 0];
    if ~isempty(twix.hdr.Meas.flPhaseOS); PEOverSamp(1) = twix.hdr.Meas.flPhaseOS;end
    if ~isempty(twix.hdr.Meas.flSliceOS); PEOverSamp(2) = twix.hdr.Meas.flSliceOS;end
    
    %Grid K-Space array size 
    NKSpaceColsGrid = twix.hdr.Meas.iRoFTLength;
    if ~isSPACE
        NKSpaceLinsGrid = twix.hdr.Config.NImageLins*(1+PEOverSamp(1))*twix.hdr.MeasYaps.sKSpace.dPhaseResolution;
        NKSpaceParsGrid = twix.hdr.Config.NImagePar*(1+PEOverSamp(2))*twix.hdr.MeasYaps.sKSpace.dSliceResolution;
    else
        NKSpaceLinsGrid = twix.hdr.Config.NImageLins*(1+PEOverSamp(1))*twix.hdr.MeasYaps.sKSpace.dPhaseResolution;
        NKSpaceParsGrid=twix.hdr.Meas.lImagesPerSlab*(1+PEOverSamp(2))*twix.hdr.MeasYaps.sKSpace.dSliceResolution %SPACE
    end

    if mod(NKSpaceLinsGrid,1)~=0; warning('KSpace array size for PE_1 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceLinsGrid,round(NKSpaceLinsGrid)); NKSpaceLinsGrid = round(NKSpaceLinsGrid);end
    if mod(NKSpaceParsGrid,1)~=0; warning('KSpace array size for PE_2 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceParsGrid,round(NKSpaceParsGrid)); NKSpaceParsGrid = round(NKSpaceParsGrid);end

    matrix_size = [size(rec.y,1),NKSpaceLinsGrid,NKSpaceParsGrid]; % ZN: for the case recon_resol < acq_resol, to be tested
    Acq_matrix_size = size(permute(rec.y,[1 3 4 2]),1:3);
    Lin_center = max(twix.image.centerLin);
    Par_center = max(twix.image.centerPar);
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
    

    %% NOISE DECORRELATION
    %%% NOISE
    if isfield(twix,'noise') % read-in noise data
        rec.N = twix.noise.unsorted; % [Col, Cha, Rep] > [Col Rep 1 Cha] 
        fprintf('De-correlating channels from noise samples.\n');
        rec.N = permute(rec.N,[1 3 4 2]); % [Col, Cha, Rep] > [Col Rep 1 Cha]
        NCha = size(rec.N,4);
        rec.N = dynInd(rec.N,1:NCha,4);
        [rec.y, rec.Par.preProcessing.deCorrNoise.covMatrix,rec.Par.preProcessing.deCorrNoise.ccm] = standardizeCoils(rec.y,rec.N);
        rec.preProcessing.deCorrNoise.Flag=1;
    end    
    
    % Remove over-sampling
    if ROOverSamp>1
        fprintf('Removing oversampling.\n');
        rec.y=resampling(rec.y,floor(size(rec.y,1)/2),2);%Remove over-encoding provided by Siemens (used e.g. in Pilot Tone) 
        if isfield(rec,'ACS');rec.ACS=resampling(rec.ACS,size(rec.ACS,1)/2,2);end
    end
    
    %% MATRIX AND RECON INFO
    sliceInfo = twix.hdr.MeasYaps.sSliceArray.asSlice{1};

    %%% DIMENSIONS IN IMAGE DOMAIN: RO-PE-SL
    NImageLins = NKSpaceLinsGrid;
    NImagePars = NKSpaceParsGrid;
    NImageColsArray=size(rec.y,1);%The FOV array kept in the readout direction
    NY=size(rec.y);NY(end+1:ND)=1;
    rec.Enc.FOVmm = [sliceInfo.dReadoutFOV , sliceInfo.dPhaseFOV , sliceInfo.dThickness];
    if ROOverSamp<=1; rec.Enc.FOVmm(1) = rec.Enc.FOVmm(1)*2; end % ZN: consider the not removing over-sampling cases
    rec.Enc.FOVmm = rec.Enc.FOVmm .* (1 + [0 PEOverSamp]); %Recale the FOV that is returned (not including the over-sampling)
    rec.Enc.FOVSize = [NY(1), NImageLins, NImagePars]; %Work on the FOV of the acquired (over-sampled) data

    rec.Enc.AcqVoxelSize=rec.Enc.FOVmm ./ rec.Enc.FOVSize;
    rec.Enc.AcqSize=[ 2^(ROOverSamp==2)*NY(1) NY(2:3)];
    if ROOverSamp<=1; rec.Enc.AcqSize(1) = NY(1); end % ZN: take care of the not over sampling removal case

    %% GEOMETRY COMPUTATION
    %%% CHANGE TO PE-RO-SL                  
    rec.y=gather(rec.y);if isfield(rec,'N');rec.N=gather(rec.N);end
    perm=1:ND; perm(1:3)=[2 1 3];
    rec.y = permute(rec.y,perm);
    if isfield(rec,'ACS'); rec.ACS = permute(rec.ACS,perm); end
    if isfield( rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
    rec.Enc.AcqVoxelSize = rec.Enc.AcqVoxelSize(perm(1:3));%Need to change spacing as well
    
    %%% GEOMETRY COMPUTATION FOR RECONSTRUCTION ARRAY
    fprintf('Computing geometry for reconstruction array.\n');
    slicePos = twix.image.slicePos(:,1);

    %SCALING
    rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

    %ROTATION 
    quaternionRaw = slicePos(4:7);
    rec.Par.Mine.Arot = eye(4);
    rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS

    %TRANSLATION
    translationRaw = slicePos(1:3); %in mm for center FOV (I think)
    % rec.Par.Mine.tablePos = TW.hdr.Config.GlobalTablePosTra; %in mm for table position
    % if ~isempty(rec.Par.Mine.tablePos);translationRaw(3) = translationRaw(3) + rec.Par.Mine.tablePos;end

    rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

    %Account for fact that tranRaw is referred to the center of FOV, not the first element in the array
    N = [ NImageLins, inf ,NImagePars];%Set inf to make sure this element is replaced
    N(2) = NImageColsArray;%N(2) will depend if you removed over-sampling
    orig = centerIdx(N(1:3));
    orig = (orig - [1 1 1])';%1 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel).
    % updated by ZN: 1 should be applied to all 3 directions - it used to only
    % applied to Partition direction somehow by substract [0 0 .5]
    %The [1 1 1] is also the difference you find between the gadgetron PACS exported data and Siemens PACS exported data
    vr = []
    if ~isempty(vr); orig(2)=orig(2)-(vr(1)-1); end%YB:orig(2) since this is readout/only vr(1)-1 elements removed from array/minus sign since should be addition and below already negative sign

    rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)...
                            - rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;

    %COMBINED MATRIX
    rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
    rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

    % % Alternative way (works): 
    % [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'translate',orig);
    % if ~isempty(vr);[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'dynInd',{':',vr,':'},NYOrig(1:3),NY(1:3));end

    %MOVE TO RAS
    rec.Par.Mine.patientPosition = twix.hdr.Dicom.tPatientPosition;
    rec.Par.Mine.PCS2RAS = getPCS2RAS(rec.Par.Mine.patientPosition);
    rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
    % [~,rec.Par.Mine.APhiRec] = mapNIIGeom([],rec.Par.Mine.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED % ZN: seems to have some issue
    rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;%Backup as how was read in

    %DEDUCE ACQUISITION ORDER
    [rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane] = acquisitionOrder(rec.Par.Mine.APhiRec);
    rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
    rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
    rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
    fprintf('Slice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
    fprintf('    Readout: %s\n', rec.Par.Scan.Mine.RO );
    fprintf('    1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
    fprintf('    2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

    rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
    rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

    %%% GEOMETRY COMPUTATION FOR AUTOCALIBRATION ARRAY
    if isfield(rec,'ACS')
        fprintf('Computing geometry for autocalibration array.\n');
        YSize = [ NImageLins, NImageColsArray ,NImagePars];%Set inf to make sure this element is replaced
        ACSSize = multDimSize(rec.ACS,1:3);
        [rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS] = mapNIIGeom([], rec.Par.Mine.APhiRec,'resampling',[],YSize,ACSSize);%TO CHECK:AD HOC - REVERSE ENGINEERED
    end

    %%% MAKE RO-PE-SL FOR ALIGNED SENSE PIPELINE CONVENTION
    perm=1:ND; perm(1:3) = [2 1 3 ];
    rec.y = permute(rec.y, perm);
    if isfield(rec,'ACS');rec.ACS = permute(rec.ACS, perm);end
    if isfield(rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
    [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'permute', perm);
    if isfield(rec,'ACS');[rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS] = mapNIIGeom(rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS,'permute', perm);end

    %%% Store the permutations performed from the PRS space
    rec.Par.Mine.permuteHist = [];rec.Par.Mine.permuteHist{1} = perm(1:4);

    %%% Make MPS consisitent with RO-PE-SL
    MPStemp = rec.Par.Scan.MPS ;
    rec.Par.Scan.MPS(1:2) = MPStemp(4:5);
    rec.Par.Scan.MPS(4:5) = MPStemp(1:2);
    

    %% SEQUENCE INFORMATION
    % General sequence parameters
    fprintf('Reading sequence parameters.\n');
    rec.Par.Labels.RepetitionTime = twix.hdr.MeasYaps.alTR{1}/1000;%In ms
    numTE=size(rec.y,8);
    rec.Par.Labels.TE = cat(2,twix.hdr.MeasYaps.alTE{1:numTE})/1000;%In ms
    rec.Par.Labels.FlipAngle = cat(2,twix.hdr.MeasYaps.adFlipAngleDegree{1:size(rec.y,10)});%In degrees
%     rec.Par.Labels.RFSpoiling = 1;
    rec.Par.Labels.dwellTime = 2 * twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;
    rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz/FOV
    rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / rec.Enc.FOVSize(1); %In Hz/voxel in the readout direction
    rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm
    rec.Par.Labels.fieldStrength = twix.hdr.Dicom.flMagneticFieldStrength; %In T
    rec.Par.Labels.ScanDuration = twix.hdr.MeasYaps.lScanTimeSec; % ZN: this is essential parameter to shotBased sequence
    
    % tse related
    if isSPACE
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
        rec.Par.Labels.inversionTime = cat(2,twix.hdr.MeasYaps.alTI{1:length(twix.hdr.MeasYaps.alTI)})/1000;%In ms
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
    if isMPRAGE
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
    rec.Enc.DISORDER=1;
    if ~isShotBased %DISORDER parameters not reliable when having a shot based sequence (To be fixed)
        rec.Enc.DISORDERInfo.tileSize = cat(2,twix.hdr.MeasYaps.sWipMemBlock.alFree{2},twix.hdr.MeasYaps.sWipMemBlock.alFree{3});
        rec.Enc.DISORDERInfo.NShots=prod(rec.Enc.DISORDERInfo.tileSize);
        rec.Enc.DISORDERInfo.segmentDuration = size(rec.Assign.z{2},2)/rec.Enc.DISORDERInfo.NShots*rec.Par.Labels.RepetitionTime/1000;%In s - Use rec.Assign.z since NY does not take into account shutter
        %if rec.Enc.DISORDER; assert(rec.Par.Labels.NShots==rec.Enc.DISORDERInfo.NShots,'DISORDER encoded shot-based sequence should have compatible shot/segment durations.');end
    end


    %%% Shutter
    rec.Enc.ellipticalShutter = ~isempty(twix.hdr.Meas.ucEnableEllipticalScanning) && twix.hdr.Meas.ucEnableEllipticalScanning==1;

    
    %% MAKE FLEXIABLE RECONSTRUCTION FOV
    % ZN: handle the scenario for PE OS or you simply would like to
    % reconstruct to a flexible size of FOV
    if 1 % ZN: for REF, we fft it to its acquired matrix size should be okay
        reconMatrixSize = rec.Enc.FOVSize; % This could be not identical with the encoded matrixSize, such as if PE OS is not 0
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
end