
function rec = TWIX2rec(TW, suppFOV, isPT, shimDir, removeOversampling, noiseFileName, fileName, resRec, useGPU, bb)

%TWIX2REC Inverts SIEMENS TWIX data to a reconstruction structure.
%   * TW is the TWIX object returned from the mapVBVD.m function.
%   * {SUPPFOV} is the support in the readout direction (defaults to empty).
%   * {ISPT} is a flag to detect PT signal in the over-sampled k-space which needs to be extracted (defaults to 0).
%   * {SHIMDIR} is directory with additional shim calibration information of the scanner.
%   * {REMOVEOVERSAMPLING} is a flag whether to remove the Siemens default oversampling (defaults to 1).
%   * {NOISENAME} is the name of the file where noise samples are stored to be used (especially useful when noise samples are corrupted in an acquisition, as in Pilot Tone)
%   * {FILENAME} is the name of the file so that we can use it for logging and saving snapshots of parsing processes.
%   * {RESREC} is the resolution at which to save the reconstruction structure (defaults to empty = native resolution).
%   ** REC is a reconstruction structure.
%

%TO DO: manage noiseFile better. Add the option to convert file
%TO DO: include under-sampling

if nargin<2 || isempty(suppFOV);suppFOV=[];end
if nargin<3 || isempty(isPT);isPT=0;end
if nargin<4 || isempty(shimDir);shimDir='/home/ybr19/Projects/B0Shimming/';end
if nargin<5 || isempty(removeOversampling);removeOversampling=1;end
if nargin<6 || isempty(noiseFileName);noiseFileName='';end
if nargin<7 || isempty(fileName);fileName=[];end
if nargin<8 || isempty(resRec);resRec=[];end
if nargin<9 || isempty(useGPU);useGPU=1;end

%TEMPORARY 
if nargin<10 || isempty(bb);bb=1;end

%% SET LOGGING INFORMATION
if ~isempty(fileName)
    exportFlag=1;logFlag=1;
    [exportFolder,acqName]=fileparts(fileName);
    logFolder = fullfile(exportFolder,'Parsing_Log'); if ~ exist(logFolder,'dir');mkdir(logFolder);end
    exportFolder = fullfile(exportFolder,'Parsing_Sn'); if ~ exist(exportFolder,'dir');mkdir(exportFolder);end
else
    exportFlag=0;logFlag=0;
end

%TEMPORARY 
exportFlag=0
%logFlag=0
coilToPlot = [];

if TW.image.NRep>1 || TW.image.NSeg>1 || TW.image.NAve>1 || TW.image.NSli>1; exportFlag=0;end%Plotting won't work yet for NY(5)>1

if logFlag
    logName = strcat( logFolder, filesep, acqName, '.txt'); 
    if exist(logName, 'file'); delete(logName) ;end; if exist(logFolder,'dir')==0; mkdir(logFolder);end
    diary(logName) %eval( sprintf('diary %s ', logName))
    tStart = tic;
end

gpu= useGPU && (gpuDeviceCount>0 && ~blockGPU);
if gpu; gpuDevice(gpuDeviceCount);end%Take last one

%% DEFINE PARAMETERS FOR DATA READOUT
%Initialise reconstruction structure
rec=[];

%Parameters regarding TWIX object
ND=length(TW.image.dataSize); %Number of dimensions that can be called in the TWIX object (see public method subsref)
TW.image.flagIgnoreSeg=1;%See mapVBVD tutorial

%Number of receiver channels
NCha=TW.image.NCha;
%NCha=min(NCha,2);

%FOV array size 
NImageCols = TW.hdr.Config.NImageCols;%TW.hdr.Config.BaseResolution
NImageLins = TW.hdr.Config.NImageLins;
NImagePars = TW.hdr.Config.NImagePar;

%K-Space array size 
if isfield(TW.hdr.Config,'NColMeas'); NKSpaceCols = TW.hdr.Config.NColMeas; else; NKSpaceCols=TW.image.NCol;end
NKSpaceLins = TW.hdr.Config.NLinMeas;
NKSpacePars = TW.hdr.Config.NParMeas;

assert(TW.hdr.Dicom.flReadoutOSFactor*NImageCols==NKSpaceCols,'TWIX2rec:: Parameters involving over-sampling dealing not consistent.')

%Acceleration
if TW.hdr.MeasYaps.sPat.ucPATMode>1;rec.Enc.UnderSampling.isAccel = 1;else; rec.Enc.UnderSampling.isAccel = 0;end %ucPATMode: No accel (1) GRAPPA (2) or mSENSE (3)
if isfield(rec.Enc, 'UnderSampling') && rec.Enc.UnderSampling.isAccel 
    %Method of acceleration
    if  TW.hdr.MeasYaps.sPat.ucPATMode==2; rec.Enc.UnderSampling.accelMethod = 'GRAPPA'; else ; rec.Enc.UnderSampling.accelMethod = 'mSENSE'; end
    %Type of method
    if isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==2
        rec.Enc.UnderSampling.accelType = 'Integrated'; 
    elseif isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==4
        rec.Enc.UnderSampling.accelType = 'Non-Integrated'; 
    else
        rec.Enc.UnderSampling.accelType = []; 
    end
    %TODO: assign autocalibration lines to rec.ACS
        %if strcmp(rec.Enc.UnderSampling,'Non-Integrated'); rec.ACS = dynInd(TW.refscan,':',ND);end %TODO: this needs to undergo similar processing as rec.y
    %Undersampling factors and corresponding k-space indices to extract
    rec.Enc.UnderSampling.R = cat(2, TW.hdr.MeasYaps.sPat.lAccelFactPE , TW.hdr.MeasYaps.sPat.lAccelFact3D );

    idxAcqPE = cell(1,2); 
    idxAcqPE{1} = min(TW.image.Lin):rec.Enc.UnderSampling.R(1):NKSpaceLins;
    idxAcqPE{2} = min(TW.image.Par):rec.Enc.UnderSampling.R(2):NKSpacePars;    
end

%Trajectory
rec.Enc.ellipticalShutter = ~isempty(TW.hdr.Meas.ucEnableEllipticalScanning) && TW.hdr.Meas.ucEnableEllipticalScanning==1;

%% DATA READOUT
typ2Rec={'y'};

if isfield(TW,'noise') && ~isPT && ( isempty(noiseFileName) || exist(strcat( noiseFileName,'.mat'),'file')~=2 )
    rec.N=TW.noise.unsorted;%rec.N=dynInd(TW.noise,':',ND); 
    %ATTENTION: This is in k-space (I assumed generated without gradients playing). Since no gradients played, the PT signal will appear as a big DC component on top of the white noise.
    typ2Rec=[typ2Rec 'N'];%When noise loaded, precessing already done
elseif isPT
    if exist(strcat(noiseFileName, '.mat'),'file')==2%.mat file exists
        ss = load(noiseFileName);fprintf('Using noise samples from separate acquisition: %s\n',noiseFileName)
        if isfield(ss.rec,'N');rec.N = ss.rec.N;%Noise field exists
        else; warning('TWIX2rec:: No noise samples found in noiseName.mat provided. Noise de-correlation disabled!');%Noise field does not exist
        end
        rec.PT.noiseFileName = noiseFileName;
        clear ss;
    elseif exist(strcat(noiseFileName, '.dat'),'file')==2%.dat file exists
        recN = dat2Rec(noiseFileName,[],[],1);%Save .mat file
        rec.N = recN.N;recN=[];
        rec.PT.noiseFileName = noiseFileName;
        fprintf('TWIX2rec:: noiseName.dat converted and noise files used from converted structure.\n');
    else
        warning('TWIX2rec:: No file noiseName.mat or noiseName.dat found. Noise de-correlation disabled!');
    end
end

if isfield(rec,'N') && gpu;rec.N=gpuArray(rec.N);end

y=[];%Not yet in GPU (high res scans would be too big with over-sampling) 
perm=1:ND;
perm(2:4)=[3:4 2]; %Re-arrange so that 1st = Readout (RO), 2nd = First PE direction (PE1 or Lines), 3rd = Second PE direction (PE2 or Partitions) and 4th = Channel

reverseStr = '';
for s=1:NCha
    msg = sprintf('Processing channel: %d/%d\n',s,NCha);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    rec.y=dynInd(TW.image,{s ':'},[2 ND]); %RO-Channel-PE1-PE2
    if bb && rec.Enc.UnderSampling.isAccel; rec.y = dynInd(rec.y, idxAcqPE , 3:4);end
        
    if gpu;rec.y=gpuArray(rec.y);end 
    
    for n=1:length(typ2Rec)
        t2r=typ2Rec{n}; 
        rec.(t2r)=permute(rec.(t2r),perm);
        NY=size(rec.(t2r));NY(end+1:ND)=1;    
        
        if exportFlag && strcmp(t2r,'y') 
            APhiRecPlot = getPCS2RAS(TW.hdr.Dicom.tPatientPosition) * dynInd( diag(ones([1 4])) , {1:3,1:3},1:2,convertNIIGeom( ones([1,3]), TW.image.slicePos(4:7,1)', 'qForm', 'sForm'));
            [~,APhiRecPlot] = mapNIIGeom([], APhiRecPlot,'permute', [2 1 3]);  
            if isPT && ismember(s,coilToPlot)
                %Plot full k-space
                plotND({angle(dynInd(rec.y,1,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .7*dynInd(abs(rec.y),ceil((NY(1:3)+1)/2),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_3D_%d',acqName,s));
                %Plot central line of k-space
                line = dynInd(rec.y,ceil((NY(2:3)+1)/2),2:3).';
                visPTSignal (line, [], [],[], 1, [], [], 200, replace(sprintf('\\textbf{k-space:}\n %s - Coil %d',acqName,s),'_',' '), fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_Line_%d',acqName,s),'k-space readout [\#]',{'y','a.u.'});
            end
        end
    
        for l=1
            rec.(t2r)= fftshiftGPU(rec.(t2r),l); 
            rec.(t2r)= fftGPU(rec.(t2r),l); %YB: this has once given inf as number close to realmax('single') 
            rec.(t2r)= ifftshiftGPU(rec.(t2r),l);%if ifft here, you need to change PT on the bottom to ifft as well (before fft)
        end  
    end
    typ2Rec={'y'}; %All the channels are aready in the rec.N so loop only needs to permute + SupportFOV once
    y=cat(4,y,gather(rec.y));
    
    if exportFlag && ismember(s,coilToPlot)
        if isPT; plotND({angle(dynInd(rec.y,1,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .7*dynInd(abs(rec.y),ceil((NY(1:3)+1)/2),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{Hybrid k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilHybrid') , sprintf('%s_%d',acqName,s)); end
    end
end
rec.y=y;y=[];


%% EXTRACT PART OF READOUT TO AVOID FFT IN BOTH PE DIRECTIONS
%%% Remove PT signal
if isPT
    fprintf('Pilot Tone signal detection activated.\n');
    %Detect peak outside FOV
    NY = size(rec.y);
    rec.PT.yProjRO =  multDimSum( abs( dynInd(rec.y,ceil((NY(2:3)+1)/2),2:3 )),4); %Not in the readout
        zeroF=ceil((NY+1)/2);
        orig=zeroF(1)-ceil(((NImageCols)-1)/2); %YB: These are the elements to extract from the bigger F matrix to resample the array
        fina=zeroF(1)+floor(((NImageCols)-1)/2);
        idxTemp = orig(1):fina(1);
    yProjRO = dynInd(rec.PT.yProjRO, idxTemp, 1, 0);%Set FOV to if gpu;rec.y=gpuArray(rec.y);end zero. This assumes the PT is set outside the FOV
    [~,rec.PT.idxRO] = max(yProjRO,[],1) ;
    rec.PT.factorFOV = ( rec.PT.idxRO - ceil((NY(1)+1)/2) )/ ( NY(1)/2 );%If 0 -> Centre FOV / .5 -> Edge FOV / 1 -> Edge oversampled FOV
    
    %Extract PT signal (can be multiband signal)
    dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9; BWVox= 1 / dwellTime/ NImageCols;%In Hz
    rec.PT.FWHMHz = 0;
    try rec.PT.FWHMHz = getFWHM( BWVox * (1:length(yProjRO)) , yProjRO');catch;end        
    rec.PT.offsetMBHz=0;%Multiband in Hz to extract 
    if ~isfield(rec.PT,'offsetMBHz'); rec.PToffsetMBHz = 2*rec.PT.FWHMHz;end
    rec.PT.offsetMBidx=round(rec.PT.offsetMBHz / BWVox);
    rec.PT.idxMB = rec.PT.offsetMBidx +1;
    rec.PT.pSliceImage = dynInd(rec.y, (rec.PT.idxRO-rec.PT.offsetMBidx):(rec.PT.idxRO+rec.PT.offsetMBidx),1);%Extract from rec.y (image domain), not from yTemp
    
    %Average Multiband frequencies
    rec.PT.isAveragedMB = 0;
    if rec.PT.isAveragedMB
        rec.PT.pSliceImage = multDimMea(rec.PT.pSliceImage,1);%Note this is averaging in RO direction
        rec.PT.idxMB=1;
    end
    
    %Reference phase to avoid phase drift
    rec.PT.isRelativePhase = 1;
    if rec.PT.isRelativePhase
        rec.PT.relativePhaseCoilIdx = 1;
        rec.PT.pSliceImage =  ... %This in hybrid k-space, as it is supposed to be for phase subtraction
            rec.PT.pSliceImage .* conj( exp(1i*angle(dynInd(rec.PT.pSliceImage,rec.PT.relativePhaseCoilIdx,4)))); %Note this is averaging in RO direction
    end
    printPTInformation(rec);
end   

%%% Remove over-sampling
if removeOversampling
    fprintf('Removing oversampling.\n');
    rec.y=resampling(rec.y,NImageCols,2);%Remove over-encoding provided by Siemens (used e.g. in Pilot Tone) 
    if isfield(rec,'N');rec.N=resampling(rec.N,NImageCols,2);end
end
NImageColsArray=size(rec.y,1);%The FOV array kept in the readout direction

%%% Reduce FOV in readout dimension
if ~isempty(suppFOV)
    fprintf('Extracting part of FOV with range [%.2f --> %.2f].\n', suppFOV(1), suppFOV(2));
    NY = size(rec.y);
    vr=max(1,round(NY(1)*suppFOV(1))):min( NY(1),round(NY(1)*suppFOV(2))); %Only in the RO direction
    rec.y=dynInd(rec.y,vr,1); %Extract supported FOV
    if isfield(rec,'N');rec.N=dynInd(rec.N,vr,1);end %No need to remove oversampling for noise I think 
else
    vr=[];
end
 

%% GO TO IMAGE DOMAIN IN PE DIMENSIONS
%%% Check if gpu can handle all channels
yMem = getSize(rec.y,1);%In bytes
if gpu
    ss=gpuDevice; 
    safetyFact=.85;%Safety margin for GPU memmory
    if yMem > (safetyFact * ss.AvailableMemory); gpu=0; rec=gatherStruct(rec); end%Disable gpu and take everything out of it
end
if gpu;rec.y=gpuArray(rec.y);end

%%% Fourier transforms in PE dimensions - move all channels to image domain
NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    rec.y=fftshiftGPU(rec.y,n);
    rec.y=fftGPU(rec.y,n);
    rec.y=ifftshiftGPU(rec.y,n);
    
    if isPT
        rec.PT.pSliceImage=fftshiftGPU(rec.PT.pSliceImage,n);
        rec.PT.pSliceImage=fftGPU(rec.PT.pSliceImage,n);
        rec.PT.pSliceImage=ifftshiftGPU(rec.PT.pSliceImage,n);
    end
end

%%% Plot coil images to inspect later
if exportFlag
    if ~exist(fullfile( exportFolder ,'CoilImage'),'dir'); mkdir(fullfile(exportFolder ,'CoilImage'));end
    for s=1:coilToPlot
        plotND({abs(dynInd(rec.y,{s,1},[4 8])),1,0},angle(dynInd(rec.y,{s,1},[4 8])),[-pi pi 0 2.5*multDimMea(abs(dynInd(rec.y,{s,1},[4 8])))],[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{Image:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilImage') , sprintf('3D_%s_%d',acqName,s));;
    end
    plotND([],RSOS(dynInd(rec.y,1,8)),[0 .9*multDimMax(RSOS(dynInd(rec.y,1,8)))],[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{RSOS:} %s',acqName),'_',' '),[],[],fullfile(exportFolder ,'RSOS') , sprintf('%s',acqName));
end


%% PLOT PT SIGNAL
if isPT
    line = multDimSum(abs(rec.PT.yProjRO),[2:4]).';
    peakId = zerosL(line); peakId(rec.PT.idxRO)=1;
    visPTSignal (line, [], [],[], 1, peakId, [], 200, replace(sprintf('\\textbf{Projected hybrid k-space:}\n %s',acqName),'_',' '), [], [], 'Image readout [\# voxels]',{'y_{hybrid}','a.u.'});
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Peak',acqName)); end
    
    if gpu;rec.PT=gatherStruct(rec.PT);end
end

%% DE-CORRELATE NOISE
if isfield(rec,'N')
    fprintf('De-correlating channels from noise samples.\n');
    [rec.y, rec.Alg.covMatrix] = standardizeCoils(rec.y,rec.N);
    rec.Alg.decorrNoise=1;
else
    warning('TWIX2rec:: No noise-decorrelation applied!. Might result in sub-optimal performance.');
    rec.Alg.decorrNoise=0;
end
if exist('covMatrix','var')
    h=figure('color','w'); imshow(abs(rec.Alg.covMatrix),[]);title('Noise covariance matrix.');set(h,'color','w','Position',get(0,'ScreenSize'));
    if exportFlag; saveFig(fullfile(exportFolder ,'Noise','Covariance',acqName)); end
end

%% CONVERT TO PRS (PE-RO-SL), CONSISTENT WITH SIEMENS ORIENTATION CONVENTION
sliceInfo = TW.hdr.MeasYaps.sSliceArray.asSlice{1};

%%% DIMENSIONS IN IMAGE DOMAIN: RO-PE-SL
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.FOVmm = [sliceInfo.dReadoutFOV , sliceInfo.dPhaseFOV sliceInfo.dThickness];
rec.Enc.FOVSize = [NY(1), NImageLins, NImagePars];

rec.Enc.AcqVoxelSize=rec.Enc.FOVmm./ rec.Enc.FOVSize;
rec.Enc.AcqSize=[ 2^(TW.hdr.Dicom.flReadoutOSFactor==2)*NY(1) NY(2:3)];

%%% CHANGE TO PE-RO-SL                  
rec.y=gather(rec.y);if isfield(rec,'N');rec.N=gather(rec.N);end
perm=1:ND; perm(1:3)=[2 1 3];
rec.y = permute(rec.y,perm);
if isfield( rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
rec.Enc.AcqVoxelSize = rec.Enc.AcqVoxelSize(perm(1:3));%need to change spacing as well


%% GEOMETRY COMPUTATION - this will not work for under-sampled data
% asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1}; %See invert7T_old to see how to use this
% sNormal=asSlice.sNormal; sPosition=asSlice.sPosition;
slicePos = TW.image.slicePos(:,1);

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

%ROTATION 
quaternionRaw = slicePos(4:7);
rec.Par.Mine.Arot = eye(4);
rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS

%TRANSLATION
translationRaw = slicePos(1:3); %in mm for center FOV (I think)
rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

% account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = [ NImageLins, inf ,NImagePars];%Set inf to make sure this element is replaced
N(2)=NImageColsArray;%N(2) will depend if you removed over-sampling
orig = ( ceil((N(1:3)+1)/2) - [0 0 .5] )';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
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
rec.Par.Mine.patientPosition = TW.hdr.Dicom.tPatientPosition;
rec.Par.Mine.PCS2RAS = getPCS2RAS(rec.Par.Mine.patientPosition);
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
[~,rec.Par.Mine.APhiRec] = mapNIIGeom([], rec.Par.Mine.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSED ENGINEERED
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;%Backup as how was read in

%DEDUCE ACQUISITION ORDER
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane ] = acquisitionOrder(rec.Par.Mine.APhiRec);
rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
fprintf('\nSlice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
fprintf('Readout: %s\n', rec.Par.Scan.Mine.RO );
fprintf('1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
fprintf('2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

%% MAKE RO-PE-SL FOR DISORDER RECONSTRUCTION
perm=1:ND; perm(1:3) = [2 1 3 ];
rec.y = permute(rec.y, perm);
if isfield(rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'permute', perm);

%%% Store the permutations performed from the PRS space
rec.Par.Mine.permuteHist = [];rec.Par.Mine.permuteHist{1} = perm(1:4);

%%% Make MPS consisitent with RO-PE-SL
MPStemp = rec.Par.Scan.MPS ;
rec.Par.Scan.MPS(1:2) = MPStemp(4:5);
rec.Par.Scan.MPS(4:5) = MPStemp(1:2);

%% SEQUENCE INFORMATION
rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
rec.Par.Labels.ZReconLength=1;

rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;%In ms
rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;%In ms
rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};%In degrees
if isfield(TW.hdr.MeasYaps.sFastImaging,'ulEnableRFSpoiling'); rec.Par.Labels.RFSpoiling = TW.hdr.MeasYaps.sFastImaging.ulEnableRFSpoiling;else;rec.Par.Labels.RFSpoiling=0;end
rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;%In seconds
rec.Par.Labels.dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;%In seconds %*2 /1e+9 since alDwellTime{1} in nanoseconds and for oversampling
%https://www.magnetom.net/t/dwelltime-in-header-dat-and-rda-files-incorrect-for-spectroscopy/3342
%https://www.magnetom.net/t/how-to-get-the-spectroscopy-bandwidth-value-from-meas-as/320/4
rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz/FOV
rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / NImageCols;%In Hz/voxel in the readout direction
rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm

rec.Par.Labels.fieldStrength = TW.hdr.Dicom.flMagneticFieldStrength; %In T

%% SAMPLING INFORMATION
NY=size(rec.y);NY=NY(1:3);

%%% Take out under-sampling intervals
if rec.Enc.UnderSampling.isAccel
    Lines = ceil(TW.image.Lin/rec.Enc.UnderSampling.R(1));
    Partitions = ceil(TW.image.Par/rec.Enc.UnderSampling.R(2));
else
    Lines = TW.image.Lin;
    Partitions = TW.image.Par;
end

%%% Flip sampling and offset (probably because ambiguous Siemens fft vs. ifft convention)
if mod(NY(2),2)==0; offset =-1;else; offset=0; end
Lines = mod( Lines-1 + offset, NY(2) ) +1;
if mod(NY(3),2)==0; offset =-1;else; offset=0; end
Partitions = mod( Partitions -1 + offset, NY(3) ) +1;

%%% Shift to make compatible with DISORDER recon
kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE 
rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices
%ATTENTION: this mod() operator gives big spikes in the sampling. Question whether this a real thing or just processing....

rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};

%%%DISORDER information
if isfield(TW.hdr.MeasYaps,'sWipMemBlock') && isfield(TW.hdr.MeasYaps.sWipMemBlock,'alFree') && TW.hdr.MeasYaps.sWipMemBlock.alFree{1}==2
    rec.Enc.DISORDER=1;
else
    rec.Enc.DISORDER=0;
end
if rec.Enc.DISORDER
    rec.Enc.DISORDERInfo.tileSize=cat(2,TW.hdr.MeasYaps.sWipMemBlock.alFree{2},TW.hdr.MeasYaps.sWipMemBlock.alFree{3});
    rec.Enc.DISORDERInfo.segmentDuration = size(rec.Assign.z{2},2)/prod(rec.Enc.DISORDERInfo.tileSize)*rec.Par.Labels.RepetitionTime/1000;%In s - Use rec.Assign.z since NY does not take into account shutter
end

%% Make the PT signal in time to later check in solveXTB (optional)
if isPT
    tt = rec.PT.pSliceImage;
    
    %[kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    %filterStruct =[];filterStruct.medFiltKernelWidthidx=400/ rec.Par.Labels.RepetitionTime; 
    %tt = filterPTSlice (tt, filterStruct, NY, kIndex);

    for n=2:3
        tt=fftshiftGPU(tt,n);
        tt=fftGPU(tt,n)/NY(n);
        tt=fftshiftGPU(tt,n);%fftshift since in solveXT there is an iffthift on timeIndex
    end

    [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);

    tt = resSub(tt, 2:3);
    rec.PT.pTimeTest = squeeze(dynInd(tt, idx, 2)).';%sorted
    
    coilPlot=1:min(5,size(rec.PT.pTimeTest,1));
    %Plot original PT signal
    visPTSignal (dynInd(rec.PT.pTimeTest,coilPlot,1), [], [],[], 2, [], [], 200, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',acqName)) ; end
    %Pre-processed signal
    tt = bsxfun(@rdivide, rec.PT.pTimeTest , sqrt(normm(rec.PT.pTimeTest,[],1)) );
    tt = bsxfun(@minus, tt, multDimMea((tt),2) );
    visPTSignal (dynInd(tt,coilPlot,1), [], [],[], 2, [], [], 201, replace(sprintf('\\textbf{Pre-processed PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','TracePreProcessed',acqName)) ; end
    
    %Plot trajectory
    permTraj=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
    kTraj=zeros([length(rec.Assign.z{2}) 2],'single');
    for n=1:2;kTraj(:,n)=rec.Assign.z{permTraj(n)}(:);end
    kTrajSS=reshape(kTraj,[length(rec.Assign.z{2}) 1 2]);
    visTrajectory(kTrajSS,0);%h=figure('color','w'); imshow(timeMat , []);title('Time indices of k-space sampling'); set(h,'color','w','Position',get(0,'ScreenSize'));
    sgtitle(replace(sprintf('\\textbf{Trajectory:} %s\n',acqName),'_',' '),'Interpreter','latex','Color',[1 1 1]*0,'FontSize',20);
    if exportFlag; saveFig(fullfile(exportFolder ,'Trajectory',acqName)); end
end

%% SHIM SETTINGS
%%% Extract shim currents
rec.Par.Labels.Shim.shimCurrents_au = zeros(13,1);
if isfield(TW.hdr.Dicom, 'lFrequency'); rec.Par.Labels.Shim.shimCurrents_au(1) = TW.hdr.Dicom.lFrequency;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetX'); rec.Par.Labels.Shim.shimCurrents_au(2) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetX;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetY'); rec.Par.Labels.Shim.shimCurrents_au(3) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetY;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetZ'); rec.Par.Labels.Shim.shimCurrents_au(4) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetZ;end
if isfield(TW.hdr.Meas, 'alShimCurrent'); rec.Par.Labels.Shim.shimCurrents_au(5:13) = TW.hdr.Meas.alShimCurrent(1:9).';end
                                  
rec.Par.Labels.Shim.shimOrder = 3; %Given the 13 parameters in the 7T Terra console - for 3T this will be fine as the non-existent coils with have zero shim current                          

%%% Load correctionFactors and crossTerms
if ~isempty(shimDir)
    correctionFactorsFileName = strcat(shimDir,'Calibration', filesep, 'Results',filesep, 'correctionFactors.mat');
    if exist(correctionFactorsFileName,'file')==2
        ss = load(correctionFactorsFileName); rec.Par.Labels.Shim.correctFact_mA2au = ss.correctFact_mA2au;%Multiplicative correction
        rec.Par.Labels.Shim.shimCurrents_mA = rec.Par.Labels.Shim.shimCurrents_au ./ rec.Par.Labels.Shim.correctFact_mA2au;
    else
       fprintf('No Calibration data (correctionFactrs.mat) found in the shim directory.\n')  
    end

    crossTermFileName = strcat(shimDir,'Calibration', filesep, 'Results',filesep, 'crossTerms.mat');
    if exist(crossTermFileName,'file')==2
        ss = load(crossTermFileName); rec.Par.Labels.Shim.crossTerms = ss.crossTerms;
    else
       fprintf('No Calibration data (crossTerms.mat) found in the shim directory.\n')  
    end
else
   fprintf('No shim directory provided and hence no Calibration data read out.\n') 
end

%% Add additional info about PT
if isPT
    centreFOV = rec.Par.Mine.APhiRec * [ceil( (rec.Enc.FOVSize(1:3)-1 -.5 +1)/2) 1]';%-1 since NIFTI convention and not MATLAB 
    centreFOV = rec.Par.Mine.APhiRec * [(rec.Enc.FOVSize(1:3)/2 -[ 0 0 .5]) 1]';%-1 since NIFTI convention and not MATLAB
    RODir = rec.Par.Mine.APhiRec * [1;0;0;1] - rec.Par.Mine.APhiRec * [0;0;0;1];

    offsetFOV = multDimSum(RODir(1:3).*centreFOV(1:3))./sqrt(normm(RODir(1:3)));%Projection on RO direction

    rec.PT.transmitFreqHz = setPTFreq(...
       rec.Par.Labels.Shim.shimCurrents_au(1) ,...
       NImageCols,...
       rec.Par.Labels.voxBandwidth ,...
       rec.PT.factorFOV, ...
       offsetFOV,...
       rec.Enc.AcqVoxelSize(1))/1e6;
end

%% RESAMPLE STRUCTURE TO OTHER RESOLUTION/UNDER-SAMPLING
if ~isempty(resRec)
    if length(resRec)==1; resRec = resRec * ones([1 3]); else; resRec(end+1:3) = rec.Enc.AcqVoxelSize((length(resRec)+1):3);end
    fprintf('Resampling reconstruction to new resolutions %.2fmm (RO) / %.2fmm (PE1) %.2fmm (PE2).\n', resRec)
    rec = resampleRec(rec,resRec); %Resamples the arrays but also changes the sampling pattern accordingly
end

if logFlag; tStop = toc(tStart); fprintf('Data converted in %.0fmin %.0fs.\n',floor(tStop/60),mod(tStop,60)); diary off; end


end
