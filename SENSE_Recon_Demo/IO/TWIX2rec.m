
function rec = TWIX2rec(TW, suppFOV, isPT, shimDir, removeOversampling, noiseFileName, fileName, resRec, useGPU, removeZeroSamples, targetR)

%TWIX2REC Inverts SIEMENS TWIX data to a (Philips) reconstruction structure.
%   * TW is the TWIX object returned from the dat2TWIX.m function.
%   * {SUPPFOV} is the support in the readout direction (defaults to empty).
%   * {ISPT} is a flag to detect PT signal in the over-sampled k-space which needs to be extracted (defaults to 0).
%   * {SHIMDIR} is directory with additional shim calibration information of the scanner.
%   * {REMOVEOVERSAMPLING} is a flag whether to remove the Siemens default oversampling (defaults to 1).
%   * {NOISENAME} is the name of the file where noise samples are stored to be used (especially useful when noise samples are corrupted in an acquisition, as in Pilot Tone)
%   * {FILENAME} is the name of the file so that we can use it for logging and saving snapshots of parsing process.
%   * {RESREC} is the resolution at which to save the reconstruction structure (defaults to empty = native resolution).
%   * {USEGPU} is a flag to enable gpu functionality (defaults to 1).
%   * {REMOVEZEROSAMPLES} is a flag to remove the un-acquired samples in regular undersampling (defaults to 1).
%   * {TARGETR} is the undersampling to change the acquisition to by changing the FOV size.
%   ** REC is the converted reconstruction structure.
%
%   Yannick Brackenier 2022-06-12

if nargin<2 || isempty(suppFOV);suppFOV=[];end
if nargin<3 || isempty(isPT);isPT=0;end
if nargin<4 || isempty(shimDir);shimDir='/home/ybr19/Projects/B0Shimming/';end
if nargin<5 || isempty(removeOversampling);removeOversampling=1;end
if nargin<6 || isempty(noiseFileName);noiseFileName='';end
if nargin<7 || isempty(fileName);fileName=[];end
if nargin<8 || isempty(resRec);resRec=[];end
if nargin<9 || isempty(useGPU);useGPU=1;end
if nargin<10 || isempty(removeZeroSamples);removeZeroSamples=1;end
if nargin<11 || isempty(targetR);targetR=[];end

assert(~isempty(TW),'TWIX2rec:: TWIX object cannot be empty.');

%% SET EXPORTING INFORMATION
if ~isempty(fileName)
    exportFlag=1;
    [exportFolder,acqName]=fileparts(fileName);
    exportFolder = fullfile(exportFolder,'Parsing_Sn'); if ~ exist(exportFolder,'dir');mkdir(exportFolder);end
else
    exportFlag=0;acqName='';
end

%%% ======= TEMPORARY 
exportFlag = 0
coilToPlot = []
%%% ======= TEMPORARY 

if TW.image.NRep>1 || TW.image.NSeg>1 || TW.image.NAve>1 || TW.image.NSli>1 || TW.image.NEco>1; exportFlag=0;end%Plotting won't work yet for NY(5)>1
Date = TW.hdr.Phoenix.tReferenceImage0; Date = strsplit(Date,'.');Date = strcat( Date{end}(7:8),'/',Date{end}(5:6),'/',Date{end}(1:4));

gpu = useGPU && (gpuDeviceCount>0 && ~blockGPU);
gpu = 0
if gpu; gpuDevice(gpuDeviceCount);end%Take last one

%% DEFINE PARAMETERS FOR DATA READOUT
%Initialise reconstruction structure
rec=[];
rec.Names.Date=Date;

%Parameters regarding TWIX object
ND=length(TW.image.dataSize); %Number of dimensions that can be called in the TWIX object (see public method subsref)
TW.image.flagIgnoreSeg=1;%See mapVBVD tutorial

%Number of receiver channels
NCha=TW.image.NCha;
%NCha=min(NCha,1);warning('Number of channels fixed to lower count.');

%Samples in both phase encoding dimensions
Lin = TW.image.Lin;
Par = TW.image.Par;

%Oversampling in both PE dimensions
ROOverSamp = TW.hdr.Dicom.flReadoutOSFactor;
PEOverSamp = [0 0];
if ~isempty(TW.hdr.Meas.flPhaseOS); PEOverSamp(1)= TW.hdr.Meas.flPhaseOS;end
if ~isempty(TW.hdr.Meas.flSliceOS); PEOverSamp(2)= TW.hdr.Meas.flSliceOS;end

%FOV array size 
%if isfield(TW.hdr.Meas,'NImageCols');NImageCols=TW.hdr.Meas.NImageCols;else NImageCols=TW.hdr.Config.BaseResolution;end%nCol is an alternative, also in Config
NImageCols = TW.hdr.Config.NImageCols;%TW.hdr.Config.BaseResolution#THIS IS PROBLEMATIC AS SOMETIMES IT IS EMPTY
% NImageLins = TW.hdr.Config.NImageLins;%Attention, this is "recon FOV", not what you acquired!
% NImagePars = TW.hdr.Config.NImagePar;%TW.hdr.Meas.lImagesPerSlab%Attention, this is "recon FOV", not what you acquired!

%Measured K-Space array size 
NKSpaceCols = TW.hdr.Config.NColMeas;%TW.hdr.Config.RawCol / NKSpaceCols=TW.image.NCol
% NKSpaceLins = TW.hdr.Config.NLinMeas;%TW.hdr.Config.RawLin  
% NKSpacePars = TW.hdr.Config.NParMeas;%TW.hdr.Config.RawPar  
assert(ROOverSamp*NImageCols==NKSpaceCols,'TWIX2rec:: Parameters involving over-sampling dealing not consistent.')

%Grid K-Space array size 
NKSpaceColsGrid = TW.hdr.Meas.iRoFTLength;
% NKSpaceLinsGrid = TW.hdr.Meas.iPEFTLength;
% NKSpaceParsGrid = TW.hdr.Meas.i3DFTLength;
NKSpaceLinsGrid = TW.hdr.Config.NImageLins*(1+PEOverSamp(1))*TW.hdr.MeasYaps.sKSpace.dPhaseResolution;
NKSpaceParsGrid = TW.hdr.Config.NImagePar*(1+PEOverSamp(2))*TW.hdr.MeasYaps.sKSpace.dSliceResolution;
if mod(NKSpaceLinsGrid,1)~=0; warning('KSpace array size for PE_1 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceLinsGrid,round(NKSpaceLinsGrid)); NKSpaceLinsGrid = round(NKSpaceLinsGrid);end
if mod(NKSpaceParsGrid,1)~=0; warning('KSpace array size for PE_2 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceParsGrid,round(NKSpaceParsGrid)); NKSpaceParsGrid = round(NKSpaceParsGrid);end

NImageLins = NKSpaceLinsGrid;%We don't care about the recon FOV. We reconstruct at the resol/FOV we acquire.
NImagePars = NKSpaceParsGrid;

%Padding to be performed in k-space
kSpacePad = zeros([2 3]);
%Padding at origin of array
kSpacePad(1,1) = centerIdx(NKSpaceColsGrid) - TW.image.centerCol(1);
kSpacePad(1,2) = centerIdx(NKSpaceLinsGrid) - TW.image.centerLin(1);
kSpacePad(1,3) = centerIdx(NKSpaceParsGrid) - TW.image.centerPar(1);
%Padding after array
kSpacePad(2,1) = NKSpaceColsGrid - TW.image.sqzSize(1) - kSpacePad(1,1);
kSpacePad(2,2) = NKSpaceLinsGrid - TW.image.sqzSize(3) - kSpacePad(1,2);
kSpacePad(2,3) = NKSpaceParsGrid - TW.image.sqzSize(4) - kSpacePad(1,3);
assert(all(kSpacePad(:)>=0),'Reconstruction pipeline can not yet deal with negative padding.');
applyPad = kSpacePad>0;
if any(applyPad(:,1))
    warning('TWIX2Rec:: Needed to pad k-space RO dimension with [ %d ] pre- and [ %d ] post-elements since size(TW.image) was not consistent with physically acquired Columns.',kSpacePad(1,1),kSpacePad(2,1));
end
if multDimSum(applyPad(:,2:3))>0;warning('TWIX2Rec:: Needed to pad k-space PE dimensions with [ %.0f %.0f ] pre- and [ %.0f %.0f ] post-elements since size(TW.image) was not consistent with physically acquired Lines/Partitions.',kSpacePad(1,2),kSpacePad(1,3),kSpacePad(2,2),kSpacePad(2,3));end
assymEcho = multDimSum( kSpacePad(:,1))>0;

%If padding in from of array, need to change the indices to extract
if applyPad(1,2); Lin = Lin+kSpacePad(1,2); end
if applyPad(1,3); Par = Par+kSpacePad(1,3); end
        
%Acceleration
if TW.hdr.MeasYaps.sPat.ucPATMode>1;rec.Enc.UnderSampling.isAccel = 1;else; rec.Enc.UnderSampling.isAccel = 0;end %ucPATMode: No accel (1) GRAPPA (2) or mSENSE (3)
if isfield(rec.Enc, 'UnderSampling') && rec.Enc.UnderSampling.isAccel 
    fprintf('Accelerated scan detected.\n');
    %Method of acceleration
    if TW.hdr.MeasYaps.sPat.ucPATMode==2; rec.Enc.UnderSampling.accelMethod = 'GRAPPA'; else ; rec.Enc.UnderSampling.accelMethod = 'mSENSE'; end
    %Type of method
    if isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==2
        rec.Enc.UnderSampling.accelType = 'Integrated'; 
    elseif isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==4
        rec.Enc.UnderSampling.accelType = 'Non-Integrated'; 
    else
        rec.Enc.UnderSampling.accelType = []; 
    end
    %Assign autocalibration lines to rec.ACS
    if strcmp(rec.Enc.UnderSampling.accelType,'Non-Integrated'); rec.ACS = dynInd(TW.refscan,':',ND);end 
    
    %Undersampling factors and corresponding k-space indices to extract
    rec.Enc.UnderSampling.R = cat(2, TW.hdr.MeasYaps.sPat.lAccelFactPE, TW.hdr.MeasYaps.sPat.lAccelFact3D );
    if isempty(targetR)
        rec.Enc.UnderSampling.RTarget = rec.Enc.UnderSampling.R./(1+PEOverSamp);
    else
        rec.Enc.UnderSampling.RTarget = targetR;
    end
    
    idxAcqPE = cell(1,2);
    offsetAcq{1} = rec.Enc.UnderSampling.R(1)*floor((min(Lin)-1)/rec.Enc.UnderSampling.R(1));
    offsetAcq{2} = rec.Enc.UnderSampling.R(2)*floor((min(Par)-1)/rec.Enc.UnderSampling.R(2));
    idxAcqPE{1} = min(Lin)-offsetAcq{1}:rec.Enc.UnderSampling.R(1):NKSpaceLinsGrid;
    idxAcqPE{2} = min(Par)-offsetAcq{2}:rec.Enc.UnderSampling.R(2):NKSpaceParsGrid;
    
    %Save in rec to access later
    rec.Enc.UnderSampling.limPE1 = [idxAcqPE{1}(1) idxAcqPE{1}(end)];
    rec.Enc.UnderSampling.limPE2 = [idxAcqPE{2}(1) idxAcqPE{2}(end)];
    
%     tt = zeros([NKSpaceLinsGrid, NKSpaceParsGrid]);
%     tt = dynInd(tt,idxAcqPE,1:2,1);
%     figure;imshow(tt,[]);
end

%% NOISE READOUT
typ2Rec={'y'};
if rec.Enc.UnderSampling.isAccel && strcmp(rec.Enc.UnderSampling.accelType,'Non-Integrated');typ2Rec{end+1}='ACS';end

if existsFileVar(strcat( noiseFileName,'.mat'))==1  || existsFileVar(strcat( noiseFileName,'.dat'))==1 %Take noise samples from other data
    %Get the data that should contain the noise files
    if existsFileVar(strcat( noiseFileName,'.mat'),'rec')==2 % Converted .mat file exists
        recN = load(noiseFileName,'rec'); recN=recN.rec;
    elseif existsFileVar(strcat( noiseFileName,'.dat'))==1 % Convert .dat file on the fly
        fprintf('Converting noise file on the fly: %s\n',noiseFileName)
        evalc('recN = dat2Rec(noiseFileName,[],[],1);')
    end 
    %Assign it to the current rec structure if it exists
    if isfield(recN,'N')
        fprintf('Using noise samples from separate acquisition: %s\n',noiseFileName)
        rec.N = recN.N;
        rec.PT.noiseFileName = noiseFileName;
    else
        warning('TWIX2rec:: No noise samples found in %s.mat provided/converted. Noise de-correlation disabled!',noiseFileName);%Noise field does not exist
    end
elseif isfield(TW,'noise') %Use noise for current file
    if ~isPT
        rec.N=TW.noise.unsorted;%rec.N=dynInd(TW.noise,':',ND); 
        typ2Rec=[typ2Rec 'N'];%Still need to pre-process the noise signal
    else
        warning('TWIX2rec:: Noise samples ignored since PT signal is present. No external noise file was provided.');
    end  
end

if isfield(rec,'N') && gpu;rec.N=gpuArray(rec.N);end

%% K-SPACE READOUT
y=[];%Not yet in GPU (high resolution scans would be too big with over-sampling) 
perm=1:ND;
perm(2:4)=[3:4 2]; %Re-arrange so that 1st = Readout (RO), 2nd = First PE direction (PE1 or Lines), 3rd = Second PE direction (PE2 or Partitions) and 4th = Channel
    
%removeZeroSamples=0
reverseStr = '';
for s=1:NCha
    msg = sprintf('Processing channel: %d/%d\n',s,NCha);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %Extract channel of interest
    rec.y=dynInd(TW.image,{s ':'},[2 ND]); %RO-Channel-PE1-PE2
    %Pad if needed
    if multDimSum(kSpacePad (1,:))>0; rec.y= padArrayND(rec.y, [kSpacePad(1,1) 0 kSpacePad(1,2) kSpacePad(1,3)],[],0,'pre'); end
    if multDimSum(kSpacePad (2,:))>0; rec.y= padArrayND(rec.y, [kSpacePad(2,1) 0 kSpacePad(2,2) kSpacePad(2,3)],[],0,'post'); end
    %Remove undersampling
    if removeZeroSamples && rec.Enc.UnderSampling.isAccel; rec.y = dynInd(rec.y, idxAcqPE, 3:4);end
    %Fourier transform in the readout direction
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
                plotND({angle(dynInd(rec.y,1,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .7*dynInd(abs(rec.y),centerIdx(NY(1:3)),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_3D_%d',acqName,s));
                %Plot central line of k-space
                line = dynInd(rec.y,centerIdx(NY(2:3)),2:3).';
                visPTSignal (line, [], [],[], 1, [], [], 200, replace(sprintf('\\textbf{k-space:}\n %s - Coil %d',acqName,s),'_',' '), fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_Line_%d',acqName,s),'k-space readout [\#]',{'y','a.u.'});
            end
        end
    
        for l=1
            %rec.(t2r)= fftshiftGPU(rec.(t2r),l); 
            rec.(t2r) = fftshiftOperator(rec.(t2r), 1, 0, l);
            rec.(t2r)= fftGPU(rec.(t2r),l); %YB: this has once given inf as number close to realmax('single') 
            %rec.(t2r)= ifftshiftGPU(rec.(t2r),l);%if ifft here, you need to change PT on the bottom to ifft as well (before fft)
            rec.(t2r) = fftshiftOperator(rec.(t2r), 2, 0, l);
        end  
    end
    typ2Rec = typ2Rec(strcmp(typ2Rec,'y'));%All the channels are aready in the rec.N so loop only needs to permute + SupportFOV once
    y=cat(4,y,gather(rec.y));
    
    if exportFlag && ismember(s,coilToPlot)
        if isPT; plotND({angle(dynInd(rec.y,1,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .7*dynInd(abs(rec.y),centerIdx(NY(1:3)),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{Hybrid k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilHybrid') , sprintf('%s_%d',acqName,s)); end
    end
end
rec.y=y;y=[];

%% EXTRACT PART OF READOUT TO AVOID FFT IN BOTH PE DIRECTIONS
%%% Remove PT signal
if isPT
    fprintf('Pilot Tone signal detection activated.\n');
    NY = size(rec.y);
    rec.PT.enforceUpperFOV=1
    rec.PT.sincInterp.Flag=0
    if rec.PT.sincInterp.Flag
        rec.PT.sincInterp.windowCent=.75;
        rec.PT.sincInterp.windowWidth=1/4;
        rec.PT.sincInterp.padFac=0;
    end
    %Parameters on original grid
    rec.PT.yProjRO = multDimMea( abs( dynInd(rec.y,centerIdx(NY(2:3)),2:3 )),4); %k-space centre averaged across coils        
    yProjRO = fillCent(rec.PT.yProjRO, 1, 0, rec.PT.enforceUpperFOV);%This assumes the PT is set outside the FOV (and optionally in the upper part of the over-sampled FOV)
    [~,rec.PT.idxRO] = max(yProjRO,[],1);
    rec.PT.factorFOV = ( rec.PT.idxRO - centerIdx(NY(1)) )/ ( NY(1)/2 );%If 0 -> Centre FOV / .5 -> Edge FOV / 1 -> Edge oversampled FOV
    dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9; BWVox= 1 / dwellTime/ NImageCols;%In Hz
    try rec.PT.FWHMHz = getFWHM( BWVox * (1:length(yProjRO)) , yProjRO');catch;rec.PT.FWHMHz = 0;end%Measure FWHM        
    
    %Extract
    if rec.PT.sincInterp.Flag
        rec.PT.pSliceImage = interpolatePT(rec.PT.pSliceImage, 1, rec.PT.sincInterp.padFac,  rec.PT.sincInterp.windowCent, rec.PT.sincInterp.windowWidth, [], rec.PT.enforceUpperFOV);
        rec.PT.offsetMBHz=0;
        rec.PT.idxMB=1;
    else%PT signal can be multiband signal       
        rec.PT.offsetMBHz=0%Multiband in Hz to extract 
        if isempty(rec.PT.offsetMBHz); rec.PT.offsetMBHz = 2*rec.PT.FWHMHz;end%By default twice the FWHM
        rec.PT.offsetMBidx=round(rec.PT.offsetMBHz / BWVox)*[ 1 1];
        rec.PT.offsetMBidx= [ min(rec.PT.offsetMBidx(1),rec.PT.idxRO-1)   min(rec.PT.offsetMBidx(2),-rec.PT.idxRO+NY(1))   ];     
        rec.PT.pSliceImage = dynInd(rec.y, rec.PT.idxRO-rec.PT.offsetMBidx(1):rec.PT.idxRO+rec.PT.offsetMBidx(2),1);%Extract from rec.y (image domain), not from yTemp
        rec.PT.idxMB = rec.PT.offsetMBidx(1)+1;%round(rec.PT.offsetMBHz / BWVox) +1;   
    end
    %Average Multiband frequencies
    rec.PT.isAveragedMB = 0;
    if rec.PT.isAveragedMB && rec.PT.sincInterp.Flag==0
        rec.PT.pSliceImage = multDimMea(rec.PT.pSliceImage,1);%Note this is averaging in RO direction
        rec.PT.idxMB=1;
    end
    
    %Reference phase to avoid phase drift
    rec.PT.isRelativePhase = 0
    rec.PT.relativePhaseCoilIdx = 1;
    if rec.PT.isRelativePhase%This in hybrid k-space, as it is supposed to be for phase subtraction
        rec.PT.pSliceImage = rec.PT.pSliceImage .* conj( exp(1i*angle(dynInd(rec.PT.pSliceImage,rec.PT.relativePhaseCoilIdx,4)))); 
    end
end   

%%% Remove over-sampling
if removeOversampling
    fprintf('Removing oversampling.\n');
    rec.y=resampling(rec.y,NImageCols,2);%Remove over-encoding provided by Siemens (used e.g. in Pilot Tone) 
    if isfield(rec,'N');rec.N=resampling(rec.N,NImageCols,2);end
    if isfield(rec,'ACS');rec.ACS=resampling(rec.ACS,NImageCols,2);end
end
NImageColsArray=size(rec.y,1);%The FOV array kept in the readout direction

%%% Reduce FOV in readout dimension
if ~isempty(suppFOV)
    fprintf('Extracting part of FOV with range [%.2f --> %.2f].\n', suppFOV(1), suppFOV(2));
    NY = size(rec.y);
    vr=max(1,round(NY(1)*suppFOV(1))):min( NY(1),round(NY(1)*suppFOV(2))); %Only in the RO direction
    rec.y=dynInd(rec.y,vr,1); %Extract supported FOV
    if isfield(rec,'N');rec.N=dynInd(rec.N,vr,1);end %No need to remove oversampling for noise I think 
    if isfield(rec,'ACS');rec.ACS=dynInd(rec.ACS,vr,1);end
else
    vr=[];
end
 
%% TRANSFORM TO IMAGE DOMAIN IN PE DIMENSIONS
%%% Check if gpu can handle all channels
yMem = getSize(rec,1);%In bytes
if gpu
    ss=gpuDevice; 
    safetyFact=.85;%Safety margin for GPU memory
    if yMem > (safetyFact * ss.AvailableMemory); gpu=0; rec=gatherStruct(rec); end%Disable gpu and take everything out of it
end
if gpu;rec.y=gpuArray(rec.y);end

%%% Fourier transforms in PE dimensions - move all channels to image domain
NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    %rec.y=fftshiftGPU(rec.y,n);%This affects sampling in k-space
    rec.y = fftshiftOperator(rec.y, 1, 0, n);
    rec.y=fftGPU(rec.y,n);
    %rec.y=ifftshiftGPU(rec.y,n);%This determines where image is --> mapping S onto it will matter
    rec.y = fftshiftOperator(rec.y, 2, 0, n);
    
    if isPT
        %rec.PT.pSliceImage=fftshiftGPU(rec.PT.pSliceImage,n);
        rec.PT.pSliceImage = fftshiftOperator(rec.PT.pSliceImage, 1, 0, n);
        rec.PT.pSliceImage=fftGPU(rec.PT.pSliceImage,n);
        %rec.PT.pSliceImage=ifftshiftGPU(rec.PT.pSliceImage,n);
        rec.PT.pSliceImage = fftshiftOperator(rec.PT.pSliceImage, 2, 0, n);
    end
    if isfield(rec,'ACS')
        %rec.ACS=fftshiftGPU(rec.ACS,n);
        rec.ACS = fftshiftOperator(rec.ACS, 1, 0, n);
        rec.ACS=fftGPU(rec.ACS,n);
        %rec.ACS=ifftshiftGPU(rec.ACS,n);
        rec.ACS = fftshiftOperator(rec.ACS, 2, 0, n);
    end
end

%%% Plot coil images to inspect later
if exportFlag
    if ~exist(fullfile( exportFolder ,'CoilImage'),'dir'); mkdir(fullfile(exportFolder ,'CoilImage'));end
    for s=1:coilToPlot
        plotND({abs(dynInd(rec.y,{s,1},[4 8])),1,0},angle(dynInd(rec.y,{s,1},[4 8])),[-pi pi 0 2.5*multDimMea(abs(dynInd(rec.y,{s,1},[4 8])))],[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{Image:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilImage') , sprintf('3D_%s_%d',acqName,s));
    end
    plotND([],RSOS(dynInd(rec.y,1,8)),[0 .9*multDimMax(RSOS(dynInd(rec.y,1,8)))],[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{RSOS:} %s',acqName),'_',' '),[],[],fullfile(exportFolder ,'RSOS') , sprintf('%s',acqName));
end

%% PLOT PT SIGNAL
if isPT
    line = multDimSum(abs(rec.PT.yProjRO),[2:4]).';
    peakId = zerosL(line); peakId(rec.PT.idxRO)=1;
    visPTSignal (line, [], [],[], 1, peakId, [], 199, replace(sprintf('\\textbf{Central line hybrid k-space (mean across coils):}\n %s',acqName),'_',' '), [], [], 'Image readout [\# voxels]',{'y_{hybrid}','a.u.'});
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Peak',acqName)); end
    
    if gpu;rec.PT=gatherStruct(rec.PT);end
end

%% DE-CORRELATE NOISE
if isfield(rec,'N')
    fprintf('De-correlating channels from noise samples.\n');
    rec.N = dynInd(rec.N,1:NCha,4);
    [rec.y, rec.Alg.covMatrix] = standardizeCoils(rec.y,rec.N);
    if isfield(rec,'ACS');[rec.ACS] = standardizeCoils(rec.ACS,rec.N);end
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
rec.Enc.FOVmm = [sliceInfo.dReadoutFOV , sliceInfo.dPhaseFOV , sliceInfo.dThickness];
rec.Enc.FOVmm = rec.Enc.FOVmm .* (1+ [0 PEOverSamp]); %Recale the FOV that is returned (not including the over-sampling)
rec.Enc.FOVSize = [NY(1), NImageLins, NImagePars]; %Work on the FOV of the acquired (over-sampled) data

rec.Enc.AcqVoxelSize=rec.Enc.FOVmm ./ rec.Enc.FOVSize;
rec.Enc.AcqSize=[ 2^(ROOverSamp==2)*NY(1) NY(2:3)];

%%% CHANGE TO PE-RO-SL                  
rec.y=gather(rec.y);if isfield(rec,'N');rec.N=gather(rec.N);end
perm=1:ND; perm(1:3)=[2 1 3];
rec.y = permute(rec.y,perm);
if isfield(rec,'ACS'); rec.ACS = permute(rec.ACS,perm); end
if isfield( rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
rec.Enc.AcqVoxelSize = rec.Enc.AcqVoxelSize(perm(1:3));%Need to change spacing as well

%% GEOMETRY COMPUTATION FOR RECONSTRUCTION ARRAY
fprintf('Computing geometry for reconstruction array.\n');
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

%Account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = [ NImageLins, inf ,NImagePars];%Set inf to make sure this element is replaced
N(2) = NImageColsArray;%N(2) will depend if you removed over-sampling
orig = ( centerIdx(N(1:3)) - [0 0 .5] )';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
%The [0 0 .5] is also the difference you find between the gadgetron PACS exported data and Siemens PACS exported data
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
[~,rec.Par.Mine.APhiRec] = mapNIIGeom([], rec.Par.Mine.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;%Backup as how was read in

%DEDUCE ACQUISITION ORDER
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane ] = acquisitionOrder(rec.Par.Mine.APhiRec);
rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
fprintf('Slice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
fprintf('    Readout: %s\n', rec.Par.Scan.Mine.RO );
fprintf('    1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
fprintf('    2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

%% GEOMETRY COMPUTATION FOR AUTOCALIBRATION ARRAY
if isfield(rec,'ACS')
    fprintf('Computing geometry for autocalibration array.\n');
    YSize = [ NImageLins, NImageColsArray ,NImagePars];%Set inf to make sure this element is replaced
    ACSSize = [size(rec.ACS,1), size(rec.ACS,2), size(rec.ACS,3)];
    [rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS] = mapNIIGeom([], rec.Par.Mine.APhiRec,'resampling',[],YSize,ACSSize);%TO CHECK:AD HOC - REVERSE ENGINEERED
end

%% MAKE RO-PE-SL FOR ALIGNED SENSE PIPELINE CONVENTION
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
fprintf('Reading sequence parameters.\n');

rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;%In ms
rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;%In ms
rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};%In degrees
if isfield(TW.hdr.MeasYaps.sFastImaging,'ulEnableRFSpoiling')
    if strcmp(TW.hdr.MeasYaps.sFastImaging.ulEnableRFSpoiling,'0x1');rec.Par.Labels.RFSpoiling=1;%Hexadecimal code
    else; rec.Par.Labels.RFSpoiling = TW.hdr.MeasYaps.sFastImaging.ulEnableRFSpoiling;
    end
else
    rec.Par.Labels.RFSpoiling=0;
end
rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;%In seconds
rec.Par.Labels.dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;%In seconds %*2 /1e+9 since alDwellTime{1} in nanoseconds and for oversampling
rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz/FOV
rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / NImageCols;%In Hz/voxel in the readout direction
rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm

rec.Par.Labels.fieldStrength = TW.hdr.Dicom.flMagneticFieldStrength; %In T

rec.Par.Labels.sequenceType = TW.hdr.Meas.tScanningSequence;
if strcmp(rec.Par.Labels.sequenceType,'GR\IR')
    isShotBased = 1;% TW.hdr.MeasYaps.sFastImaging.lTurboFactor~=1;
    if TW.hdr.MeasYaps.sKSpace.unReordering==1 %Linear --> Each shot covers all partition samples
        rec.Par.Labels.shotLayout='Linear';
        rec.Par.Labels.TFEfactor = size(rec.y,3);
    elseif TW.hdr.MeasYaps.sKSpace.unReordering==256 %Linear rotated --> Each shot covers all line samples
        rec.Par.Labels.shotLayout='Linear-Rotated';
        rec.Par.Labels.TFEfactor = size(rec.y,2);%TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor  
    end
else
   isShotBased=0; 
end
if isShotBased
    %rec.Par.Labels.TFEfactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor  
    rec.Par.Labels.NShots = length(TW.image.Par)/rec.Par.Labels.TFEfactor;
    fprintf('Shot-based sequence with %d shots (%d samples per shot).\n',rec.Par.Labels.NShots, rec.Par.Labels.TFEfactor)
    assert(mod(NShots,1)==0,'TWIX2Rec:: Number of samples per shot not an integer. Probably a corrupted scan.')
else
    rec.Par.Labels.TFEfactor=length(Par); 
end
rec.Par.Labels.ZReconLength=1;

%% SAMPLING INFORMATION
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

%%% Flip sampling and offset (probably because ambiguous Siemens fft vs. ifft convention)
offset =-1;
Lines = mod( Lines-1 + offset, NY(2) ) +1;
Partitions = mod( Partitions -1 + offset, NY(3) ) +1;

%%% Shift to make compatible with DISORDER recon
kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE 
rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices
%ATTENTION: this mod() operator gives big spikes in the sampling. This is just a processing implication.

rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};

%%%DISORDER information
if isfield(TW.hdr.MeasYaps,'sWipMemBlock') && isfield(TW.hdr.MeasYaps.sWipMemBlock,'alFree') && TW.hdr.MeasYaps.sWipMemBlock.alFree{1}==2
    rec.Enc.DISORDER=1;
    if ~isShotBased %DISORDER parameters not reliable when having a shot based sequence (To be fixed)
        rec.Enc.DISORDERInfo.tileSize=cat(2,TW.hdr.MeasYaps.sWipMemBlock.alFree{2},TW.hdr.MeasYaps.sWipMemBlock.alFree{3});
        rec.Enc.DISORDERInfo.NShots=prod(rec.Enc.DISORDERInfo.tileSize);
        rec.Enc.DISORDERInfo.segmentDuration = size(rec.Assign.z{2},2)/rec.Enc.DISORDERInfo.NShots*rec.Par.Labels.RepetitionTime/1000;%In s - Use rec.Assign.z since NY does not take into account shutter
        %if rec.Enc.DISORDER; assert(rec.Par.Labels.NShots==rec.Enc.DISORDERInfo.NShots,'DISORDER encoded shot-based sequence should have compatible shot/segment durations.');end
    end
else
    rec.Enc.DISORDER=0;
end

%%% Shutter
rec.Enc.ellipticalShutter = ~isempty(TW.hdr.Meas.ucEnableEllipticalScanning) && TW.hdr.Meas.ucEnableEllipticalScanning==1;

%%% Plot trajectory
permTraj=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
kTraj=zeros([length(rec.Assign.z{2}) 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{permTraj(n)}(:);end
kTrajSS=reshape(kTraj,[length(rec.Assign.z{2}) 1 2]);
visTrajectory(kTrajSS,0,[],[],206);%h=figure('color','w'); imshow(timeMat , []);title('Time indices of k-space sampling'); set(h,'color','w','Position',get(0,'ScreenSize'));
sgtitle(replace(sprintf('\\textbf{Trajectory:} %s\n',acqName),'_',' '),'Interpreter','latex','Color',[1 1 1]*0,'FontSize',20);
if exportFlag; saveFig(fullfile(exportFolder ,'Trajectory',acqName)); end

%% SHIM SETTINGS
fprintf('Reading shim information.\n');
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
    if existsFileVar(correctionFactorsFileName)==1
        ss = load(correctionFactorsFileName); rec.Par.Labels.Shim.correctFact_mA2au = ss.correctFact_mA2au;%Multiplicative correction
        rec.Par.Labels.Shim.shimCurrents_mA = rec.Par.Labels.Shim.shimCurrents_au ./ rec.Par.Labels.Shim.correctFact_mA2au;
    else
       fprintf('No Calibration data (correctionFactrs.mat) found in the shim directory.\n')  
    end

    crossTermFileName = strcat(shimDir,'Calibration', filesep, 'Results',filesep, 'crossTerms.mat');
    if existsFileVar(crossTermFileName)==1 
        ss = load(crossTermFileName); rec.Par.Labels.Shim.crossTerms = ss.crossTerms;
    else
       fprintf('No Calibration data (crossTerms.mat) found in the shim directory.\n')  
    end
else
   fprintf('No shim directory provided and hence no Calibration data read out.\n') 
end

%% PREPARE TEMPORAL PT SIGNAL (OPTIONAL)
if isPT
    tt = rec.PT.pSliceImage;
    
    %[kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    %filterStruct =[];filterStruct.medFiltKernelWidthidx=400/ rec.Par.Labels.RepetitionTime; 
    %tt = filterPTSlice (tt, filterStruct, NY, kIndex);

    for n=2:3
        tt=fftshiftOperator(tt,1,1,n);
        tt=fftGPU(tt,n)/NY(n);
        tt=fftshiftGPU(tt,n);%fftshift since in solveXT there is an iffthift on timeIndex
    end

    [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);

    tt = resSub(tt, 2:3);
    rec.PT.pTimeTest = permute(dynInd(tt, idx, 2),[3 2 1]);%sorted
    
    coilPlot=1:min(5,size(rec.PT.pTimeTest,1));
    plotType = 1;% (1) magn/phase (2) real/image (3) all
    %Plot original PT signal
    visPTSignal(dynInd(rec.PT.pTimeTest,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 200, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',acqName)) ; end
    visPTSignal(normm(rec.PT.pTimeTest(:,:,rec.PT.idxMB),[],1), [], [],[], plotType, [], [], 201, replace(sprintf('\\textbf{Norm original PT signal:}\n %s',acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',strcat(acqName,'_Norm'))) ; end
    %visPTSignalImage (rec.PT.pTimeTest(:,1:100,rec.PT.idxMB), [], [],[], plotType, [], [], 202, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    %if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',strcat(acqName,'_Image'))) ; end

    %Pre-processed signal
    tt = rec.PT.pTimeTest;
    tt = bsxfun(@rdivide, tt , sqrt(normm(dynInd(rec.PT.pTimeTest,{':',rec.PT.idxMB},[1 3]),[],1)) );
    if rec.PT.isRelativePhase; tt = bsxfun(@minus, tt, multDimMea((tt),2) );end%This only works well if you have referenced the phase to a channel
    visPTSignal (dynInd(tt,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 203, replace(sprintf('\\textbf{Pre-processed PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','TracePreProcessed',acqName)) ; end
    %visPTSignalImage (tt(:,1:100,rec.PT.idxMB), [], [],[], plotType, [], [], 205, replace(sprintf('\\textbf{Pre-processed PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    %if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','TracePreProcessed',strcat(acqName,'_Image'))) ; end
end

%% ADD ADDITIONAL PT SIGNAL INFORMATION (OPTIONAL)
% if isPT
%     centreFOV = rec.Par.Mine.APhiRec * [ceil( (rec.Enc.FOVSize(1:3)-1 -.5 +1)/2) 1]';%-1 since NIFTI convention and not MATLAB 
%     centreFOV = rec.Par.Mine.APhiRec * [(rec.Enc.FOVSize(1:3)/2 -[ 0 0 .5]) 1]';%-1 since NIFTI convention and not MATLAB
%     RODir = rec.Par.Mine.APhiRec * [1;0;0;1] - rec.Par.Mine.APhiRec * [0;0;0;1];
% 
%     offsetFOV = multDimSum(RODir(1:3).*centreFOV(1:3))./sqrt(normm(RODir(1:3)));%Projection on RO direction
% 
%     rec.PT.transmitFreqHz = setPTFreq(...
%        rec.Par.Labels.Shim.shimCurrents_au(1) ,...
%        NImageCols,...
%        rec.Par.Labels.voxBandwidth ,...
%        rec.PT.factorFOV, ...
%        offsetFOV,...
%        rec.Enc.AcqVoxelSize(1))/1e6;
% end

%% MAKE UNDERSAMPLING APPROPRIATE FOR MOTION CORRECTION
if isfield(rec.Enc,'UnderSampling') && isfield(rec.Enc.UnderSampling,'RTarget') && ~isequal(rec.Enc.UnderSampling.RTarget,rec.Enc.UnderSampling.R) && removeZeroSamples
    %Store original parameters if ever needed
    rec.Enc.UnderSampling.Acq.FOVSize = rec.Enc.FOVSize;
    rec.Enc.UnderSampling.Acq.AcqSize = rec.Enc.AcqSize;
    rec.Enc.UnderSampling.Acq.R = rec.Enc.UnderSampling.R;
    rec.Enc.UnderSampling.Acq.FOVmm = rec.Enc.FOVmm;
    rec.Enc.UnderSampling.Acq.PEOverSamp = PEOverSamp;
    %Assign new values
    rec.Enc.FOVSize(2:3) = round(rec.Enc.UnderSampling.RTarget .* rec.Enc.AcqSize(2:3));
    rec.Enc.FOVmm = rec.Enc.FOVSize .* rec.Enc.AcqVoxelSize;
    rec.Enc.UnderSampling.R = rec.Enc.FOVSize(2:3)./rec.Enc.AcqSize(2:3);
    %Change orientation information accordingly
    vrLin=(centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(2))  - centerIdx(rec.Enc.FOVSize(2)))+[1:rec.Enc.FOVSize(2)];
    vrPar=(centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(3))  - centerIdx(rec.Enc.FOVSize(3)))+[1:rec.Enc.FOVSize(3)];
    [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'dynInd',{':',vrLin,vrPar},rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'ipermute',rec.Par.Mine.permuteHist{1});
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'dynInd',{':',vrLin,vrPar},rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'permute',rec.Par.Mine.permuteHist{1});
    fprintf('Changing undersampling from R=%.1fx%.1f to R=%.2f-%.2f by changing the FOV size.\n', rec.Enc.UnderSampling.Acq.R(1),rec.Enc.UnderSampling.Acq.R(2),rec.Enc.UnderSampling.R(1),rec.Enc.UnderSampling.R(2));
end

%% RESAMPLE STRUCTURE TO OTHER RESOLUTION/UNDER-SAMPLING
if ~isempty(resRec)
    if length(resRec)==1; resRec = resRec * ones([1 3]); else; resRec(end+1:3) = rec.Enc.AcqVoxelSize((length(resRec)+1):3);end
    fprintf('Resampling reconstruction to new resolutions %.2fmm (RO) / %.2fmm (PE1) %.2fmm (PE2).\n', resRec)
    rec = resampleRec(rec,resRec); %Resamples the arrays but also changes the sampling pattern accordingly
end

%% VARIA
%rec.Varia.ECC = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1,1};
%rec.Varia.hdr = TW.hdr;

%% CLOSE
rec=gatherStruct(rec);


%%% RESOURCES
%FFT array sizes
%https://www.magnetom.net/t/number-of-lines-partitions-vs-fftlength/3012
%https://www.magnetom.net/t/anisotropic-in-plane-resolution-in-ima-header/797

%Bandwidth
%https://www.magnetom.net/t/how-to-get-the-spectroscopy-bandwidth-value-from-meas-as/320/4
%https://www.magnetom.net/t/dwelltime-in-header-dat-and-rda-files-incorrect-for-spectroscopy/3342

%Geometry
% asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1}; %See invert7T_old to see how to use this
% sNormal=asSlice.sNormal; sPosition=asSlice.sPosition;