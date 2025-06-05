
function rec = TWIX2recv1(TW,suppFOV,removeReadoutDelay, shimDir, removeOversampling)

%TWIX2REC   Inverts SIEMENS TWIX data to a reconstruction structure.
%   * TW is the TWIX object
%   * {SUPPFOV} is the support in the readout direction for high resolution scans
%   * {REMOVEREEADOUTDELAY} indicated to shift center of k-space due to readout delay (defaults to 0).
%   * {SHIMDIR} is directory with additional shim information of the scanner
%   ** REC is a reconstruction structure
%

ND=16; %Number of dimensions that can be called in the TWIX object (see public method subsref)
TW.image.flagIgnoreSeg=1;
if nargin<2 || isempty(suppFOV);suppFOV=[];end
if nargin<3 || isempty(removeReadoutDelay);removeReadoutDelay=0;end
if nargin<4 || isempty(shimDir);shimDir='/home/ybr19/Projects/B0Shimming/';end
if nargin<5 || isempty(removeOversampling);removeOversampling=1;end

gpu=(gpuDeviceCount>0 && ~blockGPU);

%% DATA READOUT
typ2Rec={'y'};
if isfield(TW,'noise')
    typ2Rec=[typ2Rec 'N'];
    rec.N=TW.noise.unsorted;%rec.N=dynInd(TW.noise,':',ND); 
    %ATTENTION: This is in k-space (I assumed generated without gradients
    %playing). Since no gradients played, the PT signal will appear as a big DC
    %component on top of the white noise.
    if gpu; rec.N=gpuArray(rec.N);end
end

NCha=TW.image.NCha;
NImageCols=TW.hdr.Meas.NImageCols;

y=[];
if ~isempty(suppFOV);vr=round(NImageCols*suppFOV(1)):round(NImageCols*suppFOV(2));else vr=[];end

perm=1:ND;
perm(2:4)=[3:4 2]; %Re-arrange so that 1st = Readout, 2nd = First PE direction, 3rd = Second PE direction and 4th = Channel
    
reverseStr = '';

for s=1:NCha
    msg = sprintf('Processing channel: %d/%d\n',s,NCha);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    rec.y=dynInd(TW.image,{s ':'},[2 ND]); %1st: Readout, 2nd: Channel, 3rd: First PE direction, 4rd: Second PE direction
    if gpu;rec.y=gpuArray(rec.y);end 

    %plotND([],permute(rec.y,perm),[0 0.005],[],0,{[],2},[],[],[],[],12);
    
    for n=1:length(typ2Rec); t2r=typ2Rec{n};
        rec.(t2r)=permute(rec.(t2r),perm);
        NY=size(rec.(t2r));NY(end+1:ND)=1;    
        
        for l=1
            rec.(t2r)= fftshiftGPU(rec.(t2r),l); 
            rec.(t2r)= fftGPU(rec.(t2r),l);     %*NY(l)  YB: this has once given inf as number close to realmax('single') 
            rec.(t2r)= ifftshiftGPU(rec.(t2r),l);
        end
        
    end
    typ2Rec={'y'}; % YB: all the channels are aready in the rec.N so loop only needs to permute + SupportFOV once
    y=cat(4,y,gather(rec.y));
end
rec.y=y;y=[];if gpu;rec.y=gpuArray(rec.y);end

NY=size(rec.y);NY(end+1:ND)=1;
NYOrig = NY;
for n=2:3
    rec.y=fftshiftGPU(rec.y,n);
    rec.y=fftGPU(rec.y,n);%*NY(l);
    rec.y=ifftshiftGPU(rec.y,n);
end%Now all dimensions in image domain

%Apply shift here so that we don't need to circshift and we can grab samples out of over-sampled area
if removeReadoutDelay; rec.y = rec.y .* conj( multDimRidge(rec.y, 1:4, 1, 1));end%Shift k-space to correct for readout delay  

%%% Extract PT signal
    %Detect peak outside FOV
    yProjRO =  multDimSum(abs(rec.y),[2:4]); %Not in the readout
        zeroF=ceil((NY+1)/2);
        orig=zeroF-ceil(((NY/2)-1)/2); %YB: These are the elements to extract from the bigger F matrix to resample the array
        fina=zeroF+floor(((NY/2)-1)/2);
        idxTemp = orig(1):fina(1);
    yProjRO = dynInd(yProjRO, idxTemp, 1, 0);%Set FOV to zero. This assumes the PT is set outside the FOV
    [~,rec.PT.idxRO] = max(yProjRO,[],1) ;
    rec.PT.factorFOV = ( rec.PT.idxRO - NY(1)/2 )/ ( NY(1)/2 );%If 0 -> Centre FOV / .5 -> Edge FOV / 1 -> Edge oversampled FOV

    %Extract PT signal (can be multiband signal)
    rec.PT.offsetMB=0;
    rec.PT.idxMB = rec.PT.offsetMB +1;
    rec.PT.pSliceImage = dynInd(rec.y, (rec.PT.idxRO-rec.PT.offsetMB):(rec.PT.idxRO+rec.PT.offsetMB),1);%Extract from rec.y (image domain), not from yTemp

    %Average Multiband frequencies
    rec.PT.averagedMB = 0;
    if rec.PT.averagedMB
        rec.PT.pSliceImage = multDimMea( rec.PT.pSliceImage ,1);%Note this is averaging in RO direction
        rec.PT.idxMB=1;
    end

    figure('color','w'); 

    subplot 121; plot(1:size(yProjRO,1), multDimSum(abs(rec.y),[2:4])');hold on; scatter(idxPT,yProjRO(idxPT));axis tight
    title('Untouched hybrid k-space')

    subplot 122; plot(1:size(yProjRO,1), yProjRO');hold on; scatter(idxPT,yProjRO(idxPT));axis tight
    title('FOV set to 0')

    xlabel('Readout number [#]')
    ylabel('Signal [a.u.]')
    sgtitle('Pilot Tone signal extraction (detected peak in red)')

%%% Remove over-sampling
if removeOversampling
    rec.y=resampling(rec.y,NImageCols,2);%Remove over-encoding provided by Siemens (used e.g. in Pilot Tone) 
    if isfield(rec,'N');rec.N=resampling(rec.N,NImageCols,2);end
end

%%% Reduce FOV in readout dimension
if ~isempty(vr)
    rec.y=dynInd(rec.y,vr,1); %Extract supported FOV
    if isfield(rec,'N');rec.N=dynInd(rec.N,vr,1);end 
end

%%% De-correlate noise
decorrNoise = 0;
if isfield(rec,'N') && decorrNoise;[rec.y, covMatrix] = standardizeCoils(rec.y,rec.N);else; warning('TWIX2rec:: No noise-decorrelation applied!. Might result in sub-optimal performance.');end
if exist('covMatrix','var'); figure('color','w'); imshow(covMatrix,[]);title('Noise covariance matrix.');end
    
%%% DIMENSIONS IN IMAGE DOMAIN: RO-PE-SL
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.AcqVoxelSize=[TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV/NImageCols ...%RO
                      TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV/TW.hdr.Meas.NImageLins...%PE
                      TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness/TW.hdr.Meas.NImagePar];%SL

%%% CONVERT TO PRS (PE-RO-SL), CONSISTENT WITH SIEMENS ORIENTATION CONVENTION
rec.y=gather(rec.y);if isfield(rec,'N');rec.N=gather(rec.N);end
perm=1:ND; perm(1:3)=[2 1 3];
rec.y = permute(rec.y,perm);
rec.PT = permute(rec.PT, perm);
rec.Enc.AcqVoxelSize = rec.Enc.AcqVoxelSize(perm(1:3));%need to change spacing as well

%% GEOMETRY COMPUTATION
% asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1}; %See invert7T_old to see how to use this
% sNormal=asSlice.sNormal; sPosition=asSlice.sPosition;
slicePos = TW.image.slicePos;

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

%ROTATION 
quaternionRaw = slicePos(4:7,1);
rec.Par.Mine.Arot = eye(4);
%rec.Par.Mine.Arot(1:3,1:3) = quaternionToR(quaternionRaw); %PE-RO-SL to PCS
rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom(ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS

%TRANSLATION
translationRaw = slicePos(1:3,1); %in mm for center FOV (I think)
rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

% account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = size(rec.y); N=N(1:3);
N(2) = NImageCols;%Safety for when you don't remove oversampling
orig = ( N(1:3)/2 - [0 0 .5] )';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
if ~isempty(vr); orig(2)=orig(2)-(vr(1)-1); end%YB:orig(2) since this is readout/only vr(1)-1 elements removed from array/minus sign since should be addition and below already negative sign

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)...
                        - rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;
                    
%COMBINED MATRIX
rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

%Alternative way (works): 
% [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'translate',orig);
% if ~isempty(vr);[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'dynInd',{':',vr,':'},NYOrig(1:3),NY(1:3));end

%MOVE TO RAS
rec.Par.Mine.PCS2RAS = getPCS2RAS(TW.hdr.Dicom.tPatientPosition);
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
%At this point rec.Par.Mine.APhiRec defined from PRS to RAS

rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;

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
rec.y = permute(rec.y,perm);
rec.PT = permute(rec.PT,perm);
[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'permute',perm);

%%% Store the permutations performed from the PRS space
rec.Par.Mine.permuteHist = [];rec.Par.Mine.permuteHist{1} = perm(1:4);

%%% Make MPS consisitent with RO-PE-SL
MPStemp = rec.Par.Scan.MPS ;
rec.Par.Scan.MPS(1:2) = MPStemp(4:5);
rec.Par.Scan.MPS(4:5) = MPStemp(1:2);

%% SEQUENCE INFORMATION
rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
rec.Par.Labels.ZReconLength=1;

rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;
rec.Par.Labels.dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;%In s %*2 /1e+9 since alDwellTime{1} in nanoseconds and for oversampling
%https://www.magnetom.net/t/dwelltime-in-header-dat-and-rda-files-incorrect-for-spectroscopy/3342
%https://www.magnetom.net/t/how-to-get-the-spectroscopy-bandwidth-value-from-meas-as/320/4
rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz
rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / NImageCols;%In Hz/voxel in the readout direction
rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm
NY=size(rec.y);NY=NY(1:3);

%% SAMPLING INFORMATION
%%% Flip sampling and offset (probably because ambiguous Siemens fft vs. ifft convention)
offset = -1;
Lines = mod( TW.image.Lin-1 + offset, NY(2) ) +1;
Partitions = mod( TW.image.Par -1 +offset, NY(3) ) +1;

%%% Shift to make compatible with DISORDER recon
kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE  %Minus comes from convention in invert7T.m
rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices
 
%% Make the PT signal in time to later check in solveXTB (optional)
        tt = rec.PT.pSliceImage;
        for n=2:3
            %tt=fftshiftGPU(tt,n);
            tt=fftGPU(tt,n);%*NY(l);
            tt=ifftshiftGPU(tt,n);
        end

        [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
        
        tt = resSub(tt, 2:3);
        rec.PT.pTimeTest = squeeze(dynInd(tt, idx, 2));%sorted

        %ttP =  medfilt1(abs(gather(tt)),1,[],2);
        %figure('color','w');plot(1:size(tt,1),abs(ttP(:,5)));
        %xlabel('Time (PE readouts) [TR]'); ylabel('Pilot Tone signal [a.u]');
        %axis tight

%% SHIM SETTINGS
%%% Extract shim currents
rec.Par.Labels.Shim.shimCurrents_au = zeros(13,1);
if isfield(TW.hdr.Dicom, 'lFrequency'); rec.Par.Labels.Shim.shimCurrents_au(1) = TW.hdr.Dicom.lFrequency;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetX'); rec.Par.Labels.Shim.shimCurrents_au(2) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetX;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetY'); rec.Par.Labels.Shim.shimCurrents_au(3) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetY;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetZ'); rec.Par.Labels.Shim.shimCurrents_au(4) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetZ;end
if isfield(TW.hdr.Meas, 'alShimCurrent'); rec.Par.Labels.Shim.shimCurrents_au(5:13) = TW.hdr.Meas.alShimCurrent(1:9).';end
                                  
rec.Par.Labels.Shim.shimOrder = 3; %Given the 13 parameters in the 7T Terra console                           

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



