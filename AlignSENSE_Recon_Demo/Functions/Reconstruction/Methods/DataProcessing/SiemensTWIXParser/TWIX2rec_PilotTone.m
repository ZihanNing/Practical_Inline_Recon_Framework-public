function rec=invert7T_PilotTone(TW,suppFOV)

%INVERT7T   Inverts SIEMENS data
%   * TW is the TWIX object
%   * SUPPFOV is the support in the readout direction for high resolution
%   scans
%   ** REC is a reconstruction structure
%

ND=16; % YB: Number of dimensions that can be called in the TWIX object (see public method subsref)
TW.image.flagIgnoreSeg=1;
if nargin<2;suppFOV=[];end

gpu=(gpuDeviceCount>0 && ~blockGPU);

%% DATA READOUT
typ2Rec={'y'};
if isfield(TW,'noise')
    typ2Rec=[typ2Rec 'N'];
    %rec.N=dynInd(TW.noise,':',ND);
    rec.N=TW.noise.unsorted;
    if gpu; rec.N=gpuArray(rec.N);end
end

NCha=TW.image.NCha;
NImageCols=TW.hdr.Meas.NImageCols;

y=[];
if ~isempty(suppFOV);vr=round(NImageCols*suppFOV(1)):round(NImageCols*suppFOV(2));else vr=[];end

for s=1:NCha
    rec.y=dynInd(TW.image,{s ':'},[2 ND]);
    if gpu;rec.y=gpuArray(rec.y);end %YB: Here, 1st = Readout, 2nd =  channel, 3rd = First PE direction and 4rd = Second PE direction

    perm=1:ND;
    perm(2:4)=[3:4 2]; %YB: Re-arrange so that 1st = Readout, 2nd = First PE direction, 3rd = Second PE direction and 4th = Channel
    for n=1:length(typ2Rec); t2r=typ2Rec{n};
        rec.(t2r)=permute(rec.(t2r),perm);
        NY=size(rec.(t2r));NY(end+1:ND)=1;       
        
        %%% PILOT TONE CORRECTION
        fprintf('Coil %d\n',s);
           th1=1e+37;
           if  multDimSum(isinf(th1*rec.(t2r)))>0
                fprintf('high values detected in %s\n',t2r);
                rec.(t2r)(isinf(abs(rec.(t2r)))) = 1e+30;%just high
                Ypre = rec.(t2r);
                 
                %peak identification
                h = ones([1,8,8])/((8)^2-1);h(:,4,4)=0;%81-1 since center value = 0
                th = 20;
                meanPE = imfilter(abs(rec.(t2r)),h); 
                idx = (abs(rec.(t2r)) > th * meanPE);
                Yidx = reshape(idx,size(rec.(t2r)));
                plotND([], abs(( cat(5, Yidx, rec.(t2r)))),[0 1e-3],[],0,{[],2},[],[],[],s);
                
                %interpolation kernel %interpolate neighbourhing samples
                h_interpolation = ones([1,4,4])/(4^2-1);h_interpolation(:,2,2)=0;
                linePRE = dynInd(rec.(t2r), [51, 81],2:3);
                meanPE_interpolation = imfilter(rec.(t2r),h_interpolation);

                rec.(t2r) (idx) =  meanPE_interpolation ( idx ) ;
                linePOST  = dynInd(rec.(t2r), [51, 81],2:3);

                %figure;
                %subplot 121; plot(abs(linePRE.'))  ; hold on; 
                %subplot 122; plot(abs(linePOST.'))  ; hold on; 
                %title(sprintf('Coil %d',s))
                
                Ypost = rec.(t2r);
                Yidx = reshape(idx,size(rec.(t2r)));
                plotND([], abs(( cat(5, Yidx, Ypre, Ypost))),[0 1e-4],[],0,{[],2},[],[],[],s);
                title(sprintf('Coil %d',s))
                
           end
                 
        for l=1
            rec.(t2r)= fftshiftGPU(rec.(t2r),l); 
            rec.(t2r)= fftGPU(rec.(t2r),l);     %*NY(l)  YB: this has once given inf as number close to realmax('single') 
            rec.(t2r)= ifftshiftGPU(rec.(t2r),l);
        end

        rec.(t2r)=resampling(rec.(t2r),NImageCols,2);
        if ~isempty(vr);rec.(t2r)=dynInd(rec.(t2r),vr,1);end %Extract supported FOV
    end
    typ2Rec={'y'}; % YB: all the channels are aready in the rec.N so loop only needs to permute + SupportFOV once
    y=cat(4,y,gather(rec.y));
end
rec.y=y;y=[];if gpu;rec.y=gpuArray(rec.y);end

if isfield(rec,'N');rec.y=standardizeCoils(rec.y,rec.N);end
fprintf('Number of inf in data: %d\n',multDimSum(isinf(abs(rec.y))));

NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    rec.y=fftshiftGPU(rec.y,n);
    rec.y=fftGPU(rec.y,n);%*NY(l);
    rec.y=ifftshiftGPU(rec.y,n);
end

% YB: now all dimensions are in the image domain (RO-PE-SL)
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.AcqVoxelSize=[TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV/NImageCols TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV/TW.hdr.Meas.NImageLins TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness/TW.hdr.Meas.NImagePar];

%CONVERT TO PRS (PE-RO-SL)
rec.y=gather(rec.y);rec.N=gather(rec.N);
rec.y = permute(rec.y,[2 1 3 4]);
rec.Enc.AcqVoxelSize = permute(rec.Enc.AcqVoxelSize,[2 1 3 4]);%need to permute spacing as well

%% GEOMETRY COMPUTATION
% asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1}; %See invert7T_old to see how to use this
% sNormal=asSlice.sNormal; sPosition=asSlice.sPosition;
slicePos = TW.image.slicePos;

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize' 1]);

%ROTATION 
quaternionRaw = slicePos(4:7,1);
translationRaw =  slicePos(1:3,1); %in mm for center FOV (I think)

rec.Par.Mine.Arot = eye(4);
rec.Par.Mine.Arot(1:3,1:3) = quaternion_rotation_matrix(quaternionRaw); %PE-RO-SL to PCS

%TRANSLATION
rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

% account for fact that tranRaw is referred to the center of FOV, not the first element in the array
orig = ceil((size(rec.y)+1)/2 - 1)'; orig=orig(1:3);
if ~isempty(vr); orig(1)=orig(1)+vr(1); end

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)...
                        - rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;
                    
%COMBINED MATRIX
rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

%MOVE TO RAS
rec.Par.Mine.PCS2RAS = diag([-1 -1 1 1]);%For HeadFirstSupine. Can be generalised from TW.hdr.Dicom.tPatientPosition  
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
%At this point rec.Par.Mine.APhiRec defined from PRS to RAS

rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;

%DEDUCE ACQUISITION ORDER
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane ] = acquisitionOrder(rec.Par.Mine.APhiRec);
fprintf('%s slices.\n', rec.Par.Scan.Mine.slicePlane);
rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir= rec.Par.Scan.MPS(5);%Positive RO

%% MAKE RO-PE-SL for DISORDER recon
perm = [2 1 3 4];
rec.y = permute(rec.y,perm);
rec.Par.Mine.APhiRec = rec.Par.Mine.APhiRec(:,perm);
rec.Enc.AcqVoxelSize = permute(rec.Enc.AcqVoxelSize,perm);%need to permute spacing as well

%%% Store the permutations performed from the PRS space
rec.Par.Mine.permuteHist = [];rec.Par.Mine.permuteHist{1} = perm;

%[rec.Par.Scan.MPS ] = acquisitionOrder(rec.Par.Mine.APhiRec);
%% SEQUENCE INFORMATION

rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be changed for shot-based sequences       
rec.Par.Labels.ZReconLength=1;

rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;

NY=size(rec.y);NY=NY(1:3);

offset = -1;
Lines = mod( TW.image.Lin -1 +offset, NY(2) ) +1;
Partitions = mod( TW.image.Par -1 +offset, NY(3) ) +1;

rec.Assign.z{2}= (NY(2) - (Lines-1))-floor((NY(2))/2) -1; % 2nd PE  %Minus comes from convention in invert7T.m
rec.Assign.z{3}= (NY(3) - (Partitions-1))-floor((NY(3))/2)-1; % 3rd PE = slices
% rec.Assign.z{2}=TW.image.Lin-floor(NY(2)/2)-1; %YB: 2nd PE
% rec.Assign.z{3}=TW.image.Par-floor(NY(3)/2)-1; %YB: 3rd PE - slices


%%% FOR DEBUGGING ONLY
% xRaw = sqrt( mean(rec.y.^2, 4)); %PE-RO-SL since taken from TWIX data 
% 
% dimShow = 3; %Inspect dimension by plotting over it
% for i = 1:4:size(xRaw,dimShow)
%     id = [50 50 50]; id(dimShow)=i;
%     plot_(abs(xRaw),[],id);
%     pause(0.01)
% end
