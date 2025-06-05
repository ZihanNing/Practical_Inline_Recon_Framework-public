function rec=invert7T_PilotTone_old(TW,suppFOV)

%INVERT7T   Inverts SIEMENS data
%   * TW is the TWIX object
%   * SUPPFOV is the support in the readout direction for high resolution
%   scans
%   ** REC is a reconstruction structure
%

if nargin<2;suppFOV=[];end

ND=16; % YB: Number of dimensions that can be called in the TWIX object (see public method subsref)
TW.image.flagIgnoreSeg=1;

gpu=(gpuDeviceCount>0 && ~blockGPU);

typ2Rec={'y'};
if isfield(TW,'noise')
    typ2Rec=[typ2Rec 'N'];
    %rec.N=dynInd(TW.noise,':',ND);
    rec.N=TW.noise.unsorted;
    if gpu;rec.N=gpuArray(rec.N);end
end

NCha=TW.image.NCha;
NImageCols=TW.hdr.Meas.NImageCols;

y=[];
if ~isempty(suppFOV);vr=round(NImageCols*suppFOV(1)):round(NImageCols*suppFOV(2));else vr=[];end

for s=1:NCha
    rec.y=dynInd(TW.image,{s ':'},[2 ND]);
    if gpu;rec.y=gpuArray(rec.y);end % YB: Here, 1st = Readout (HF), 2nd =  channel, 3rd = LR (fastest PE direction) and 4rd = AP (slowest readout)

    perm=1:ND;
    perm(2:4)=[3:4 2]; % YB:Re-arrange so that 1st = Readout (HF), 2nd = LR (fastest PE direction), 3rd = AP (slowest readout) and 4th = Channel
    for n=1:length(typ2Rec); t2r=typ2Rec{n};
        rec.(t2r)=permute(rec.(t2r),perm);
        NY=size(rec.(t2r));NY(end+1:ND)=1;


           fprintf('Coil %d\n',s);
           th1=1e+37;
           if  multDimSum(isinf(th1*rec.(t2r)))>0
               fprintf('high values detected in %s\n',t2r);
                 rec.(t2r)(isinf(abs(rec.(t2r)))) = 1e+30;%just high
                           %for dd=1:1:size(rec.y,2)
                 id = ceil((size(rec.y)+1)/2);
                 id(2) = 51;
                 %plot_(isinf(abs(rec.y)),[],id)
                 %plot_((abs(rec.y)),[0 0.001],id, isinf(abs(rec.y)))
                 %plot_((abs(rec.y)),[0 0.001],id, isnan(abs(rec.y)))
                 %plot_((abs(meanPE_interpolation)),[0 0.001],id, isinf(abs(meanPE_interpolation)))
                 
                 pause(0.1)
            %end 
            
            
            
            %peak identification
            h = ones([1,8,8])/(11^2-1);h(:,4,4)=0;%81-1 since center value = 0
            th = 9;
            meanPE = imfilter(abs(rec.(t2r)),h); 
            idx = (abs(rec.(t2r)) > th * meanPE);
%             plot_((abs(meanPE)),[0 0.001],id, idx)
%             plot_((abs(meanPE)),[0 0.001],id, idx)
            %interpolation kernel
            %%%%h_interpolation = ones([3, 1,1])/(3-1);h_interpolation(2,1,1)=0
            h_interpolation = ones([1,4,4])/(9-1);h_interpolation(:,2,2)=0;
            %interpolate neighbourhing samples
            linePRE = dynInd(rec.(t2r), [51, 81],2:3);
            meanPE_interpolation = imfilter(rec.(t2r),h_interpolation);
       %plot_((abs(meanPE_interpolation)),[0 0.001],id, isinf(abs(meanPE_interpolation)))
       %plot_((abs(idx)),[0 0.001],id, isinf(abs(meanPE_interpolation)))
            rec.(t2r) (idx ) =  meanPE_interpolation ( idx ) ;
            linePOST  = dynInd(rec.(t2r), [51, 81],2:3);
            
            figure;
            subplot 121; plot(abs(linePRE.'))  ; hold on; 
            subplot 122; plot(abs(linePOST.'))  ; hold on; 
            %plot(1:length(line),dynInd(idx.*abs(line), [51, 81],2:3))
            title(sprintf('Coil %d',s))
            pause(0.01)
           end
                 
%             %for dd=1:1:size(rec.y,2)
%                  id = ceil((size(rec.y)+1)/2);
%                  id(2) = 51;
%                  %plot_(isinf(abs(rec.y)),[],id)
%                  plot_((abs(rec.y)),[0 0.001],id, isinf(abs(rec.y)))
%                  pause(0.1)
%             %end 
            
multDimSum(isinf(abs(rec.y)))
        for l=1
            rec.(t2r)=fftshiftGPU(rec.(t2r),l); %Note that the natural order of fftshift is opposite to standard convention here, but should only be different with  
            rec.(t2r)=ifftGPU(rec.(t2r),l)*NY(l);%YB: this has once given inf as number close to realmax('single') 
            rec.(t2r)=ifftshiftGPU(rec.(t2r),l);
        end
multDimSum(isinf(abs(rec.y)))
        rec.(t2r)=resampling(rec.(t2r),NImageCols,2);
        if ~isempty(vr);rec.(t2r)=dynInd(rec.(t2r),vr,1);end %Extract supported FOV
    end
    typ2Rec={'y'}; % YB: all the channels are aready in the rec.N so loop only needs to permute + SupportFOV once
    y=cat(4,y,gather(rec.y));
end
rec.y=y;y=[];if gpu;rec.y=gpuArray(rec.y);end

multDimSum(isinf(abs(rec.y)))
if isfield(rec,'N');rec.y=standardizeCoils(rec.y,rec.N);end
multDimSum(isinf(abs(rec.y)))

NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    rec.y=fftshiftGPU(rec.y,n);
    rec.y=ifftGPU(rec.y,n)*NY(n);
    rec.y=ifftshiftGPU(rec.y,n);
end
% YB: now all dimensions are in the image domain
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.AcqVoxelSize=[TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV/NImageCols TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV/TW.hdr.Meas.NImageLins TW.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness/TW.hdr.Meas.NImagePar];

rec.y=gather(rec.y);rec.N=gather(rec.N);

%GEOMETRY COMPUTATION
asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1};
sNormal=asSlice.sNormal;
sPosition=asSlice.sPosition;

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

%ROTATION 
EV=zeros(1,4);
if isfield(sNormal,'dSag');EV(1)=sNormal.dSag;end
if isfield(sNormal,'dCor');EV(2)=sNormal.dCor;end
if isfield(sNormal,'dTra');EV(3)=sNormal.dTra;end % YB: These are the direction cosines
if isfield(asSlice,'dInPlaneRot'),EV(4)=asSlice.dInPlaneRot;end 

[~,max_cosd] = max(EV(1:3));
if max_cosd ==1; 
    fprintf('Sagital slices.\n'); sl = 's';
    rec.Par.Scan.MPS='HF PA LR';%YB: This assumes readout in HF direction
elseif max_cosd ==2; 
    fprintf('Coronal slices.\n'); sl = 'c';%YB why t??????????
    rec.Par.Scan.MPS='HF LR PA';%YB: This assumes readout in HF direction
else
    fprintf('Transversal slices.\n Attention: Transversal slices don''t have readout in HF! Check!.\n'); sl = 't';
    rec.Par.Scan.MPS='AP LR HF'; error('YB: Not sure if rec.Par.Scan.MPS is correct. Probably not DISORDER.')
end
rec.Par.Mine.Arot=vox2ras_rsolveAA(EV(1:3),EV(4),sl);rec.Par.Mine.Arot(4,4)=1; % YB: This finds the rotation matrix based on the direction cosines
rec.Par.Mine.Arot = rec.Par.Mine.Arot(:,[2 1 3 4]);

%TRANSLATION
rec.Par.Mine.Atra=eye(4);
if isfield(sPosition,'dSag');rec.Par.Mine.Atra(1,4)=sPosition.dSag;end
if isfield(sPosition,'dCor');rec.Par.Mine.Atra(2,4)=sPosition.dCor;end
if isfield(sPosition,'dTra');rec.Par.Mine.Atra(3,4)=sPosition.dTra;end
or=-([NImageCols TW.hdr.Meas.NImageLins TW.hdr.Meas.NImagePar]/2)';
if ~isempty(vr);or(1)=or(1)+vr(1);end

RAS2LPS =eye(4); 
%RAS2LPS(1,1) = -1; RAS2LPS(2,2) = -1; 
rec.Par.Mine.Arot = RAS2LPS * rec.Par.Mine.Arot; 

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)+rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*or;

rec.Par.Mine.MTT=[-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1];
rec.Par.Mine.MTT=inv(rec.Par.Mine.MTT);

%COMBINED MATRIX
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

%NOTE LR VS RL NEEDS FURTHER REVERSE ENGINEERING!
