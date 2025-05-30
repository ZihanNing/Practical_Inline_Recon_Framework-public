function rec=reconInvert_inline(rec,typ,field)

%RECONINVERT   Inverts data
%   REC=RECONINVERT(REC,TYP,{FIELD})
%   * REC is a recon structure
%   * TYP is the type of data to invert, one of the following, 1->body, 
%2->surface, 3->data
%   * {FIELD} is the field with header info, it defaults to 'image', other
%   options are 'refscan'
%   ** REC is a recon structure
%

if nargin<3 || isempty(field);field='image';end

typV=['B' 'S' 'x'];
t=typV(typ);
TW=rec.TW.(t);%This is the data to invert
if typ<=2;supportReadout=rec.Alg.supportReadoutS;else supportReadout=rec.Alg.supportReadout;end

ND=16;
gpu=useGPU;

%ACQUIRED GRID SIZES
rec.Enc.(t).AcqN = [size(rec.TWD.(t),1),size(rec.TWD.(t),3),size(rec.TWD.(t),4)];
if size(rec.TWD.(t),4)>1; rec.Enc.(t).ThreeD = 1; end
fprintf('Acquired size (no PF):%s\n',sprintf(' %d',rec.Enc.(t).AcqN));

%NUMBER OF COILS
rec.Enc.(t).CoilN=size(rec.TWD.(t),2);
fprintf('Number of coils: %d\n',rec.Enc.(t).CoilN);

%OVERSAMPLING FACTORS
rec.Enc.(t).Overs=ones(1,3);%Readout is treated indepedently in function invert readout
if ~isempty(TW.hdr.Meas.flPhaseOS);rec.Enc.(t).Overs(2)=1+TW.hdr.Meas.flPhaseOS;end
if ~isempty(TW.hdr.Meas.flSliceOS);rec.Enc.(t).Overs(3)=1+TW.hdr.Meas.flSliceOS;end
fprintf('Oversampling factors:%s\n',sprintf(' %.3f',rec.Enc.(t).Overs));

rec.Enc.(t).Accel=ones(1,3);%Readout is treated indepedently

%ACCELERATION FACTORS
rec.Enc.(t).RegridMode=TW.hdr.Meas.alRegridMode(1);
fprintf('Regrid mode: %d\n',rec.Enc.(t).RegridMode);
rec.Enc.(t).Accel(2)=1/TW.hdr.MeasYaps.sKSpace.dPhaseResolution;
rec.Enc.(t).Accel(3)=1/TW.hdr.MeasYaps.sKSpace.dSliceResolution;
fprintf('Acceleration factors:%s\n',sprintf(' %.3f',rec.Enc.(t).Accel));
if rec.Enc.(t).RegridMode<=1;rec.Enc.(t).AcqN(1)=rec.Enc.(t).AcqN(1)/TW.hdr.Dicom.flReadoutOSFactor;end

%PARTIAL FOURIER
if strcmp(field,'image')
%     rec.Enc.(t).AcqNNoPF=rec.Enc.(t).AcqN;
%     if rec.Enc.(t).ThreeD;rec.Enc.(t).AcqN=[rec.Enc.(t).AcqN(1) TW.hdr.Config.N0FPEFTLength TW.hdr.Config.N0F3DFTLength];
%     else rec.Enc.(t).AcqN=[rec.Enc.(t).AcqN(1) TW.hdr.Config.N0FPEFTLength rec.Enc.(t).AcqN(3)];
%     end
%     pad=rec.Enc.(t).AcqN-rec.Enc.(t).AcqNNoPF;
%     rec.Enc.(t).APF=cell(1,3);
%     for d=1:3
%         rec.Enc.(t).APF{d}=ones(rec.Enc.(t).AcqNNoPF(d),1,'single');
%         if useGPU;rec.Enc.(t).APF{d}=gpuArray(rec.Enc.(t).APF{d});end
%         rec.Enc.(t).APF{d}=padarray(rec.Enc.(t).APF{d},pad(d),0,'post');
%     end   
%     pad=[pad(1) 0 pad(2:3)];
%     rec.TWD.(t)=padarray(rec.TWD.(t),pad,0,'post');
%     if rec.Enc.(t).ThreeD;dims=[3 4];else dims=3;end
%     cshift=zeros(1,4);
%     for d=dims
%         if d==3;cshift(3)=ceil((rec.Enc.(t).AcqN(2)+1)/2)-TW.(field).centerLin(1);
%         else cshift(4)=ceil((rec.Enc.(t).AcqN(3)+1)/2)-TW.(field).centerPar(1);
%         end                
%         rec.Enc.(t).APF{d-1}=circshift(rec.Enc.(t).APF{d-1},cshift(d));
%     end
%     rec.TWD.(t)=circshift(rec.TWD.(t),cshift);
%     for d=1:3;rec.Enc.(t).APF{d}=shiftdim(rec.Enc.(t).APF{d},d-1);end

    % easier tool for zero-padding due to PF
    % ZN: taking care both partial Fourier and (reconstruction resolution < acquisition resolution)
    % ZN: only pad the Lin & Par, supposing that asymmetric echo is handled by AsymmetricEcho gadget
    matrix_size = rec.TW.(t).sampling_description.encoded_matrix; % to be padded matrix size
    Acq_matrix_size = rec.Enc.(t).AcqN; % current/acquired matrix size
    Lin_center = max(rec.TW.(t).image.centerLin);
    Par_center = max(rec.TW.(t).image.centerPar);
    Lin = rec.TW.(t).image.sort_Lin(:,1);
    Par = rec.TW.(t).image.sort_Par(:,1);
    % padding
    [rec.TWD.(t),kspace_Pad,Lin,Par] = zeroPadKSpace(rec.TWD.(t),matrix_size,Acq_matrix_size,Lin_center,Par_center,Lin,Par,1);
    % modify parameters accordingly
    rec.Enc.(t).AcqN(2:3) = size(rec.TWD.(t),3:4);
    A = rec.TWD.(t); A = A(:,:,:,:,1,1); A = permute(A, [1 3 4 2]); acc_A = sum(A, 4);           
    rec.Enc.(t).APF = findMatrixByNonZero(acc_A);
end
fprintf('Acquired size:%s\n',sprintf(' %d',rec.Enc.(t).AcqN));

%RECONSTRUCTED GRID SIZES
rec.Enc.(t).RecN=round(rec.Enc.(t).AcqN.*rec.Enc.(t).Accel);%Oversampling will be removed at the end
fprintf('Reconstructed size:%s\n',sprintf(' %d',rec.Enc.(t).RecN));

%OUTPUT GRID SIZES
rec.Enc.(t).OutN=round(rec.Enc.(t).AcqN.*rec.Enc.(t).Accel./rec.Enc.(t).Overs);%Oversamling will be removed at the end
fprintf('Output size:%s\n',sprintf(' %d',rec.Enc.(t).OutN));%This is not used

% aux=TW.hdr.MeasYaps.sSliceArray.asSlice{1};
% %SLICE SEPARATION
% rec.Enc.(t).SliceSeparation=0;
% if ~rec.Enc.(t).ThreeD
%     posSl1=[aux.sPosition.dSag aux.sPosition.dCor aux.sPosition.dTra];
%     aux=TW.hdr.MeasYaps.sSliceArray.asSlice{2};
%     posSl2=[aux.sPosition.dSag aux.sPosition.dCor aux.sPosition.dTra];
%     rec.Enc.(t).SliceSeparation=norm(posSl1-posSl2);
%     fprintf('Slice separation: %.2f\n',rec.Enc.(t).SliceSeparation);
% end

% aux=TW.hdr.MeasYaps.sSliceArray.asSlice{1};
%FULL ACQUIRED FOV
% if rec.Enc.(t).ThreeD
%     rec.Enc.(t).AcqFOV=[aux.dReadoutFOV aux.dPhaseFOV aux.dThickness];rec.Enc.(t).SliceThickness=0;
% else 
%     rec.Enc.(t).AcqFOV=[aux.dReadoutFOV aux.dPhaseFOV rec.Enc.(t).SliceSeparation*rec.Enc.(t).AcqN(3)];rec.Enc.(t).SliceThickness=aux.dThickness;
%     fprintf('Slice thickness: %.2f\n',rec.Enc.(t).SliceSeparation);
% end
rec.Enc.(t).AcqFOV = TW.sampling_description.recon_fov; % used the FOV and matrix size provided by the header
rec.Enc.(t).AcqFOV=rec.Enc.(t).AcqFOV.*rec.Enc.(t).Overs;
fprintf('Acquired field of view:%s\n',sprintf(' %.2f',rec.Enc.(t).AcqFOV));

%ACQUIRED VOXEL SIZE
rec.Enc.(t).AcqDelta=rec.Enc.(t).AcqFOV./rec.Enc.(t).AcqN(1:3);
fprintf('Sampled resolution:%s\n',sprintf(' %.2f',rec.Enc.(t).AcqDelta));

%TO EXTRACT FOV
if ~isempty(supportReadout);vr=round(rec.Enc.(t).AcqN(1)*supportReadout(1)):round(rec.Enc.(t).AcqN(1)*supportReadout(2));else vr=1:rec.Enc.(t).AcqN(1);end
%if ~isempty(supportReadout);vr=round(rec.Enc.(t).RecN(1)*supportReadout(1)):round(rec.Enc.(t).RecN(1)*supportReadout(2));else vr=1:rec.Enc.(t).RecN(1);end

%SLICE AND PHASE ORDER
NX=size(rec.TWD.(t));
iSlicePositionsMin=1;
if ~isempty(rec.Alg.maxNumberRepeats);rec.TWD.(t)=dynInd(rec.TWD.(t),1:min(size(rec.TWD.(t),9),rec.Alg.maxNumberRepeats),9);end

%%%%%%%%%%%%%%%%%%%%%%%GEOMETRY%%%%%%%%%%%%%%%%%%%%%%%
% N=rec.Enc.(t).AcqN;
% MS=rec.Enc.(t).AcqDelta;
% S=diag([MS 1]);
% rec.Geom.(t).MS=MS;
% slicePos=TW.(field).slicePos(:,iSlicePositionsMin);
% quaternionRaw = slicePos(4:7);
% R=eye(4);
% R(1:3,1:3)=convertNIIGeom(ones([1,3]),quaternionRaw','qForm','sForm');%PE-RO-SL to PCS
% R=R(:,[2 1 3 4]);
% T=eye(4);T(1:3,4)=slicePos(1:3);%in mm for center FOV (I think)
% %MT=T*R(:,[2 1 3 4])*S;%Of the center of the FOV
% MT=T*R*S;%Of the center of the FOV
% Nsub=N/2;
% if ~rec.Enc.(t).ThreeD;Nsub(3)=1/2;end
% %Nsub(1)=Nsub(1)-2;%---it is the solution
% T(:,4)=MT*[-floor(Nsub)';1];
% %T(:,4)=MT*[-Nsub';1];
% %MT=T*R(:,[2 1 3 4])*S;%Of the origin
% MT=T*R*S;
% rec.Geom.(t).patientPosition=TW.hdr.Dicom.tPatientPosition;
% fprintf('Patient position: %s\n',rec.Geom.(t).patientPosition);
% rec.Geom.(t).PCS2RAS=getPCS2RAS(rec.Geom.(t).patientPosition);
% rec.Geom.(t).MT=rec.Geom.(t).PCS2RAS*MT;
% % tI=rec.Geom.(t).MT*[vr(1)-1;(rec.Enc.(t).OutN(2:3)-rec.Enc.(t).RecN(2:3))';1]; % ZN: debug, should be here, but comment out
% % %[vr(1)-1;(rec.Enc.(t).OutN(2:3)-rec.Enc.(t).RecN(2:3))']
% % %tI
% % %tI=rec.Geom.(t).MT*[vr(1)-1;0;0;1];
% % rec.Geom.(t).MT(1:3,4)=tI(1:3);

% for the twix-like structure, some geom related headers are still missing
% we use another way to compute the geom
if typ == 2 % refscan
    % In the new geom computation method, the geom of ACS is compute based
    % on the geom of the high-res image
    
    % matrix info about high-res image
    matrixSize_x = rec.TW.S.sampling_description.encoded_matrix; % Assume the NImageCols for refscan & high-res is the same
    AcqVoxelSize_x=rec.Enc.S.AcqFOV./matrixSize_x; % the AcqFOV of refscan & high-res should be the same
    ACSSize = rec.Enc.S.AcqN;
    TablePosTra = TW.hdr.Dicom.lGlobalTablePosTra;
    
    % include some hdr information
    image_hdr = rec.TW.S.image; % small hdrs of each kspace lines
    hdr = rec.TW.S.hdr;
    
    % geom computation
    geom = computeGeom(AcqVoxelSize_x,matrixSize_x,hdr,image_hdr,TablePosTra,1,ACSSize);
    
    % rec structure
    rec.Geom.S.PCS2RAS = geom.PCS2RAS; % ACS
    rec.Geom.S.MT = geom.APhiACS;
    rec.Geom.S.MS=rec.Enc.S.AcqDelta;
    rec.Geom.x.PCS2RAS = rec.Geom.S.PCS2RAS; % image
    rec.Geom.x.MT = geom.APhiRec;
    rec.Geom.x.MS=AcqVoxelSize_x;
    
end

%%%%%%%%%%%%%%%%%%%%%%%INVERSION%%%%%%%%%%%%%%%%%%%%%%%
%INVERT NOISE
% if isfield(TW,'noise')
%     %READ NOISE
%     if rec.Alg.removeOversamplingReadout;TW.noise.flagRemoveOS = true;end
%     rec.N.(t)=TW.noise.unsorted;
%     if gpu;rec.N.(t)=gpuArray(rec.N.(t));end    
%     %INVERT READOUT
%     rec.N.(t)=gather(invertReadout(rec.N.(t),rec.Enc.(t).GridMatrix));
% else
%     rec.N.(t)=[];
% end

%INVERT PHASE
if strcmp(field,'image');fieldPhase='phasecor';
elseif strcmp(field,'refscan');fieldPhase='refscanPC';
else error('Unrecognized field: %s',field);
end
% if isfield(TW,fieldPhase) && ~isfield(rec,'PS')
%     rec.P.(t)=dynInd(TW.(fieldPhase),{':'},ND);   
%     %rec.P.(t)=dynInd(TW.refscanPC,{':'},ND);
%     rec.P.(t)=sum(rec.P.(t),3);%Sum over the phase encodes, assumes that typically one, but sometimes with zeros 
%     rec.P.(t)=sum(rec.P.(t),6)./(sum(single(rec.P.(t)~=0),6)+eps);%Sum over the averages (changed recently to include the eps, necessary in case fMRI_2023_07_04)
%     rec.P.(t)=mean(rec.P.(t),[5 9]);%Average over the repeats       
%     %rec.P.(t)=prod(rec.P.(t),2).^(1/size(rec.P.(t),2));   
%     rec.P.(t)=mean(rec.P.(t),2); 
%     if gpu;rec.P.(t)=gpuArray(rec.P.(t));end
%     rec.P.(t)=invertReadout(rec.P.(t),rec.Enc.(t).GridMatrix);
% %    visGhost(dynInd(rec.P.(t),1,9))
%     rec.P.(t)=gather(ridgeDetection(rec.P.(t),1,32));  
% %    visGhost(dynInd(rec.P.(t),1,9))
% %1
% else
%     rec.P.(t)=[];
% end
% rec.P.echoOrder=TW.(field).Seg(TW.(field).Sli==1 & TW.(field).Rep==1);
% rec.P.corrLine=TW.(field).Lin(TW.(field).Sli==1 & TW.(field).Rep==1);

%INVERT READOUT DATA
rec.(t)=[];
blkSz=8;%Depends on available GPU memory
%rec.TWD.(t)=dynInd(rec.TWD.(t),1:min(size(rec.TWD.(t),9),1),9);%Only one dynamic
rec.Enc.(t).GridMatrix=[];
for s=1:blkSz:rec.Enc.(t).CoilN;vS=s:min(s+blkSz-1,rec.Enc.(t).CoilN);%Inversion by channel groups
    %EXTRACT
    x=dynInd(rec.TWD.(t),vS,2);
    if gpu;x=gpuArray(x);end        
    %INVERT READOUT
    x=invertReadout(x,rec.Enc.(t).GridMatrix); 
    %ASSIGN
    rec.(t)=cat(4,rec.(t),gather(x));
end
rec.TWD.(t)=[];%Releasing memory

%GENERATE MASKS FOR KY AND MB
rec.Ay=sum(abs(rec.(t)).^2,setdiff(1:ND,2));
rec.Ay=single(rec.Ay>1e-12);
if strcmp(field,'image');rec.Ay=rec.Ay+(1-rec.Enc.(t).APF{2});end

% if (~isempty(rec.P.(t)) || isfield(rec,'PS')) && strcmp(field,'image')
%     for n=1:max(rec.P.echoOrder)
%         in=rec.P.corrLine(rec.P.echoOrder==n);   
%         in=unique(in);%Recent addition to prevent situations where indexes are repeated, a bit risky
%         rec.Ay(in)=rec.Ay(in)*n;
%     end
% end

%MULTIBAND
rec.Enc.(t).MultiBandFactor=1;
rec.Ay=ifftshiftGPU(rec.Ay,2);

%INVERT
blkSz=blkSz*2;
for s=1:blkSz:rec.Enc.(t).CoilN;vS=s:min(s+blkSz-1,rec.Enc.(t).CoilN);%Inversion by channel groups
    %EXTRACT
    x=dynInd(rec.(t),vS,4);
    if gpu;x=gpuArray(x);end        
    %INVERT
    for n=2:3;x=ifftInvert(x,n);end    
    %ASSIGN 
    rec.(t)=dynInd(rec.(t),vS,4,gather(x));
end

%REORDER, MAP SLICES TO THIRD DIMENSION, AND REMOVE UNUSED SLICES FOR MB
perm=1:ND;perm([3 5])=[5 3];


%%%%%%%%%%%%%%%%%%%%%%%AUXILIARY FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%
function x=invertReadout(x,gridMatrix)
    NF=size(x,1);%Input size    

    %PERMUTE COILS TO 4TH DIMENSION
    perm=1:ND;perm(2:4)=[3:4 2];
    x=permute(x,perm);

    %REGRIDDING MATRIX
    if ~isempty(gridMatrix);A=gridMatrix;else A=eye(NF,'like',x);end

    %INVERSION MATRIX
    F=build1DFTM(NF,0,isa(x,'gpuArray'));
    F=ifftshift(fftshift(F,2),1);

    %REMOVE OVERSAMPLING AND CROP TO DEFINED AREA
    F=dynInd(resampling(F,rec.Enc.(t).AcqN(1),3),vr,1);
    
    %APPLY
    x=aplGPU(F*A,x,1);
end

function x=ifftInvert(x,m)
    NF=size(x,m);%Input size
    %m

    %INVERSION MATRIX
    F=build1DFTM(NF,0,isa(x,'gpuArray'));
    F=ifftshift(fftshift(F,2),1);
    %F=circshift(F,ceil((rec.Enc.(t).OutN(m)-rec.Enc.(t).RecN(m))/2),1);
    %APPLY
    x=aplGPU(F,x,m);
    %m
    %any(isnan(x(:)))
end

end
