pathCode='/home/lcg13/Work/DefinitiveImplementationDebug05';
addpath(genpath(pathCode));
pathDe='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';

%path='2017_11_24/HO_20211';%Full run fMRI MB3S1.4+DWI
%file='Dy-Di/ConsistencyTests/ho_24112017_1039418_24_2_zs1dmbst142sensep1o0V4';

%path='2017_11_22/FE_19710';%Full run fMRI MB3S1.4+DWI
%file='Dy-Di/ConsistencyTests/fe_22112017_1515251_29_2_zs1dmbst142sensep1o0V4';

%path='2017_12_14/ba_25112';
%file='Dy-Di/ba_14122017_1807318_3_2_zs3dmbst142sensep1o0V4';

%path='2018_01_03/AD_27310';%Full run fMRI MB3S1.4+DWI-FD on
%file='Dy-Di/ad_03012018_1107544_23_2_zs3dmbst142sensep1o0V4';

%direFile='dhcp_fetal142_c.txt';

path='2018_01_05/WH_27710';%Full run fMRI MB3S1.4+DWI-FD on
file='Dy-Di/wh_05012018_1030587_24_2_zs3dmbst140fdfidsensep1o0V4';

direFile='dhcp_fetal140.txt';

gpu=1;

suff{1}='_Vo';suff{2}='_Ma';suff{3}='_Tr';
for n=1:length(suff)
    %suff{n}=strcat(suff{n},'ESD');    
    %suff{n}=strcat(suff{n},'NonUndistorted');
    %suff{n}=strcat(suff{n},'NonDenoisedNonUndistorted');
end

nii=load_untouch_nii(sprintf('%s/%s/%s%s.nii',pathDe,path,file,suff{2}));
MS=nii.hdr.dime.pixdim(2:4);
MT=eye(4);MT(1,:)=nii.hdr.hist.srow_x;MT(2,:)=nii.hdr.hist.srow_y;MT(3,:)=nii.hdr.hist.srow_z;
M=nii.img;
M(M<1e-6)=0;M(M>=1e-6)=1;

load(sprintf('%s/%s/%s%s.mat',pathDe,path,file,suff{3}));%T

phFile=sprintf('%s/Data/DynamicFiles/%s',pathCode,direFile);
phData=load(phFile);
bval=phData(:,1:4);
bun=unique(bval(:,4));

nii=load_untouch_nii(sprintf('%s/%s/%s%s.nii',pathDe,path,file,suff{1}));
%To save the spin echo
%niftiIm=make_nii(nii.img(:,:,:,1:2:end),MS);
%niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1;         
%save_nii(niftiIm,sprintf(sprintf('%s/%s/%s%sSPINECHO.nii',pathDe,path,file,suff{1})));
%return

x=nii.img;
NX=size(x);
x=dynInd(x,1:2:NX(4),4);
NX=size(x);
x=reshape(x,[prod(NX(1:3)) NX(4)]);
med=median(x(logical(M),:),1);
%x=bsxfun(@times,x,1./med);
x=reshape(x,NX);
M=morphFourier(M,15,[1 1 1]);
%M(:)=0;
%M(65,57,27:29)=1;

ROI=computeROI(M);
x=extractROI(x,ROI,1);
M=extractROI(M,ROI,1);

NX=size(x);
T=shiftdim(T,1);
T=repmat(T,[1 1 NX(3) 1 1]);

%TRANSFORM GRADIENTS
if gpu;T=gpuArray(T);x=gpuArray(x);M=gpuArray(M);bval=gpuArray(bval);end

%Here we may need the following coordinate transform (now it was all -1)
%MR.Transform('rec', 'xyz')

bval=permute(bval,[2 3 4 1]);
bval=repmat(bval,[1 1 NX(3) 1]);
%bval(3,:,:,:)=-bval(3,:,:,:);
%T(:)=0;
bval=dynInd(bval,1:3,1,transformGradients(dynInd(bval,1:3,1),T));
bval=permute(bval,[5 2 3 4 1]);

%Here not sure whether the ODF information should be "multiplied" by MT
order=6;
no=[2 2];%Fitting norm for b0 and other b-vals
%MT
%pause
%%%BE CAREFUL WITH CONVERSION OF COORDINATES INSIDE TO SPHERICAL HARMONICS!!!
x=solveODF(x,bval,order,M,no);


% for s=2:4:size(x,3)
% figure
% TensorODF=dynInd(x,{s,2:7},3:4);
% ImageODF=dynInd(x,[s 1],3:4);
% TensorODF=permute(TensorODF,[4 1 2 3]);
% subplot(1,2,1);
% plotTensors(gather(TensorODF),1,[81 1],gather(ImageODF));
% 
% %subplot(1,2,2);
% %TensorODF=dynInd(x,{s,8:13},3:4);
% %ImageODF=dynInd(x,[s 1],3:4);
% %TensorODF=permute(TensorODF,[4 1 2 3]);
% %plotTensors(gather(TensorODF),1,[81 1],gather(ImageODF));
% end

for m=1:length(x)
    x{m}=extractROI(x{m},ROI,0);
    niftiIm=make_nii(gather(x{m}),MS);
    niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1;         
    save_nii(niftiIm,sprintf(sprintf('%s/%s/%s%sB%.2f-N%d-P%d-%d.nii',pathDe,path,file,suff{1},bun(m)/1000,order,no(1),no(2))));
end
