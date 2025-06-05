pathCode='/home/lcg13/Work/DefinitiveImplementationDebug05';
addpath(genpath(pathCode));
pathDe='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';

%path='2017_11_24/HO_20211';%Full run fMRI MB3S1.4+DWI
%file='Dy-Di/ho_24112017_1039418_24_2_zs1dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2017_11_22/FE_19710';%Full run fMRI MB3S1.4+DWI
%file='Dy-Di/fe_22112017_1515251_29_2_zs1dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2017_12_14/ba_25112';
%file='Dy-Di/ba_14122017_1807318_3_2_zs3dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2018_01_03/AD_27310';%Full run fMRI MB3S1.4+DWI-FD on
%file='Dy-Di/ad_03012018_1107544_23_2_zs3dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2018_01_05/WH_27710';%Full run fMRI MB3S1.4+DWI-FD on
%file='Dy-Di/wh_05012018_1030587_24_2_zs3dmbst140fdfidsensep1o0V4';
%direFile='dhcp_fetal140.txt';

%path='2018_01_10/BE_28910';%Full run fMRI MB3S1.4+DWI-FD on
%file='Dy-Di/be_10012018_1224282_23_2_zs3dmbst142fdfsensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2018_01_12/DO_29810';
%file='Dy-Di/do_12012018_0957487_11_2_sf4dmbst90tr9sensep1o0V4';
%direFile='dhcp_fetal90.txt';

%path='2018_01_12/DO_29810';
%file='Dy-Di/do_12012018_1051031_26_2_sf4dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2018_01_15/AT_30010';
%file='Dy-Di/at_15012018_1003231_11_2_sf4dmbst90tr9sensep1o0V4';
%direFile='dhcp_fetal90.txt';

%path='2018_01_15/AT_30010';
%file='Dy-Di/at_15012018_1047581_22_2_sf4dmbst142sensep1o0V4';
%direFile='dhcp_fetal142_c.txt';

%path='2018_01_19/LE_30910';
%file='Dy-Di/le_19012018_1003592_13_2_sf4dmbst90tr9sensep1o0V4';
%direFile='dhcp_fetal90.txt';

%path='2018_01_19/LE_30910';
%file='Dy-Di/le_19012018_1054366_26_2_sf4dmbmtrst140fdsensep1o0V4';
%direFile='dhcp_fetal140.txt';

%path='2018_01_24/DE_31910';
%file='Dy-Di/de_24012018_0959172_12_2_sf4dmbst90tr9sensep1o0V4';
%direFile='dhcp_fetal90.txt';

%path='2018_01_24/DE_31910';
%file='Dy-Di/de_24012018_1054074_27_2_sf4dmbmtrst140fdsensep1o0V4';
%direFile='dhcp_fetal140.txt';

path='2018_01_29/EN_33310';
file='Dy-Di/en_29012018_1051182_30_2_sf7dmbfinalsensep1o0V4';
direFile='dhcp-fetal-final.txt';

suff{1}='_Un';%ESD';
suff{2}='_Vo';%ESD';
%suff{3}='_ExFrac1.0';
suffM='_Ma';%ESD';
nii=load_untouch_nii(sprintf('%s/%s/%s%s.nii',pathDe,path,file,suffM));
MS=nii.hdr.dime.pixdim(2:4);
MT=eye(4);MT(1,:)=nii.hdr.hist.srow_x;MT(2,:)=nii.hdr.hist.srow_y;MT(3,:)=nii.hdr.hist.srow_z;
M=nii.img;M(M<1e-6)=0;M(M>=1e-6)=1;

phFile=sprintf('%s/Data/DynamicFiles/%s',pathCode,direFile);
phData=load(phFile);
bval=phData(:,4);
nb=unique(bval);

x=cell(1,length(suff));
for d=1:length(suff)%Distorted/Undistorted       
    nii=load_untouch_nii(sprintf('%s/%s/%s%s.nii',pathDe,path,file,suff{d}));
    x{d}=nii.img;
    N=size(x{d});
    x{d}=reshape(x{d},[prod(N(1:3)) N(4)]);    
    xM=bsxfun(@times,x{d},1./median(x{d}(M==1,:),1));
    xM=reshape(xM,N);
    x{d}=reshape(x{d},N);
    niftiIm=make_nii(gather(xM),MS);
    niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1; 
    save_nii(niftiIm,sprintf(sprintf('%s/%s/%s%sNorm.nii',pathDe,path,file,suff{d})));        
    NV=size(x{d},4);    
    for l=1:2
        xE=dynInd(nii.img,l:2:NV,4);    
        for n=1:length(nb)
            xB=dynInd(xE,bval==nb(n),4);
            xB=mean(xB,4);
            niftiIm=make_nii(gather(xB),MS);
            niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1;         
            save_nii(niftiIm,sprintf(sprintf('%s/%s/%s%sB%04d-%d.nii',pathDe,path,file,suff{d},nb(n),l))); 
        end   
    end
end





% x=cell(1,2);
% for n=1:length(file)
%     nii=load_untouch_nii(sprintf('%s%s%s',path,file{n},'_B0.nii'));
%     x{n}=nii.img;
%     %if n==2;x{n}=flip(flip(x{n},2),1);end
% end
% x=cat(4,x{:});if rec.Dyn.GPU;x=gpuArray(x);end
% rec.M=sum(abs(x),4);
% rec.M(rec.M>1e-1)=1;rec.M(rec.M<1e-1)=0;  
% 
% %ROI EXTRACTION
% rec.Enc.ROI=computeROI(rec.M);
% 
% fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
% %Types={'M','x','u'};
% Types={'M','u'};
% for n=1:length(Types);datTyp=Types{n};
%     rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1);
% end
% 
% rec=volumeAlignment(rec);
% %if ~rec.Alg.parU.corrMotion
% %    if ~rec.Alg.parU.useUndist;rec.v=rec.x;rec.e=rec.x;else rec.v=rec.u;rec.e=rec.u;end
% %elseif rec.Alg.parU.corrMotion
% %    if ~rec.Alg.parU.useUndist;rec.e=rec.x;else rec.e=rec.u;end
% %end
% 
% suffType={'YeCoVoNoMed.nii','YeCoExNoMed.nii','NoCoNoMed.nii'};
% 
% suffOu={'_Ma','_Vo','_Ex','_Di'};
% Types={'M','v','e','d'};
% for n=1:length(Types);datTyp=Types{n};
%     rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0);
%     niftiIm=make_nii(gather(abs(rec.(datTyp))),MS);
%     niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1; 
%     if rec.Alg.parU.corrMotion==0;sufN=3;else sufN=rec.Alg.parU.corrMotion;end
%     save_nii(niftiIm,sprintf('%s%s%s%s',path,file{1},suffOu{n},suffType{sufN}));
% end
% if rec.Dyn.Typ2Wri(17);T=gather(rec.T);save(sprintf('%s%s%s%s',path,file{1},sprintf('_Tr%d',rec.Alg.parU.corrMotion,suffType{rec.Alg.parU.corrMotion})));end
% return
% 
% 
% pathIn='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_10_19/OM_10210/Dy-Di';
% pathOu='/home/lcg13/Articulos/02EnProceso/18_ISMRM_SAFE/Figs/Fig4';
% 
% file{1}{1}='DiB0000-80';
% file{1}{2}='DiB0400-50';
% file{1}{3}='DiB0700-42';
% file{1}{4}='DiB1000-38';
% file{2}{1}='UnB0000-80';
% file{2}{2}='UnB0400-50';
% file{2}{3}='UnB0700-42';
% file{2}{4}='UnB1000-38';
% file{3}{1}='FiB0000';
% file{3}{2}='FiB0400';
% file{3}{3}='FiB0700';
% file{3}{4}='FiB1000';
% 
% textA={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
% 
% l2=[350 350];
% 
% C=[];
% cont=1;
% xau=20;
% yau=40;
% for n=1:length(file)
%     B=[];
%     for m=1:length(file{1})        
%         A=imread(sprintf('%s/%s.png',pathIn,file{n}{m}));        
%         A=permute(A(:,1+l2(1):end-l2(2),:),[2 1 3]);
%         figure
%         imshow(A,[])
%         hold on
%         text(xau,yau,textA{cont},'Color','white','Fontsize',40)
%         %set(gca,'color','w')        
%         set(gcf,'color','k');        
%         set(gcf, 'InvertHardCopy', 'off');
% 
%         pause(1)
%         print(gcf,'-dtiffnocompression',sprintf('%s/Fig4-%02d.tiff',pathOu,cont));
%         print(gcf,'-dpng',sprintf('%s/Fig4-%02d.png',pathOu,cont))
%         close all
%         cont=cont+1;
%     end
% end
