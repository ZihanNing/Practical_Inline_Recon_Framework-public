addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug05'));
%path='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_11_03/CO_14310/Dy-Di/';
%file{1}='co_03112017_1003126_13_2_sf1dmbrsgrbvalapsensep1o0V4';
%file{2}='co_03112017_1007202_14_2_sf1dmbrsgrbvalpasensep1o0V4';

%path='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_11_06/JA_15410/Dy-Di/';
%file{1}='ja_06112017_1047010_24_2_sf1dmbrsgrsensep1o0V4';
%file{2}='ja_06112017_1051389_25_2_sf1dmbrsgrpasensep1o0V4';

path='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_11_24/HO_20211/Dy-Di/';
file{1}='ho_24112017_1039418_24_2_zs1dmbst142sensep1o0V4';

%path='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_11_22/FE_19710/Dy-Di/';%%%PROBLEM IN TEMPLATE FOR SVR
%file{1}='fe_22112017_1515251_29_2_zs1dmbst142sensep1o0V4';

%%path='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_11_17/MC_18710/Dy-Di/';
%%file{1}='mc_17112017_1111018_30_2_ztip3dmbrstbvalsensep1o0V4';

%suff{1}='_Aq.nii';
%suff{2}='_Un.nii';
suff{1}='_Un.nii';

rec.Dyn.GPU=1;
rec.Fail=0;
rec.Alg.parU.useUndist=1;
rec.Par.Mine.AdHocArray=[101 0 0 2];
rec.Par.Scan.AcqVoxelSize=2*ones(1,3);
rec.Alg.parU.corrMotion=2;
rec.Alg.parU.iterMask=1;
rec.Dyn.Typ2Rec=[];
rec.Dyn.Typ2Wri=zeros(1,20);
rec.Alg.parU.Lambda=1;

%indB0=[1 15 29 37 51];
indB0=[1 2 3 15 24 31 43 52 59 71 80 87 99 108 118 128 132];indB0=2*indB0-1;
for d=1:length(suff)%Distorted/Undistorted
    x=cell(1,2);
    for n=1:length(file)%AP-PA          
        nii=load_untouch_nii(sprintf('%s%s%s',path,file{n},suff{d}));
        if n==1 && d==1;
            MS=nii.hdr.dime.pixdim(2:4);
            MT=eye(4);MT(1,:)=nii.hdr.hist.srow_x;MT(2,:)=nii.hdr.hist.srow_y;MT(3,:)=nii.hdr.hist.srow_z;
        end
        x{n}=dynInd(nii.img,indB0,4);
        if n==2;x{n}=flip(flip(x{n},2),1);end
    end  
    %if d==1;
    %    rec.x=cat(4,x{:});
    %    if rec.Dyn.GPU;rec.x=gpuArray(rec.x);end
    %else
        rec.u=cat(4,x{:});
        if rec.Dyn.GPU;rec.u=gpuArray(rec.u);end
    %end   
end
x=cell(1,2);
for n=1:length(file)
    nii=load_untouch_nii(sprintf('%s%s%s',path,file{n},'_B0.nii'));
    x{n}=nii.img;
    %if n==2;x{n}=flip(flip(x{n},2),1);end
end
x=cat(4,x{:});if rec.Dyn.GPU;x=gpuArray(x);end
rec.M=sum(abs(x),4);
rec.M(rec.M>1e-1)=1;rec.M(rec.M<1e-1)=0;  

%ROI EXTRACTION
rec.Enc.ROI=computeROI(rec.M);

fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
%Types={'M','x','u'};
Types={'M','u'};
for n=1:length(Types);datTyp=Types{n};
    rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1);
end

rec=volumeAlignment(rec);
%if ~rec.Alg.parU.corrMotion
%    if ~rec.Alg.parU.useUndist;rec.v=rec.x;rec.e=rec.x;else rec.v=rec.u;rec.e=rec.u;end
%elseif rec.Alg.parU.corrMotion
%    if ~rec.Alg.parU.useUndist;rec.e=rec.x;else rec.e=rec.u;end
%end

suffType={'YeCoVoNoMed.nii','YeCoExNoMed.nii','NoCoNoMed.nii'};

suffOu={'_Ma','_Vo','_Ex','_Di'};
Types={'M','v','e','d'};
for n=1:length(Types);datTyp=Types{n};
    rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0);
    niftiIm=make_nii(gather(abs(rec.(datTyp))),MS);
    niftiIm.hdr.hist.srow_x=MT(1,:);niftiIm.hdr.hist.srow_y=MT(2,:);niftiIm.hdr.hist.srow_z=MT(3,:);niftiIm.hdr.hist.sform_code=1; 
    if rec.Alg.parU.corrMotion==0;sufN=3;else sufN=rec.Alg.parU.corrMotion;end
    save_nii(niftiIm,sprintf('%s%s%s%s',path,file{1},suffOu{n},suffType{sufN}));
end
if rec.Dyn.Typ2Wri(17);T=gather(rec.T);save(sprintf('%s%s%s%s',path,file{1},sprintf('_Tr%d',rec.Alg.parU.corrMotion,suffType{rec.Alg.parU.corrMotion})));end
return


pathIn='/home/lcg13/Data/rawDestin/ReconstructionsDebug05/2017_10_19/OM_10210/Dy-Di';
pathOu='/home/lcg13/Articulos/02EnProceso/18_ISMRM_SAFE/Figs/Fig4';

file{1}{1}='DiB0000-80';
file{1}{2}='DiB0400-50';
file{1}{3}='DiB0700-42';
file{1}{4}='DiB1000-38';
file{2}{1}='UnB0000-80';
file{2}{2}='UnB0400-50';
file{2}{3}='UnB0700-42';
file{2}{4}='UnB1000-38';
file{3}{1}='FiB0000';
file{3}{2}='FiB0400';
file{3}{3}='FiB0700';
file{3}{4}='FiB1000';

textA={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};

l2=[350 350];

C=[];
cont=1;
xau=20;
yau=40;
for n=1:length(file)
    B=[];
    for m=1:length(file{1})        
        A=imread(sprintf('%s/%s.png',pathIn,file{n}{m}));        
        A=permute(A(:,1+l2(1):end-l2(2),:),[2 1 3]);
        figure
        imshow(A,[])
        hold on
        text(xau,yau,textA{cont},'Color','white','Fontsize',40)
        %set(gca,'color','w')        
        set(gcf,'color','k');        
        set(gcf, 'InvertHardCopy', 'off');

        pause(1)
        print(gcf,'-dtiffnocompression',sprintf('%s/Fig4-%02d.tiff',pathOu,cont));
        print(gcf,'-dpng',sprintf('%s/Fig4-%02d.png',pathOu,cont))
        close all
        cont=cont+1;
    end
end
