function svr=svrSetUp(svr)

%SVRSETUP   Sets up the SVR information
%   SVR=SVRSETUP(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   ** SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

gpu=svr.rec{1}.Dyn.GPU;
svr.NV=length(svr.y);
if gpu%GPU arrays
    svr.y=cellfun(@abs,cellfun(@gpuArray,svr.y,'UniformOutput',false),'UniformOutput',false);
    svr.El=cellfun(@gpuArray,svr.El,'UniformOutput',false);
end

%To perform quick checks on geometry
%svr.y=[];svr.rec=[];svr.El=[];
%save('/home/lcg13/Work/DataDefinitiveImplementationDebug07/svr.mat','svr','-v7.3');

%WE USE THE ALIGNED CENTERS TO DECIDE ON FOVS---fOV OF IMAGING CUBOIDS, 
%LIMITS OF CUBOIDS, INDEXES OF LIMITS OF CUBOIDS, MATERIAL FOV LIMITS, 
%MATERIAL FOV SPACING, ORIENTATION AND MATRIX SIZE

%1) BRING ALL FOV'S TO ACQUISITION CENTERED COORDINATE FRAME
svr.cFOV=(svr.NY-1)/2;
svr.cFOV(4,:,:)=1;
svr.cFOV=matfun(@mtimes,svr.MT,svr.cFOV);%Center of the FOV in spatial coordinates
svr.cFOV=cat(3,svr.cFOV,dynInd(svr.cFOV,1,3));%To bias towards axial view
svr.cFOVm=median(svr.cFOV,3);
svr.cFOV=bsxfun(@minus,svr.cFOV,svr.cFOVm);
svr.cFOV=dynInd(svr.cFOV,1:svr.NV,3);
if ~svr.ParSVR.Prerun;svr.MT(:,4,:)=bsxfun(@minus,svr.MT(:,4,:),svr.cFOV);end

%2) BRING ALL FOV'S TO BRAIN CENTERED COORDINATE FRAME
svr.cEll=permute(svr.Elp(:,1:3,:),[2 1 3]);
svr.cEll=svr.cEll-1;%Center of the brain-Material coordinates in pixels
svr.cEll(4,:,:)=1;
svr.cEll=matfun(@mtimes,svr.MT,svr.cEll);%Brain center in spatial coordinates
%svr.cEll=cat(3,svr.cEll,dynInd(svr.cEll,1,3));%To bias towards axial view
svr.cEllm=mean(svr.cEll,3);
svr.cEll=bsxfun(@minus,svr.cEll,svr.cEllm);
%svr.cEll=dynInd(svr.cEll,1:svr.NV,3);
if ~svr.ParSVR.Prerun;svr.MT(:,4,:)=bsxfun(@minus,svr.MT(:,4,:),svr.cEll);end

%%INCORPORATE THESE IDEAS TO WRITEDATA
%svr.MT=num2cell(svr.MT,1:2);
%svr.MS=permute(svr.MS,[2 1 3]);
%svr.MS=num2cell(svr.MS,1:2);
%svrWriteData(svr,'Al',svr.y,1,[],'',0);
%svr.MS=cat(3,svr.MS{:});svr.MS=permute(svr.MS,[2 1 3]);
%svr.MT=cat(3,svr.MT{:});

%3) GENERIC DEFINITION OF RECONSTRUCTION FOV
svr.FOV=single(zeros(3,2,svr.NV));
svr.FOV(:,2,:)=svr.NY-1;
indV=1:2^3;
subV=ind2subV(2*ones(1,3),indV);
svr.FOV=indDim(svr.FOV,subV',2);%FOV corners of different stacks
svr.FOV(4,:,:)=1;
svr.FOV=matfun(@mtimes,svr.MT,svr.FOV);%Dimensions 4-8-NV
svr.MSS=min(svr.MS(:));%Resolution as the absolute minimum of all resolutions
if ~isempty(svr.ParSVR.MS);svr.MSS(:)=svr.ParSVR.MS;end%Otherwise defined as a parameter
svr.FOVLim=cat(2,min(svr.FOV,[],2),max(svr.FOV,[],2));%FOV limits of different stacks
svr.FOVLim=svr.FOVLim(1:3,:,:);
for v=1:svr.NV;fprintf('FOV center stack %d:%s\n',v,sprintf(' %.2f',mean(svr.FOVLim(:,:,v),2)));end
svr.FullFOVLim=cat(2,min(svr.FOVLim(:,1,:),[],3),max(svr.FOVLim(:,2,:),[],3));%Minimum FOV that contains all stacks
svr.MTT=svr.MSS*eye(4);%Orientation of reconstructions
svr.MTT(4,4)=1;
svr.MSS=svr.MSS*ones(3,1);%Spacing of reconstructions
svr.FullFOVExt=svr.FullFOVLim(:,2)-svr.FullFOVLim(:,1);%Extent of the FOV

%4) DEFINITION OF RECONSTRUCTION FOV CONSIDERING THE PARAMETERS
svr.cEll=permute(svr.Elp(:,1:3,:)-1,[2 1 3]);
svr.cFOV=svr.cEll;
svr.cFOV(4,:,:)=1;%Center of coordinates in pixel coordinates acquisition space
svr.cFOV=matfun(@mtimes,svr.MT(:,:,1),svr.cFOV(:,:,1));%Center of coordinates in world coordinates
if isfield(svr.ParSVR,'FOVSize');svr.maxFOVEncode=svr.ParSVR.FOVSize*ones(3,1);else svr.maxFOVEncode=svr.FullFOVExt;end
svr.FullFOVExtEncode=min(svr.maxFOVEncode,svr.FullFOVExt);%Extension of reconstruction, think not used
svr.FullFOVLimEncode=bsxfun(@plus,svr.cFOV(1:3,:),cat(2,-svr.FullFOVExtEncode/2,svr.FullFOVExtEncode/2));
if isfield(svr.ParSVR,'FOVSize');svr.maxFOV=svr.ParSVR.OverEncode*svr.ParSVR.FOVSize*ones(3,1);else svr.maxFOV=svr.FullFOVExt;end%Maximum FOV
svr.FullFOVExt=min(svr.maxFOV,svr.FullFOVExt);%Extension of the FOV for prescribed reconstruction FOV
svr.FullFOVLim=bsxfun(@plus,svr.cFOV(1:3,:),cat(2,-svr.FullFOVExt/2,svr.FullFOVExt/2));


%5) SIZES OF RECONSTRUCTION SPACE AND UPDATED HOMOGENEOUS COORDINATES
svr.NN=round((svr.FullFOVLim(:,2)-svr.FullFOVLim(:,1))./svr.MSS);%Size of the reconstruction FOV
svr.NNEncode=round((svr.FullFOVLimEncode(:,2)-svr.FullFOVLimEncode(:,1))./svr.MSS);
svr.MTT(1:3,4)=svr.FullFOVLim(:,1);
svr.MSS=(svr.FullFOVLim(:,2)-svr.FullFOVLim(:,1))./svr.NN;
svr.MTT(1:3,1:3)=svr.MSS(1)*eye(3);
fprintf('Reconstructed FOV:%s\n',sprintf(' %.2f',svr.FullFOVExt));
fprintf('Reconstructed resolution:%s\n',sprintf(' %.2f',svr.MSS));
fprintf('Reconstructed grid size:%s\n',sprintf(' %d',svr.NN));

%TRANSFORMS FROM SPATIAL TO MATERIAL
svr.PreTr=isfield(svr,'TV');%Whether a transformation for stacks already exists
svr.Nm=round(bsxfun(@rdivide,bsxfun(@times,svr.NY,svr.MS),svr.MSS));%Size of acquired data when resampled to reconstructed resolution
svr.vA=cell(svr.NV,3);
if ~svr.PreTr;svr.TV=single(zeros([1 1 1 svr.NV 6]));end%Volume transforms
perma=[2 3 4 5 1];

%1) Rotation matrix, angles of rotation and rotation grid
svr.AcqToRec=matfun(@mldivide,svr.MTT,svr.MT);
svr.AcqToRecRot=svr.AcqToRec;
svr.AcqToRecRot(1:3,4,:)=0;
svr.AcqToRecRotNorm=bsxfun(@rdivide,svr.AcqToRecRot,sqrt(sum(svr.AcqToRecRot.^2,1)));
svr.AcqToRecRotNorm=svr.AcqToRecRotNorm(1:3,1:3,:);

%2) Coordinates of acquired data in reconstructed resolution
svr.cFOVRec=matfun(@mldivide,svr.MTT,svr.cFOV);
svr.cFOVRec=svr.cFOVRec(1:3)+1;%+1;%Center of coordinates in pixels, reconstruction space
svr.cFOVAcq=matfun(@mldivide,svr.MT,svr.cFOV);%Center of coordinates in pixel coordinates, acquisition space
svr.cFOVAcq=svr.cFOVAcq(1:3,:,:);%Center of coordinates in pixels, acquisition space
svr.cFOVAcq=bsxfun(@times,svr.cFOVAcq,bsxfun(@rdivide,svr.MS,svr.MSS));%Center of coordinates in pixels at the reconstructed resolution
svr.cFOVAcq=svr.cFOVAcq+1;%+1;
svr.IndAcq=svr.cFOVAcq-svr.cFOVRec;%Starting index on the acquisition space (referred to 0)
svr.IndAcqRound=round(svr.IndAcq);
svr.IndAcqShift=svr.IndAcq-svr.IndAcqRound;%This shift should be applied as a translation after rotation, but we haven't coded it, we let the method learn it
svr.IndAcqShift=bsxfun(@times,svr.IndAcqShift,svr.MSS);
svr.IndAcqShift=matfun(@mldivide,svr.AcqToRecRotNorm,svr.IndAcqShift);
svr.NmReal=zeros(3,1,svr.NV);svr.NYReal=zeros(3,1,svr.NV);

for v=1:svr.NV
    %3) CROP TO VALID INDEXES
    for s=1:3
        svr.vAcq{v}{s}=svr.IndAcqRound(s,1,v)+1:svr.IndAcq(s,1,v)+svr.NN(s);
        svr.vRec{v}{s}=find(svr.vAcq{v}{s}>=1 & svr.vAcq{v}{s}<=svr.Nm(s,1,v));
        svr.vAcq{v}{s}=svr.vAcq{v}{s}(svr.vRec{v}{s});
        if isfield(svr,'Mtr')
            if svr.id(v,s)~=3
                svr.vAcqY{v}{s}=round((svr.vAcq{v}{s}(1)-1)*(svr.NY(s,1,v)/svr.Nm(s,1,v))+1);
                svr.vAcqY{v}{s}=svr.vAcqY{v}{s}+round(svr.NY(s,1,v)*length(svr.vAcq{v}{s})/svr.Nm(s,1,v));
                svr.NmReal(s,1,v)=length(svr.vAcq{v}{s});
                svr.vAcq{v}{s}=1:svr.NmReal(s,1,v);
            else
                svr.vAcqY{v}{s}=1:NY(s,1,v);
                svr.NmReal(s,1,v)=svr.Nm(s,1,v);
            end
        else
            svr.NmReal(s,1,v)=svr.Nm(s,1,v);
        end
    end
    
    if isfield(svr,'Mtr')
        svr.y{v}=dynInd(svr.y{v},svr.vAcqY{v},1:3);        
        svr.El{v}=dynInd(svr.El{v},svr.vAcqY{v},1:3);
        svr.NYReal(:,1,v)=size(svr.y{v})';
    else
        svr.NYReal(:,1,v)=svr.NY(:,1,v);
    end
    
    %4) TRANSFORMS
    if ~svr.PreTr
        svr.TP{v}=zeros(1,1,1,1,svr.P(v),6);%Package transforms
        for p=1:svr.P(v);svr.TE{v}{p}=zeros(1,1,1,1,1,svr.slPerPack{v}(p),6);end%Slice transforms
        svr.TV(1,1,1,v,4:6)=-convertRotation(permute(SpinCalc('DCMtoEA321',svr.AcqToRecRotNorm(:,:,v)',1e-3,0)',perma),'deg','rad');%Convert axes rotation matrix to Euler angles  
        svr.TV(1,1,1,v,1:3)=permute(svr.IndAcqShift(:,1,v),perma);
        %svr.TV(1,1,1,v,4:6)=-convertRotation(permute(SpinCalc('DCMtoEA321',svr.AcqToRecRotNorm(:,:,v),1e-3,0)',perma),'deg','rad');%Convert axes rotation matrix to Euler angles-This would be the rotation when applied the opposite way  
    end
end

%SLICE WEIGHTS/SLICE PROFILES/APODIZATION/VALID REGIONS
svr.NNEncode=svr.NNEncode';svr.NN=svr.NN';svr.NY=permute(svr.NY,[2 1 3]);svr.Nm=permute(svr.Nm,[2 1 3]);svr.NYReal=permute(svr.NYReal,[2 1 3]);svr.NmReal=permute(svr.NmReal,[2 1 3]);svr.MS=permute(svr.MS,[2 1 3]);svr.MSS=svr.MSS';svr.cFOVRec=svr.cFOVRec';

%Slice weights
if ~isfield(svr,'W')
    svr.W=cell(1,svr.NV);
    for v=1:svr.NV
        NW=ones(1,3);NW(svr.id(v,3))=svr.NY(1,svr.id(v,3),v);
        svr.W{v}=real(ones(NW,'like',svr.y{v}));
    end
end

svr.H=cell(1,svr.NV);%Apodization info
%svr.xValid=zeros([svr.NN svr.NV],'like',svr.y{1});
for v=1:svr.NV     
    %Slice profiles
    if svr.ParSVR.SlBefore==0;NSlPr=svr.NY(1,svr.id(v,3),v);else NSlPr=svr.Nm(1,svr.id(v,3),v);end        
    sigma=svr.ParSVR.FWHM/(2*sqrt(2*log(2)));%For FWHM as a ratio of slice thickness and slice separation    
    if svr.ParSVR.SlBefore==1;sigma=sigma*svr.Nm(1,svr.id(v,3),v)/svr.NY(1,svr.id(v,3),v);end
    NSl=ones(1,3);NSl(svr.id(v,3))=2*NSlPr;
    kGrid=generateGrid(NSl,gpu,pi,ceil((NSl+1)/2));
    svr.SlPr{v}=exp(-(kGrid{svr.id(v,3)}.^2)*(2*sigma^2));    
    svr.SlPr{v}=ifftshift(svr.SlPr{v},svr.id(v,3));     
    svr.SlPr{v}=sqrt(2*NSlPr)*svr.SlPr{v}/norm(svr.SlPr{v}(:));    
    svr.SlPr{v}=dynInd(svr.SlPr{v},1:NSlPr,svr.id(v,3));    
    
    %Apodization
    svr.H{v}=fftshift(buildFilter(svr.NYReal(1,:,v),'tukey',ones(1,3),gpu,svr.ParSVR.UseApo*ones(1,3)));%Not fully symmetric   
    
    %%Valid regions-NOT USED!!!
    %for d=1:3;svr.xValid(:,:,:,v)=dynInd(svr.xValid(:,:,:,v),svr.vRec{v}{d},d,1);end%Areas where a given stack does overlap 
    %dist=4;%4mm apodization
    %svr.xValid=morphFourier(morphFourier(svr.xValid,-dist*ones(1,3),svr.MSS,zeros(1,3)),dist*ones(1,3),svr.MSS,zeros(1,3),1);%First erosion then distance function for soft masking including only voxels inside
end

%BUILD THE TRANSFORM GRIDS
%%%THIS SHOULD BE DEFINED IN THE REAL DOMAIN FOR FASTER PERFORMANCE!!!
[svr.rGrid,svr.kGrid,svr.rkGrid,~,svr.cGrid]=generateTransformGrids(svr.NNEncode.*svr.MSS,gpu,svr.NNEncode,svr.cFOVRec/svr.ParSVR.OverEncode,1);
[svr.FT,svr.FTH]=buildStandardDFTM(svr.NNEncode,0,gpu);
[svr.rGridV,svr.kGridV,svr.rkGridV,~,svr.cGridV]=generateTransformGrids(svr.NN.*svr.MSS,gpu,svr.NN,svr.cFOVRec,1);
[svr.FTV,svr.FTHV]=buildStandardDFTM(svr.NN,0,gpu);
svr.FY=cell(1,svr.NV);svr.FYH=cell(1,svr.NV);
for v=1:svr.NV;[svr.FY{v},svr.FYH{v}]=buildStandardDFTM(svr.NYReal(1,:,v),0,gpu,1);end

%RECONSTRUCTED DATA AND STRUCTURES FOR RECONSTRUCTION
if isfield(svr,'x');svr.x=resampling(svr.x,svr.NN,0,ones(1,3));else svr.x=zeros(svr.NN,'like',svr.y{1});end
svr.xx=[];%MARKING THAT WE'VE JUST CAME FROM HERE
svr.yy=svr.y;
svr.E=svr.y;%Residuals
svr.My=svr.y;

%REGULARIZER
%%%THIS SHOULD BE DEFINED IN THE REAL DOMAIN!!!
if ~isempty(svr.ParSVR.regFracOrd)
    NG=length(svr.ParSVR.regFracOrd);
    for g=1:NG
        svr.F{g}=buildFilter(2*svr.NN,'FractionalFiniteDiscreteIsoNorm',svr.ParSVR.spti(g)*ones(1,3),gpu,svr.ParSVR.regFracOrd(g),1);
        if svr.ParSVR.regFracOrd(g)~=0;svr.F{g}(1)=0;end%NEW, THIS WAS INTRODUCING LOTS OF INSTABILITIES
    end
end

%SHEARLET
if svr.ParSVR.tiSh~=0    
    J=3;
    svr.Sh=buildShearlet(svr.NN,J,0,ones(1,J),'kos');
    NSh4=size(svr.Sh.S,4);
    fprintf('Size shearlets:%s %d\n',sprintf(' %d',svr.NN),NSh4);
end

%PARAMETERS FOR MOTION CORRECTION
svr.M=svr.x;svr.M(:)=1;
svr.MHist=[];
svr.MM=svr.yy;

%%%THIS SHOULD BE DEFINED IN THE REAL DOMAIN!!!
if svr.ParSVR.fracOrd~=0
    svr.G=cell(1,svr.NV);
    for v=1:svr.NV
        NG=svr.NYReal(1,:,v);
        NG(svr.id(v,3))=1;%%%THIS IS NEW---IT WAS IN 3D!!!!
        svr.G{v}=buildFilter(NG,'FractionalFiniteDiscreteIso',ones(1,3),gpu,svr.ParSVR.fracOrd);
    end%For fractional finite difference motion estimation    
end

%TO WRITE SOME INFORMATION FOR ILLUSTRATION
if svr.ParSVR.Prerun
    %SOFT MASK
    svr=svrGenerateSoftMask(svr);

    %WRITE AVERAGED INFORMATION FOR DEBUGGING
    svr=svrEmpiricalPseudoInverseAveraging(svr);
end
