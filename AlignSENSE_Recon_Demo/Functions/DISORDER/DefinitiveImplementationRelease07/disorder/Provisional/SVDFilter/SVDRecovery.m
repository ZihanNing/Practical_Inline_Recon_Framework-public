function [rec,gamma]=SVDRecovery(rec)

% SVDRECOVERY recovers the signal based on SVD shrinkage in local patches
%   REC=SVDRECOVERY(REC,{EXTRUN})
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   ** REC is a reconstruction structure with recovered signal information 
%   ** GAMMA is an estimate of the patch size
%

beta=[];
if rec.Fail || ~rec.Alg.SVDRecover;buildSoftMask;return;end
gpu=isa(rec.x,'gpuArray');

%To check the value of Beta that is being used:
if rec.Dyn.Debug>=1;fprintf('Used shape factor:%s\n',sprintf(' %.2f',rec.Alg.parR.Gamma));end

%DATA DIMENSIONS, SIZES, STRUCTURE AND PERHAPS NOISE GENERATION
N=size(rec.x);N(end+1:4)=1;
if prod(N(4:end))<30;rec.Alg.SVDRecover=0;buildSoftMask;return;end%Not enough samples to filter
NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
N(4:5)=[NE N(4)/NE];
rec.x=reshape(rec.x,N);
NPE=length(rec.Par.Mine.pedsUn);
if rec.Alg.PlugNoise==6
    rec.x=plugNoise(rec.x);
    rec.G(:)=1;
end
%REMOVE PARTIAL FOURIER AND REGRID FROM PARTIAL FOURIER. WE ASSUME BOTH ECHOES HAVE THE SAME PF STRUCTURE
cov=[];
perm=1:5;perm(4)=5;perm(5)=4;

if rec.Alg.parR.illustration 
    sl=25;
end

if ismember(rec.Alg.parR.HalfScanCorrection,[3 5])
    for n=[8 12];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    rec=PEregrid(rec,0);    
    for n=[8 12];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
    if NPE>1
        for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.Nat,:);end
    end
    if rec.Alg.MargosianFilter
        [rec.r,indM]=margosianFilter(rec.x,Enc,1,0);
    end    
    if rec.Alg.parR.illustration==1
        aux=dynInd(rec.r,[sl 1 2],3:5);
        fig01(aux,1);
    end        
    
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end    
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,18);
    rec=PEregrid(rec);
    rec.Dyn.Typ2Rec(ismember(rec.Dyn.Typ2Rec,18))=[];
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    g=abs(dynInd(rec.G,2,5));
    if rec.Alg.parR.illustration==1
        aux=dynInd(g,sl,3);
        fig01(aux,2);
    end
elseif ismember(rec.Alg.parR.HalfScanCorrection,[1 4])
    for n=[8 12];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    rec=PEregrid(rec,0);    
    for n=[8 12];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
    if NPE>1
        for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.Nat,:);end
    end
    if rec.Alg.MargosianFilter        
        [rec.r,indM]=margosianFilter(rec.x,Enc,1,0);                
        M=single(abs(margosianRegrid(dynInd(rec.M,{rec.Par.Mine.pedsUn,1},4:5),indM))>0.5);
        rec.r=margosianRegrid(rec.r,indM);
        g=M.*abs(margosianRegrid(dynInd(rec.G,2,5),indM,1,1));        
    else
        rec.r=rec.x;
        g=abs(dynInd(rec.G,2,5));
    end
end
if ~ismember(rec.Alg.parR.HalfScanCorrection,[1 3 4 5])
    rec.r=rec.x;
    g=abs(dynInd(rec.G,2,5));
end
mg=mean(g(g>1e-3));%Mean noise
rec.x=gather(rec.x);

NRes=size(rec.r);NRes(end+1:5)=1;
NDS=numDims(g);NDS=min(NDS,3);
voxsiz=rec.Par.Scan.AcqVoxelSize(1:NDS);
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
voxsiz=voxsiz.*N(1:NDS)./NRes(1:NDS);
if rec.Dyn.Debug>=1;fprintf('Voxel size:%s\n',sprintf(' %.2f',voxsiz));end

%COVARIANCE COMPUTATION
if ismember(rec.Alg.parR.HalfScanCorrection,1:6)    
    Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
    if isfield(rec.Par.Mine,'StrFactorMax')
        for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(1,:);end
    end    
    Enc.AcqSize=repmat(Enc.AcqSize,[NPE 1]);    
    for p=1:NPE
        permr=(1:2)';
        if ~isempty(rec.Par.Mine.Signs);permr(1:2)=rec.Par.Mine.Signs(:,:,rec.Par.Mine.pedsUn(p))*permr(1:2);end
        if NPE>1;Enc.AcqSize(p,1:2)=Enc.AcqSize(rec.Par.Mine.Nat,abs(permr)');else Enc.AcqSize(p,1:2)=Enc.AcqSize(1,abs(permr)');end
    end    
    Enc.N=N;
    gaux=g;        
    if rec.Alg.parR.GFactorCorrection==1;gaux(:)=1;end      
    if ismember(rec.Alg.parR.HalfScanCorrection,[2 6]);cov=margosianCovariance(Enc,gaux,'ramp',gpu);
    elseif ismember(rec.Alg.parR.HalfScanCorrection,[3 5]);cov=margosianCovariance(Enc,gaux,'zefi',gpu);
    elseif ismember(rec.Alg.parR.HalfScanCorrection,[1 4]);cov=margosianCovariance(Enc,gaux,'iden',gpu);
    end    
end

%ESTIMATE AND CORRECT FOR LINEAR PHASE
if rec.Alg.parR.LinearPhaseCorrection
    fl=ridgeDetection(dynInd(rec.r,1,4),1:2);
    rec.r=bsxfun(@times,rec.r,conj(fl));
    fl=gather(fl);
end
if rec.Alg.parR.illustration==1
    aux=dynInd(rec.r,[sl 1 2],3:5);
    fig01(aux,3);
    aux=dynInd(fl,[sl 1 2],3:5);
    fig01(aux,4);
end

%NON LINEAR PHASE CORRECTION
f=rec.r;%Magnitude denoising
if rec.Alg.parR.FilteredPhaseCorrection==1;f=filterLow(f);%LPF
elseif rec.Alg.parR.FilteredPhaseCorrection==2;f=filterShearlet(f);%Shearlets
elseif rec.Alg.parR.FilteredPhaseCorrection==3;f(:)=1;%Only linear
end
f=sign(f);
rec.r=bsxfun(@times,rec.r,conj(f));
if ~rec.Alg.parR.UseComplexData;rec.r=real(rec.r);end%Real information (it can be magnitude if FilterPhaseCorrection=0 - Rician bias / real part if FilterPhaseCorrection=3: Imaginary residual bias)
f=gather(f);

%STANDARDIZE NOISE
M=g;M=single(M>1e-3);
indPEV=[];
if NPE>1
    for p=1:NPE;indPEV=vertcat(indPEV,find(rec.Par.Mine.diInfo(1:N(5),5)==rec.Par.Mine.pedsUn(p)));end
    rec.r=dynInd(rec.r,indPEV,5);        
    rec.r=reshape(rec.r,[NRes(1:3) NRes(5)/NPE NPE*NRes(4)]);
    rec.r=permute(rec.r,perm);
end
if rec.Alg.parR.GFactorCorrection==1
    g(g<=1e-3)=mg;
    rec.r=bsxfun(@times,rec.r,real(M./g));
end

%ADD NOISE OUTSIDE THE MASK SO THE PROBLEM BECOMES INDEPENDENT FROM MASKING
if rec.Alg.parR.GFactorCorrection~=2 || ~strcmp(rec.Alg.parR.NoiEstMeth,'None')
    if rec.Alg.parR.GFactorCorrection==1;no=plugNoise(rec.r);else no=mg*plugNoise(rec.r);end
    if ismember(rec.Alg.parR.HalfScanCorrection,[2:3 5:6])
        Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
        for p=1:NPE
            permr=(1:2)';
            if ~isempty(rec.Par.Mine.Signs);permr(1:2)=rec.Par.Mine.Signs(:,:,rec.Par.Mine.pedsUn(p))*permr(1:2);end
            if NPE>1;Enc.AcqSize(1,1:2)=rec.Enc.AcqSize(rec.Par.Mine.Nat,abs(permr)');else Enc.AcqSize(1,1:2)=rec.Enc.AcqSize(1,abs(permr)');end
            for l=1:length(Enc.kRange);Enc.kRange{l}=rec.Enc.kRange{l}(rec.Par.Mine.pedsUn(p),:);end            
            if rec.Alg.MargosianFilter
                if ismember(rec.Alg.parR.HalfScanCorrection,[2 6]);no=dynInd(no,NE*(p-1)+1:p*NE,4,margosianFilter(dynInd(no,NE*(p-1)+1:p*NE,4),Enc,1,1,'ramp'));
                else no=dynInd(no,NE*(p-1)+1:p*NE,4,margosianFilter(dynInd(no,NE*(p-1)+1:p*NE,4),Enc,1,1,'zefi'));
                end
            end
        end
    end
    if ~rec.Alg.parR.UseComplexData;no=abs(no);end
    rec.r=bsxfun(@times,rec.r,M)+bsxfun(@times,no,1-M);
    no=[];
end

%SVD PATCH BASED RECOVERY
if ismember(rec.Alg.parR.HalfScanCorrection,[1 5:6])
    for l=1:NPE        
        indL=NE*(l-1)+1:l*NE;
        covL=cov;covL.dimPE=cov.dimPE(l);covL.dimNoPE=cov.dimNoPE(:,l);covL.Lambda=dynInd(cov.Lambda,l,4);
        [r,s,p,a,b]=patchSVShrinkage(dynInd(rec.r,indL,4),3,voxsiz,rec.Alg.parR,covL,[],[],rec.Alg.parV);
        rec.r=dynInd(rec.r,indL,4,r);gamma(l)=b;
        if ~isfield(rec,'s');rec.s=s;rec.p=p;rec.a=a;else rec.s=cat(4,rec.s,s);rec.p=cat(4,rec.p,p);rec.a=cat(4,rec.a,a);end
%        if ~isfield(rec,'s');rec.s=dynInd(rec.r,[1 1],4:5);rec.p=rec.s;rec.a=rec.s;else rec.s=cat(4,rec.s,dynInd(rec.r,[1 1],4:5));rec.p=cat(4,rec.p,dynInd(rec.r,[1 1],4:5));rec.a=cat(4,rec.a,dynInd(rec.r,[1 1],4:5));end
    end
elseif rec.Alg.parR.HalfScanCorrection==4
    assert(NE==1,'Paired PE denoising not implemented for multiple echo data');
    assert(mod(NPE,2)==0,'Paired PE denoising not implemented for odd number of PEs');
    rec.r=dynInd(rec.r,2:2:NPE,4,flip(dynInd(rec.r,2:2:NPE,4),2));
    for s=1:2;cov.Lambda=dynInd(cov.Lambda,2:2:NPE,4,flip(dynInd(cov.Lambda,2:2:NPE,4),s));end  
    for l=1:NPE/2
        indL=2*l-1:2*l;
        covL=cov;covL.dimPE=cov.dimPE(indL);covL.dimNoPE=cov.dimNoPE(:,indL);covL.Lambda=dynInd(cov.Lambda,indL,4);
        [r,s,p,a,b]=patchSVShrinkage(dynInd(rec.r,indL,4),3,voxsiz,rec.Alg.parR,covL,[],[],rec.Alg.parV);
        rec.r=dynInd(rec.r,indL,4,r);gamma(l)=b;
        s=repmat(s,[1 1 1 2]);p=repmat(p,[1 1 1 2]);a=repmat(a,[1 1 1 2]);
        if ~isfield(rec,'s');rec.s=s;rec.p=p;rec.a=a;else rec.s=cat(4,rec.s,s);rec.p=cat(4,rec.p,p);rec.a=cat(4,rec.a,a);end
%        if ~isfield(rec,'s');rec.s=dynInd(rec.r,{1:2 1},4:5);rec.p=rec.s;rec.a=rec.s;else rec.s=cat(4,rec.s,dynInd(rec.r,{1:2,1},4:5));rec.p=cat(4,rec.p,dynInd(rec.r,{1:2,1},4:5));rec.a=cat(4,rec.a,dynInd(rec.r,{1:2,1},4:5));end
    end
    for n=18:21;datTyp=rec.Plan.Types{n};
        rec.(datTyp)=dynInd(rec.(datTyp),2:2:NPE,4,flip(dynInd(rec.(datTyp),2:2:NPE,4),2));
    end
else%THIS IS USED FOR NEONATES
    [rec.r,rec.s,rec.p,rec.a,gamma]=patchSVShrinkage(rec.r,3,voxsiz,rec.Alg.parR,cov,[],[],rec.Alg.parV);
    %rec.s=abs(dynInd(rec.r,[1 1],[4 5]));rec.p=rec.s;rec.a=rec.s;
end
rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,(18:21)');rec.Dyn.Typ2Wri(18:21)=1;

%DESTANDARDIZE NOISE
if rec.Alg.parR.GFactorCorrection==1;rec.r=bsxfun(@times,rec.r,g);
elseif rec.Alg.parR.GFactorCorrection~=2 || ~strcmp(rec.Alg.parR.NoiEstMeth,'None');rec.r=bsxfun(@times,rec.r,M);
end

if NPE>1
    rec.r=permute(rec.r,perm);
    rec.r=reshape(rec.r,NRes);
    rec.r=dynInd(rec.r,indPEV,5,rec.r);
end

%REINCORPORATE NON LINEAR PHASE
if gpu;f=gpuArray(f);end
rec.r=bsxfun(@times,rec.r,f);f=[];
%REINCORPORATE LINEAR PHASE
if rec.Alg.parR.illustration==1     
    aux=dynInd(rec.r,[sl 1 2],3:5);
    fig01(aux,5);
end  
if rec.Alg.parR.LinearPhaseCorrection
    if gpu;fl=gpuArray(fl);end
    rec.r=bsxfun(@times,rec.r,fl);fl=[];
end
if rec.Alg.parR.illustration==1
    aux=dynInd(rec.r,[sl 1 2],3:5);
    fig01(aux,6);
end
if gpu;rec.x=gpuArray(rec.x);end

%%PARTIAL FOURIER, REGRID, FILTER AND OBTAIN RESIDUALS 
if ismember(rec.Alg.parR.HalfScanCorrection,[1 4])
    for n=18:21;datTyp=rec.Plan.Types{n};
        if n==18       
            rec.(datTyp)=margosianRegrid(rec.(datTyp),indM,0);
            Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
            if NPE>1
                for m=1:length(Enc.kRange);Enc.kRange{m}=Enc.kRange{m}(rec.Par.Mine.Nat,:);end
            end
            if rec.Alg.MargosianFilter;rec.(datTyp)=margosianFilter(rec.(datTyp),Enc);end
        else       
            rec.(datTyp)=margosianRegrid(rec.(datTyp),indM,0,2);
            rec.(datTyp)=real(rec.(datTyp));
        end
    end
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    rec=PEregrid(rec);
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
elseif ismember(rec.Alg.parR.HalfScanCorrection,[3 5])
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    rec.Dyn.Typ2Rec(ismember(rec.Dyn.Typ2Rec,19:21))=[];
    rec=PEregrid(rec,0);
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end    
    Enc=rec.Enc;for l=1:length(Enc.kRange);Enc.kRange{l}=Enc.kRange{l}(rec.Par.Mine.pedsUn,:);Enc.rRange{l}=Enc.rRange{l}(:,:,rec.Par.Mine.pedsUn);end;Enc.MPSOffcentres=Enc.MPSOffcentres(:,:,rec.Par.Mine.pedsUn);
    if NPE>1
         for m=1:length(Enc.kRange);Enc.kRange{m}=Enc.kRange{m}(rec.Par.Mine.Nat,:);end
    end
    if rec.Alg.MargosianFilter;rec.r=margosianFilter(rec.r,Enc);end
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
    rec=PEregrid(rec);
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,(19:21)');
    for n=[8 12 18];datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end            
end

rec.n=rec.r-rec.x;
rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,22);rec.Dyn.Typ2Wri(22)=1;


%BUILDSOFTMASK
buildSoftMask;

%FINISH ROI EXTRACTION
if rec.Alg.parR.ExtractROI
    typ2Rec=rec.Dyn.Typ2Rec;
    for n=typ2Rec';datTyp=rec.Plan.Types{n};
        rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,2);
    end
end

function buildSoftMask
    %BUILD SOFT MASK
    NDims=min(numDims(rec.M),3);
    distan=10*ones(1,NDims);
    mirr=ones(1,NDims);
    voxsiz=rec.Par.Scan.AcqVoxelSize(1:NDims);
    rec.M=dynInd(rec.M,2,5);
    rec.M=morphFourier(rec.M,distan,voxsiz,mirr,1);%Soft masking
    rec.M(rec.M>1-1e-6)=1;rec.M(rec.M<1e-6)=0;
end

function x=filterLow(x)
    x=extractROI(x,rec.Enc.ROI,0,1);
    NN=size(x);
    H=buildFilter(NN(1:2),'tukeyIso',1/2,gpu,1);%Hardcoded to 1/6th of the resolution on the basis of no artifacts on the adult scan
    x=filtering(x,H);
    x=extractROI(x,rec.Enc.ROI,1,1);
end

function x=filterShearlet(x)
    ups=1;
    J=3+log2(ups);
    x=extractROI(x,rec.Enc.ROI,0,1);
    g=extractROI(g,rec.Enc.ROI,0,1);
    NN=size(x);NN(end+1:5)=1;
    sH=buildShearlet(ups*NN(1:2),J,gpu,ones(1,J),'kos');
    g=resampling(g,ups*NN(1:2),2);
    for qq=1:NN(5)
        for pp=1:NN(4)
            xV=resampling(dynInd(x,[pp qq],4:5),ups*NN(1:2),2);
            xV=shearletFilter(xV,sH,g,0.05,[16 32]);            
            x=dynInd(x,[pp qq],4:5,resampling(xV,NN(1:2),2));
        end        
    end
    g=resampling(g,NN(1:2),2);
    x=extractROI(x,rec.Enc.ROI,1,1);
    g=extractROI(g,rec.Enc.ROI,1,1);
end

end
