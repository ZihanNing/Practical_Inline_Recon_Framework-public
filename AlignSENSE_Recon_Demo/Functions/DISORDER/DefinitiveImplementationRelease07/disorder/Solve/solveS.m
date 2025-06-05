function rec=solveS(rec)

%SOLVES   Estimates the sensitivities
%   REC=SOLVES(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and mask rec.M
%

if isfield(rec,'Names') && isfield(rec.Names,'Name');name=rec.Names.Name;else name='NULL';end
gpu=isa(rec.y,'gpuArray');%It was rec.Dyn.GPU
if gpu;gpuF=2;else gpuF=0;end
if ~isfield(rec,'Alg') || ~isfield(rec.Alg,'parS')
    %PARAMETERS FOR COIL AND MASK ESTIMATION
    rec.Alg.parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
    rec.Alg.parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
    rec.Alg.parS.Otsu=[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)    
    rec.Alg.parS.lambda=20;%This parameter has been constantly changing, toddlers have been processed with it set to 1, but it was introducing noise for neonates, 2 may suffice, but 5 gives a bit of margin
    rec.Alg.parS.order=2;%Regularization order for coil estimation-Something around 0.2 is optimal for order 2
    rec.Alg.parS.nErode=0;%2;%Erosion for masking (in mm)
    rec.Alg.parS.nDilate=3;%Dilation for masking (in voxels then in mm)
    rec.Alg.parS.conComp=2;%Whether to get the largest connected component after erosion (1) and to fill holes (2)
    rec.Alg.parS.GibbsRingi=[0 0];%Gibbs ringing for coil profiles
    rec.Alg.parS.ResolRatio=[1 1];%Resolution ratio for coil profiles
    rec.Alg.parS.TolerSolve=1e-5;%S-solver tolerance
    rec.Alg.parS.nIt=300;%S-solver maximum number of iterations
    rec.Alg.parS.softMaskFactor=10;%Factor to normalize the body coil data for soft masking
end
parS=rec.Alg.parS;

if ~isfield(rec,'Plan') || ~isfield(rec.Plan,'NDims');rec.Plan.NDims=14;end
if ~isfield(rec,'Par');rec.Par=[];end
if ~isfield(rec.Par,'Labels') || ~isfield(rec.Par.Labels,'CoilNrsPerStack');rec.Par.Labels.CoilNrsPerStack=[];end
if ~isfield(rec.Par,'Mine') || ~isfield(rec.Par.Mine,'pedsUn');rec.Par.Mine.pedsUn=1;end
if ~isfield(rec.Par.Mine,'Modal');rec.Par.Mine.Modal=2;end
if ~isfield(rec.Par.Labels,'SliceGaps');rec.Par.Labels.SliceGaps=0;end
if ~isfield(rec,'Enc') || ~isfield(rec.Enc,'AcqVoxelSize');rec.Enc.AcqVoxelSize=ones(1,3);end
if ~isfield(rec.Enc,'FOVSize');rec.Enc.FOVSize=size(rec.y);rec.Enc.FOVSize(end+1:3)=1;rec.Enc.FOVSize=rec.Enc.FOVSize(1:3);end
if ~isfield(rec.Enc,'RecSize')
    if isfield(rec,'x');rec.Enc.RecSize=size(rec.x);else rec.Enc.RecSize=size(rec.y);end
    rec.Enc.RecSize(end+1:3)=1;rec.Enc.RecSize=rec.Enc.RecSize(1:3);
end
if ~isfield(rec,'Dyn') || ~isfield(rec.Dyn,'Typ2Rec');rec.Dyn.Typ2Rec=[];end
if ~isfield(rec.Alg,'OverDec');rec.Alg.OverDec=ones(1,3);end
    
%DECIDE ON BODY/SURFACE COIL ROLS
if ~isempty(rec.x) && isfield(rec,'y')    
    D=rec.x;Z=rec.y;%Body and surface coils
    NDY=numDims(rec.y);    
    if NDY>4
        perm=1:rec.Plan.NDims;perm(4)=NDY;perm(NDY)=4;
        D=permute(D,perm);
        perm=1:rec.Plan.NDims;perm(NDY)=5;perm(5)=NDY;
        D=permute(D,perm);Z=permute(Z,perm);
    end
    if isfield(rec,'S');S=rec.S;rec=rmfield(rec,'S');end
    %This is to recompute the mask!
    if isfield(rec,'M');rec=rmfield(rec,'M');end
elseif isfield(rec,'y') && length(rec.Par.Labels.CoilNrsPerStack)==2
    Z=dynInd(rec.y,rec.Par.Labels.CoilNrsPerStack{1}+1,4);
    D=dynInd(rec.y,rec.Par.Labels.CoilNrsPerStack{2}+1,4);
else
    rec.Par.Labels.CoilNrsPerStack{1}=0:size(rec.y,4)-1;
    D=compressCoils(rec.y,1);
    Z=rec.y;
    fprintf('NO BODY COIL IDENTIFIED FOR FILE %s',name)
end

%AVERAGE DATA AND COMPRESS THE BODY COIL
N=size(Z);N(end+1:rec.Plan.NDims)=1;
M=size(D);M(end+1:rec.Plan.NDims)=1;
if any(N(5:end)~=M(5:end))
    fprintf('Outer dimensions of surface coil (%s) are not matched with body coil dimensions (%s). Using the first element in the outer array\n',sprintf('%d ',N(5:end)),sprintf('%d ',M(5:end)));
    Z=Z(:,:,:,:,1);M(5:rec.Plan.NDims)=1;
    D=D(:,:,:,:,1);N(5:rec.Plan.NDims)=1;
end
D=resPop(D,5:rec.Plan.NDims,prod(M(5:rec.Plan.NDims)),5);Z=resPop(Z,5:rec.Plan.NDims,prod(N(5:rec.Plan.NDims)),5);
N=size(Z);N(end+1:5)=1;
M=size(D);M(end+1:5)=1;
assert(N(5)==M(5),'Outer dimensions of surface coil (%d) are not matched with body coil dimensions (%d) for file %s',N(5),M(5),name);

%WE USE THE MEDIAN OVER THE REPEATS TO ESTIMATE
NPE=min(length(rec.Par.Mine.pedsUn),N(5));
if NPE==1
    D=multDimMed(D,5:rec.Plan.NDims);
    for s=1:N(4);Z=dynInd(Z,[s 1],4:5,multDimMed(dynInd(Z,s,4),5:rec.Plan.NDims));end
    Z=dynInd(Z,1,5);
end
for p=1:NPE
    Dp=dynInd(D,p,5);Zp=dynInd(Z,p,5);    
    if exist('S','var')%Necessary for SE reference for fMRI in neonates for instance
        if size(S,5)==4;Sp=dynInd(S,rec.Par.Mine.pedsUn(p),5);
        elseif size(S,5)==NPE;Sp=dynInd(S,p,5);
        end
    end
    Dp=compressCoils(Dp,ones(1,3),[],[]);
    if rec.Par.Mine.Modal==2;rec.x=Dp;end

    if ~exist('S','var');Sp=repmat(Dp,[1 1 1 size(Z,4)]);Sp(:)=0;end

    %DECIDE ON RESOLUTION TO OPERATE AT: WE OPERATE AT THE MINIMUM ACQUIRED RESOLUTION BUT ONLY FOR INITIALIZATION
    sp=rec.Enc.AcqVoxelSize;
    sp(3)=sp(3)+rec.Par.Labels.SliceGaps(1);
    N=size(Dp);

    %NORMALIZATION  
    normal=max(abs(Dp(:)));
    %normal=prctile(abs(Dp(:)),95);
    Dp=Dp/normal;Zp=Zp/normal;

    %INITIALIZE COILS
    if ~exist('S','var');Sp=bsxfun(@rdivide,Zp,Dp);end
    
    %MASK
    W=refineMask(Dp,parS,sp(:)'/max(sp(:)));

    %INFORMATION FOR PSEUDOINVERSE COMPUTATION AND ACTUAL COMPUTATION
    E.Xf=Dp;%Encoder
    E.Uf=cell(1,3);for n=1:length(E.Uf);E.Uf{n}.NY=rec.Enc.RecSize(n);E.Uf{n}.NX=rec.Enc.FOVSize(n);end
    [E.UAf,EH.UAb]=buildFoldM(rec.Enc.FOVSize(1:length(E.Uf)),rec.Enc.RecSize(1:length(E.Uf)),gpu);

    EH.Xb=conj(Dp);%Decoder
    EH.Mc=single(W);
    EH.Ub=cell(1,3);for n=1:length(EH.Ub);EH.Ub{n}.NY=rec.Enc.RecSize(n);EH.Ub{n}.NX=rec.Enc.FOVSize(n);end
    
    %Constrain
    C=[];    
    if exist('S','var')%Modification to have a similar phase
        aux=EH.Mc.*EH.Xb;
        EH.Xb=EH.Xb*exp(-1i*angle(mean(aux(:))));
        E.Xf=E.Xf.*exp(1i*angle(mean(aux(:))));
    end    

    %Preconditioner
    A=[];

    if rec.Dyn.Debug>=1;fprintf('Resolution of the reference:%s\n',sprintf(' %.2f',sp));end
    %Regularizer    
    if ~exist('S','var');GF=abs(bsxfun(@rdivide,Dp,Zp));
    else            
        GF=abs(Dp);  
    end        
    %GF=abs(GF).^2;%IT IS GIVING GOOD RESULTS AS IT IS, ALTHOUGH IT IS 'DIMENSIONALLY' ARGUABLE
    NGF=size(GF);
    if ~exist('S','var');H=buildFilter(2*NGF(1:3),'tukeyIso',0.125*sp/5,gpu,1,1);else H=buildFilter(2*NGF(1:3),'tukeyIso',0.25*sp/5,gpu,1,1);end
    GF=abs(filtering(GF,H,1));

    %R.Fo.la=bsxfun(@rdivide,parS.lambda*(mean(sp)/5)^3,GF+1e-2);
    R.Fo.la=parS.lambda;
    %visReconstruction(R.Fo.la)
    R.Fo.Fi=buildFilter(2*N,'FractionalFiniteDiscreteIso',sp,gpu,parS.order,1);R.Fo.mi=1;%[1 1 1];    

    %THIS FORCES SAME RESULTS AS IN THE INPUT
    %R.Tp.la=10000;%bsxfun(@rdivide,10*parS.lambda*(mean(sp)/5)^3,GF+1e-2);
    %R.Tp.x0=Sp;

    %SOLVE
    Sp=CGsolver(Zp,E,EH,A,C,R,Sp,parS.nIt,[],parS.TolerSolve);

    %ASSIGN THE COMPUTED INFORMATION TO THE STRUCTURE
    if ~isfield(rec,'S');rec.S=Sp;else rec.S=cat(5,rec.S,Sp);end    
    %%%THIS IS PROBLEMATIC, IT MAY NEED REFINEMENT
    if ~isfield(rec,'M');rec.M=W;else rec.M=cat(5,rec.M,W);end
end

indM=7:8;
if rec.Par.Mine.Modal==2;indM=cat(2,indM,12);end
for m=indM
    if ~any(ismember(rec.Dyn.Typ2Rec,m));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,m);end
    rec.Dyn.Typ2Wri(m)=1;
end
if ismember(12,indM)
    rec.x=rec.x.*rec.M;
    rec.x=removeOverencoding(rec.x,rec.Alg.OverDec);
    M=buildFilter(2*size(rec.x),'tukeyIso',0.5,gpu,1,1);
    %rec.x=abs(filtering(abs(rec.x),M,1)/20);
    rec.x=abs(filtering(abs(rec.x),M,1)/parS.softMaskFactor);
end


