
function rec=solveXTB_PT(rec)

tsta=tic;
clc

%SET DEFAULT PARAMETERS
rec=disorderAlgorithm(rec);
rec=pilotToneAlgorithm(rec);

%FILENAME
folderSnapshots=strcat(rec.Names.pathOu, filesep,'An-Ve_Sn');%YB: Anatomical - Volumetric
[~,fileName]=fileparts(generateNIIFileName(rec));

%LOGGING INFO
if rec.Dyn.Log
    logDir = strcat(rec.Names.pathOu,filesep,'An-Ve_Log', filesep);
    logName =  strcat( logDir, fileName,rec.Plan.Suff,rec.Plan.SuffOu, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~exist(logDir,'dir'); mkdir(logDir);end
    diary(logName) %eval( sprintf('diary %s ', logName))
end
fprintf('Running file: %s \n', strcat( fileName, rec.Plan.Suff, rec.Plan.SuffOu) )
c = clock; fprintf('Date of reconstruction: %d/%d/%d \n \n', c(1:3));

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.AcqVoxelSize;%Acquired voxel size
parXT=rec.Alg.parXT;parXB=rec.Alg.parXB;
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end 
gpu=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);%YB: set gpu to single(gpuDeviceCount && ~rec.Dyn.BlockGPU)

if size(rec.y,5)>=parXT.maximumDynamics;fprintf('Problem too big to fit in memory (%d dynamics for a limit of %d)\n',size(rec.y,5),parXT.maximumDynamics);rec.Fail=1;return;end
on=cell(1,5);for n=1:5;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
[resPyr,L,estT,resIso,estB]=pyramidPlan(voxSiz,parXT.resolMax,parXT.NWend,parXT.accel);% YB: resPy is used for resAni, which in its turn is used in downSampingOperators.m
resPyrMax = [ .25 .5  1  ];
resPyr=rec.Alg.resPyr;%resPyr = resPyrMax(end+1-min(length(rec.Alg.resPyr),length(resPyr)):end) ;
resIso=sqrt(2)*((prod(voxSiz).^(1/3))./resPyr);%resIso = [resIso(end-(length(resPyr)-1):end)]
L = length(resIso);
estT = parXT.PT.estT;%estT = ones(size(resIso));
estB = ones(size(resIso));

rec = makePTParamsCompatible(rec, resPyr, resIso, estT, estB);

fprintf('RECONSTRUCTION PLAN\n')
fprintf('   Levels:             %s \n',sprintf(' %d ',1:L))
fprintf('   Effective levels:   %s \n',sprintf(' %d ',(1:L) - sum(resPyr~=1)) )
fprintf('   Resolution:         %s \n',sprintf(' %g ',resPyr))
fprintf('   Motion estimation:  %s \n',sprintf(' %d ',estT))
fprintf('   B0 estimation:      %s \n',sprintf(' %d ',estB))
fprintf('   PT usage flag:      %s \n',sprintf(' %d ',parXT.PT.usePT))
fprintf('   PT subdividing:     %s \n\n',sprintf(' %d ',parXT.PT.subDiv))

%ROI COMPUTATION AND EXTRACTION IN THE READOUT DIRECTION
rec.Enc.ROI=computeROI(rec.M,[],[0 0 0],[1 0 0]); %YB: At this point non-logical values
if rec.Dyn.Debug>=2
    fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    fprintf('Acquisition directions: %s\n',rec.Par.Scan.MPS);
end
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if ~ismember(n,[5 12]);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1);end%YB: 1 at the end means 'forrward'. At the end also inverse operation
end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:3)=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
for n=typ2Rec'        
    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};%YB: 5 and 12 are noise and reconstruction. rec.x does not exist yet, so don't permute.
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end
if multDimSum(parXT.PT.usePT)>0 && isfield(rec,'PT') && ~isempty(rec.PT.pSliceImage)%PTHandling
    rec.PT.pSlice = permute(rec.PT.pSliceImage,perm);
    if gpuIn; rec.PT.pSlice = gpuArray(rec.PT.pSlice);end
elseif multDimSum(parXT.PT.usePT)>0
    parXT.PT.usePT=zerosL(parXT.PT.usePT);
    warning('No pilot tone signal detected in data. PT usage disabled.');
end
voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));
if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = [];end
rec.Par.Mine.permuteHist{end+1} = perm(1:4);

MS = voxSiz; 
[~,MT] = mapNIIGeom([], rec.Par.Mine.APhiRec,'permute',rec.Par.Mine.permuteHist{end} );%Don't change MS since already permuted in voxSiz
isValidGeom(MS,MT);

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
NX=size(rec.M);
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,rec.y);
NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    if parXT.perc(1)<1;fprintf('Number of compressed coil elements at%s%%: %d (%s )\n',sprintf(' %0.2f',parXT.perc*100),NS(4),sprintf(' %d',eivaS));else fprintf('Number of compressed coil elements: %d (%s )\n',NS(4),sprintf(' %d',eivaS));end
    fprintf('Number of coils used for intermediate image estimation: %d\n', eivaS(1));
    fprintf('Number of coils used for motion estimation: %d\n', eivaS(2));
    fprintf('Number of coils used for final image estimation: %d\n', eivaS(3));
end

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS %YB: apodisation means windowing
NY=size(y);NY(end+1:rec.Plan.NDims)=1;
y=bsxfun(@times,y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpuIn,parXT.apod))); %Apodize % YB: y in image domain
[y,NY]=resSub(y,5:rec.Plan.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
y=gather(y);

%ACQUISITION STRUCTURE
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
NProfs=numel(rec.Assign.z{2});
if NProfs<=1 || NProfs~=numel(rec.Assign.z{3});fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n');rec.Fail=1;return;end
if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med'))
    fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(4:5));
else; fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(7:8));
end
%YB: z{2} = PE1 and z{3} = PE2
kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{perm(n)}(:);end
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if mod(NShots,1)~=0;fprintf('NUMBER OF SHOTS: %.2f, PROBABLY AN INTERRUPTED SCAN\n',NShots);rec.Fail=1;return;end
NRepeats=NY(5);
if NShots<NRepeats
    NShots=NRepeats;
    NProfs=NShots*NEchos;
    kTraj=repmat(kTraj,[NRepeats 1]);
    NEchos=NProfs/NShots;
end
if rec.Dyn.Debug>=2
    fprintf('Number of echoes: %d\n',NEchos);
    fprintf('Number of shots: %d\n',NShots);
end
kTrajSS=reshape(kTraj,[NEchos NShots 2]);% YB: SS for snapshot I think
%PLOT TRAJECTORIES
if rec.Alg.WriteSnapshots;visTrajectory(kTrajSS,2,strcat(folderSnapshots,filesep,'Trajectories'),strcat(fileName,rec.Plan.SuffOu));end
kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.kRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kSize=(diff(kRange,1,2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
rec.Par.Mine.DISORDER.kTraj=gather(kTraj);
rec.Par.Mine.DISORDER.kRange=gather(kRange);
rec.Par.Mine.DISORDER.kShift=gather(kShift);
rec.Par.Mine.DISORDER.kSize=gather(kSize);
rec.Par.Mine.DISORDER.kIndex=gather(kIndex);

%STEADY-STATE CORRECTIONS
isSteadyState=0;
if mod(NShots,NRepeats)==0
    NShotsUnique=NShots/NRepeats;  
    isSteadyState=(NShotsUnique==1);%No shots have been identified
end
fprintf('Steady-state: %d\n',isSteadyState);
if isSteadyState
    kTrajs=abs(fct(kTraj)); % YB: By doin this frequency analysis, you don't actually need the sampling patter procided to the scanner!
    N=size(kTrajs);
    [~,iAs]=max(dynInd(kTrajs,1:floor(N(1)/2),1),[],1);    
    iM=min(iAs,[],2);
    NSamples=NProfs/(iM-1);
        if isfield(rec.Alg.parXT,'sampleToGroup');NSamples=rec.Alg.parXT.sampleToGroup;end
    NSamplesPerAcquisitionSweep=NSamples; %YB added to comapre sample grouping to original acquisition
    NSweepsAcquisition=round(NProfs/(NSamplesPerAcquisitionSweep));
    fprintf('Number of acquisition sweeps: %d\n',NSweepsAcquisition);
    fprintf('Temporal resolution per acquisition sweep: %.2f sec\n',NSamplesPerAcquisitionSweep*rec.Par.Labels.RepetitionTime(1)/1000);
    NSweeps=round(NProfs/(NSamples*parXT.groupSweeps));% YB: Groupsweeps only for the first estimation, as afterwards the grouping will be based on the motoincompression
    NSamples=NProfs/NSweeps; %YB: Not necessary a round number as profiles per segment might slightly deviate    
    fprintf('Number of acquisition sweeps to group: %d\n',parXT.groupSweeps);
    fprintf('Number of sweeps after grouping: %d\n',NSweeps);
    fprintf('Number of samples per sweep: %.1f\n',NSamples);   
    fprintf('Temporal resolution per sweep: %.2f sec\n',NSamples*rec.Par.Labels.RepetitionTime(1)/1000);  
else
    NSamples=NEchos;
    NSweeps=NShots;
end
NStates=NSweeps;
%SWEEP SUBDIVISION AND TIME INDEX
sweepSample=ceil(NSweeps*(((1:NProfs)-0.5)/NProfs));% YB: I think the 0.5 is just a small number (compared to NProfs) so that ceil returns the right numbers
stateSample=sweepSample;
timeIndex=zeros([NY(1:2) NY(5)]);
if gpu;timeIndex=gpuArray(timeIndex);end
for n=1:size(kIndex,1)
    for s=1:NY(5) % YB: over repeats, so if filled for one, break as they all contain same timing
        if timeIndex(kIndex(n,1),kIndex(n,2),s)==0
            timeIndex(kIndex(n,1),kIndex(n,2),s)=n;
            break;
        end
    end
end
for m=1:2;timeIndex=ifftshift(timeIndex,m);end%YB: This goes to downSamplingOperators.m which expects a sampling operator A in the non-shifted fft domain
if isSteadyState
    timeSample=(0:NProfs-1)*rec.Par.Labels.RepetitionTime(1)/1000;
else
    TPerShot=rec.Par.Labels.ScanDuration/NShots;    
    timeSample=(sweepSample-1)*TPerShot;
    if strcmp(rec.Par.Scan.Technique,'TSE') || strcmp(rec.Par.Scan.Technique,'TIR')
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*(rec.Par.Labels.TE(1)/(1000*NEchos/2));
    else
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*rec.Par.Labels.RepetitionTime(1)/1000;
    end
end
sweepToSample=cell(1,NSweeps);
sampleToSweep=1:NSweeps;
for n=1:NSweeps;sweepToSample{n}=n;end%YB: 

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
rec.x=zeros(NX,'single');%Uncorrected non-regularized
rec.d=zeros([NX L-sum(resPyr~=1)],'single');%Corrected non-regularized  %YB: For all the last (full resolution levels), the reconstructions are stored
if parXT.computeCSRecon;rec.r=zeros([NX L-sum(resPyr~=1)],'single');end %Corrected regularized
x=zeros(NX,'single');if gpuIn;x=gpuArray(x);end
typ2RecI=[12;16]; % 12-> Reconstruction 16-> Volumetric alignment
typ2RecI=[16]
warning('YB: disabled saving Acq NIFTI files since saved in SWCC and since it is the same for all the options for MoCo ...')
if parXT.computeCSRecon;typ2RecI=[typ2RecI;18];end %18-> Filtered reconstruction
if parXB.dephaseCorrection>0 || parXB.useTaylor;typ2RecF=29;else typ2RecF=[];end %29-> Linear coefficients of B0 fields with motion
for n=typ2RecI'
    if ~any(rec.Dyn.Typ2Rec==n);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end
end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZE TRANSFORM, NUMBER OF EXTERN ITERATIONS, TOLERANCES, IDENTIFY 
%STEADY STATE SEQUENCES 
rec.T=zeros([on{4} NSweeps 6],'single');

if ~isfield(parXB,'redFOV'); parXB.redFOV = parXT.redFOV;end
if ~isfield(parXB,'intravoxDephasing'); E.intraVoxDeph=0;else E.intraVoxDeph = parXB.intraVoxDeph;end
if  parXB.useSH>0
    SHNumCoef = size(SH_basis(2*on{3}, parXB.SHorder),2); 
    E.Db.c = zeros([on{4} NSweeps SHNumCoef],'like',real(x));E.Db.cr=E.Db.c;E.NMs = NSweeps;
    E.Db.useSH = 1;E.Db.SHorder = parXB.SHorder;
    E.Db.chist = [];
    E.Db.TE = rec.Par.Labels.TE * 1e-3;%in seconds
    if parXB.useSH==2; E.Dbs = E.Db; E = rmfield(E,'Db');end%define fields in scanner frame
    E = updateBasis(E,parXB, NX);
end
if  isfield(parXB,'useB1') && parXB.useB1
    E = updateBasisB1(E,parXB, NX);
    E.B1m.c = zeros([on{4} NSweeps size(E.B1m.B,2) ],'like',real(x));E.B1m.cr=E.B1m.c;E.NMs = NSweeps;
else
    parXB.useB1=0;
end
if isfield(parXT,'corrFact'); E.corrFact = parXT.corrFact;end

subT=ones(1,NStates);%To perform subdivisions of T
nExtern=99;
%nExtern=2;
tolType={'Energy','RelativeResidual2Error'};%For CG--without vs with motion correction
tol=[parXT.tolerSolve 0];%For CG
nIt=[300 1];%For CG
nowithin=0;%To stop estimation of within motion

%ENERGY
En=cell(1,L);Re=cell(1,L);EnX=[];
mSt=cell(1,L);iSt=cell(1,L);
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time arranging: %.3f s\n',tend);end

%SOLVE WITHOUT MOTION
tsta=tic;
resAni=single(on{3});
[~,indMinRes]=min(voxSiz);        
resAni(indMinRes)=voxSiz(indMinRes)/resPyr(end);
indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
resAni=voxSiz./resAni;
BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
if ~isfield(rec.Dyn,'GPUbS'); rec.Dyn.GPUbS = [2 4];end
E.bS=round([rec.Dyn.GPUbS].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
E.Je=1;%To use joint encoding/decoding
E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NY(n);E.Uf{n}.NX=NX(n);end%Folding
EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NY(n);EH.Ub{n}.NX=NX(n);end%Unfolding
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)
R=[];%Regularizer
ncx=1:eivaS(1);

%MASKING
CX.Ma=rec.M;
if rec.Alg.UseSoftMasking==0 && ~all(ismember(CX.Ma(:),[0 1])) % Need to convert rec.M into discrete mask
    %maxint = max(abs(rec.M(:)));  level = graythresh(abs(rec.M ) / maxint);
    maxint =1;  level = 0.45 * mean(abs(rec.M(:)));
    CX.Ma = single( abs(rec.M)> level*maxint);
    %CX.Ma = imopen(CX.Ma, strel('sphere', 3));
    CX.Ma = imclose(CX.Ma, strel('sphere', 15));
    CX.Ma = imdilate(CX.Ma, strel('sphere', 5));
    rec.M = CX.Ma;
end
plotND([], RSOS(multDimMea(rec.y,5:16)),[0 .51e6],[],0,[],MT,[],rec.M ,{2},2);title('Mask used for reconstruction')
if ~exist(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask'),'dir'); mkdir(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask'));end
%saveFig(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask',rec.Names.Name)) ;

if gpuIn;CX.Ma=gpuArray(CX.Ma);end
if rec.Alg.UseSoftMasking==2;CX.Ma(:)=1;end %YB: 2 means not masking
discMa=all(ismember(CX.Ma(:),[0 1])); %YB: discrete masking - if every element is rather 0 or 1
if ~discMa %Soft masking
    fprintf('Using soft masking\n');
    M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpuIn,1,1);
    CX.Ma=abs(filtering(CX.Ma,M,1));M=[];
    if size(S,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction
        NMa=[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX];
        for m=1:3
            if ~isempty(E.Uf{m});CX.Ma=fold(CX.Ma,n,E.Uf{m}.NX,E.Uf{m}.NY,E.UAf{m});end
        end         
        CX.Ma=resampling(CX.Ma,NMa,2);
    end
    if isfield(rec,'W');CX.Ma=bsxfun(@times,CX.Ma,rec.W);end
    Ti.la=CX.Ma;
    R.Ti.la=1./(abs(Ti.la).^2+0.01);%YB: This is an example where the regularisation is precomputed as R'R and not specifically computed as R'R in regularize.m  
    CX=[];
end

%SENSITIVITIES
E.Sf=dynInd(S,ncx,4);E.dS=[1 ncx(end)];
if gpuIn;E.Sf=gpuArray(E.Sf);end

%%% B0 Shim
if ~isfield(parXB, 'modelShim');parXB.modelShim=0;end
if ~parXT.exploreMemory && isfield(rec.Par.Labels,'Shim') && ( isfield(parXB, 'modelShim') && parXB.modelShim==1 )% Also computed when est(B)==0
    E = updateShim(E, rec.Par.Labels.Shim, NX, rec.Par.Mine.APhiRec*Tpermute(rec.Par.Mine.permuteHist{end}));
    E.Shim.TE = rec.Par.Labels.TE * 1e-3;%in seconds
    plotND({abs(x),1,0}, E.Shim.B0,[-200 200],[],1,{[],2}, rec.Par.Mine.APhiRec*Tpermute(rec.Par.Mine.permuteHist{end}),{'Shim B0 field'},[],[],9);
end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

yX=mean(dynInd(y,ncx,4),5);
if gpuIn;yX=gpuArray(yX);end
gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
if gibbsRing~=0 && gibbsRing<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpuIn,gibbsRing));end
nX=nIt(1);
if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

if nX~=0
    [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,x,nX,tolType{1},tol(1));
    rec.x=gather(xRes);
end
if ~isempty(EnX)%Print and store final energy
    EnX=EnX/numel(yX);
    if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
    EnX=[];
end;yX=[];

if gpuIn;y=gpuArray(y);end
Residuals=encode(rec.x,E)-y;
for m=1:2;Residuals=fftGPU(Residuals,m)/sqrt(NY(m));end
Residuals=normm(Residuals,[],3:4);
if rec.Alg.WriteSnapshots
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NX(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NX(3)):NX(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);end
    visSegment(flip(permute(dynInd(rec.x,E.nF,3),[2 1 3]),2),[],2,1,[],[],'Uncorrected',strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,'_Aq'));
    visResiduals([],[],[],Residuals,kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',0)));
end
E.Sf=[];

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing non-corrected reconstruction: %.3f s\n',tend);end

x=resampling(xRes,NX); %YB: Start first level with the uncorrected guess

%%% PILOT TONE
if multDimSum(parXT.PT.usePT)>0 %PTHandling
    %Initialise variables
    rec.PT.NChaPT= size(rec.y,4);% # Non-compressed channels
    ssPTConversion=[]; ssPTConversion.Geom.APhiRecOrig = rec.Par.Mine.APhiRecOrig; ssPTConversion.Geom.permuteHist=rec.Par.Mine.permuteHist;%For transformation conversion
    if parXT.PT.useRASMotion==2;ssPTConversion.Geom.useLog=1;end
    %Shift to be consistent with sampling
    for n=1:2; rec.PT.pSlice=fftshiftGPU(rec.PT.pSlice,n);end
    %Create calibration matrix
    rec.PT.Af = zeros([rec.PT.NChaPT, 6+parXT.PT.calibrationOffset],'like',x);%6 motion parameters and possible offset
    rec.PT.Ab = zeros([6, rec.PT.NChaPT+parXT.PT.calibrationOffset],'like',x);%6 motion parameters and 1 offset
    %MB handling
    if rec.Alg.parXT.PT.combineMB==1
        rec.PT.pSlice = dynInd(rec.PT.pSlice, rec.PT.idxMB,3); %Extract peak in RO direction
    elseif rec.Alg.parXT.PT.combineMB==2
        rec.PT.pSlice = multDimMea(rec.PT.pSlice, 3);%Average frequencies
    end
    %Normalise
    if parXT.PT.normaliseCoils; for pp=1:NChaPT; rec.PT.pSlice = dynInd(rec.PT.pSlice,pp,4, dynInd(rec.PT.pSlice,pp,4)/sqrt(normm(dynInd(rec.PT.pSlice,pp,4))));end;end
    %Filter
    %rec.PT.pSlice = filterPT (rec.PT.pSlice , parXT.PT,NY, kIndex, 'pSlice');
    plotPTSlice(rec.PT.pSlice, 182, sprintf('Original PT signal filtersed using median filter with width %.1f ms', parXT.PT.medFiltKernelWidthms) , NY, kIndex);
end

if isfield(rec.Par.Labels ,'H0') && rec.Par.Labels.H0 ==3;MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];else MTT = eye(4);end
rec.Par.Mine.APhiRec = MTT*rec.Par.Mine.APhiRec;
MT = MTT * MT;

%ITERATING THROUGH THE DIFFERENT LEVELS
groupSweepsRecon = parXT.groupSweeps;
for l=1:L    
    %WE RESET MOTION CONVERGENCE INSIDE
    nReset=0;THist=[]; % YB: setting nReset to zero ensures resetting motion convergence
    leff=l-sum(resPyr~=1);

    %UPSAMPLE THE TRANSFORM AND GENERATE MOTION PARAMETERS   
    timesToSubDiv = max(1,log2(max(subT(:)))); %SubT is set as a power of 2 
    subT = min(subT,2);%YB: if subT has value 4, you want to set it to 2 and timesToSubDiv will take care of the multiple subdivision
    for ii = 1:timesToSubDiv
        rec.T=repelem(rec.T,1,1,1,1,subT,1); %YB: so intra-shot motion is not an extra dimension but are repeated elements in the 5th dimension!
        if parXB.useSH==1; E.Db.c = repelem(E.Db.c,1,1,1,1,subT,1);elseif parXB.useSH==2; E.Dbs.c = repelem(E.Dbs.c,1,1,1,1,subT,1);end
        subTT=repelem(subT,1,subT);
        sampleToSweep=repelem(sampleToSweep,1,subT);
        
        subdivided = subT(1)==2; 
        if subdivided; groupSweepsRecon=groupSweepsRecon/2;end
        if subdivided && groupSweepsRecon >= 1
            NSweeps=NSweeps*2;%Increase but don't let it exceed NSweepsAcquisition
        end
        
        sweepToSample=cell(1,NSweeps);
        c=0;
        for n=1:length(subT)
            if subT(n)==2
                indState=find(stateSample==(n+c));
                Ni=length(indState);
                indState=indState(floor(Ni/2)+1); % YB: halfway since you split the segment in half
                stateSample(indState:end)=stateSample(indState:end)+1; % YB: so max(stateSample) can change over resolutions
                sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c:n+c+1];%YB: Not exactly similar to sweepSample (not used anymore) and relates the states to the sweeps(stays constant)
                c=c+1;
            else
                sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c];
            end
        end
        subT = subTT;%YB new for PT
    end

    timeState=regionprops(stateSample,timeSample,'MeanIntensity'); %YB: if find(subT(:)==1) is not empty, max( stateSample) > NSweeps
    timeStateLimitInf=regionprops(stateSample,timeSample,'MinIntensity');
    timeStateLimitSup=regionprops(stateSample,timeSample,'MaxIntensity');
    timeState={timeState.MeanIntensity};
    timeStateLimitInf={timeStateLimitInf.MinIntensity};
    timeStateLimitSup={timeStateLimitSup.MaxIntensity};
    timeState=cat(1,timeState{:});
    timeStateLimitInf=cat(1,timeStateLimitInf{:});
    timeStateLimitSup=cat(1,timeStateLimitSup{:});
    %if leff<=1;timeSweep=timeState;end
    if groupSweepsRecon >= 1%YB: added for when you subdivide sooner than leff>0 %WEISSUE
        timeSweep=timeState;
    end
    timeStateoutlWePT = timeState;
    
    %NOTE THIS ONLY WORKS IF THERE IS NO INTERSHOT MOTION
    if leff<=1
        if ~isSteadyState;timeState=cat(3,timeStateLimitInf,timeStateLimitSup);else timeState=cat(3,timeStateLimitInf,timeState,timeStateLimitSup);end
    end
    
    %RELATIVE SPATIAL RESOLUTION CALCULATIONS AT THIS LEVEL
    resAni=single(on{3});
    [~,indMinRes]=min(voxSiz);        
    resAni(indMinRes)=voxSiz(indMinRes)/resPyr(l);% If no rounding issues, this should be the 4mm of the coarses resolution
    indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
    resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));  %YB: Taking is maximum just says that if you have a resolution >4m you should increase dimension of array! see downsampleOperators.m
    resAni=voxSiz./resAni;%YB: goal of this is when you take -log2(resAni) that you obtain the 2- subdivisions per dimension so that each has ~=4mm in the coarsest scale and is isotropic 

    NT=size(rec.T);    
    NStates=NT(5);%YB: NStates changes over resolutions but NSweeps does not!
    subT=ones(1,NStates);%To perform subdivisions of T 
    %YB: this is reset in a way that for the new NStates, all set to 1
    cT=zeros(NT(1:5),'single');%Flag to indicate whether a given motion state can be assumed to have converged     
    w=parXT.winit*ones(NT(1:5),'single');%Weights of LM iteration            
    NSamplesPerState=accumarray(stateSample(:),ones(NProfs,1))';
    voxSizRes = voxSiz .* ( NX ./  ceil(NX(1:3)./(2.^(-log2(resAni)))) );
    if rec.Dyn.Debug>=2
        fprintf('\n========================== Resolution level: %d ==========================\n',l);
        fprintf('Voxel size: %.2fx%.2fx%.2f mm\n',voxSizRes);
        fprintf('Number of motion states: %d\n',NStates);
        %fprintf('Minimum/Maximum number of echoes: %d/%d\n',min(NSamplesPerState),max(NSamplesPerState));
        if any(NSamplesPerState(:)<parXT.echoesMin);fprintf('Trying to operate below the minimum required number of reads, motion estimates will not be obtained for stability\n');end
    end
    if estT(l)
        cT(NSamplesPerState(:)<parXT.echoesMin)=1;
        cT=reshape(cT,NT(1:5));
    end

    %ENCODING STRUCTURES AND BLOCK SIZES FOR GPU COMPUTATION AT THIS LEVEL
    fprintf('Downnsampling k-space and image data.\n');
    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gpuArray(S);rec.M=gpuArray(rec.M);end
    if ~discMa%Soft masking
        [xRes,SRes,yRes,timeIndexRes,R.Ti.la,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,Ti.la,gibbsRing,0,gpuIn); %YB: 0 confirms that y is in image domain
        R.Ti.la=1./(abs(R.Ti.la).^2+0.01);
    else
        [xRes,SRes,yRes,timeIndexRes,CX.Ma,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,rec.M,gibbsRing,0,gpuIn);%YB:attention: discrete mask provided by CX.Ma not used! 
    end%YB: These E.Fof and E.Fob are the forward and backward fourier operators used for the transformations! Not for the actual fft per segment (carried out by HarmonicSampling)

    NXres=size(xRes);NXres(end+1:3)=1;NYres=size(yRes);NYres(end+1:5)=1;
    if multDimSum(parXT.PT.usePT)>0;rec.PT.pSliceRes = resampling(rec.PT.pSlice, NYres(1:2));end%PTHandling -  %Resampling in the image domain
    
    [MSRes,MTRes] = mapNIIGeom(MS, MT, 'resampling',[], NX, NXres);

    %B1 FIELD INITIALISATION
    if ~parXT.exploreMemory && parXB.useB1
        E = updateBasisB1(E, parXB, NXres);
        if strcmp(rec.Par.Labels.FatShiftDir,'F');E.B1m.Optim.Mbg=1:floor((1-parXB.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.B1m.Optim.Mbg=floor(1+parXB.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.B1m.Optim.Mbg=1:NXres(3);end
    end
    
    %B0 FIELD INITIALISATION
    %%% Basis functions in head frame
    if ~parXT.exploreMemory && parXB.useSH==1 % Also computed when est(B)==0
        E = updateBasis(E, parXB, NXres);
        if strcmp(rec.Par.Labels.FatShiftDir,'F');E.Db.Optim.Mbg=1:floor((1-parXB.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.Db.Optim.Mbg=floor(1+parXB.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Db.Optim.Mbg=1:NXres(3);end
        if ~isfield(E.Db.Optim,'deTaylor'); E.Db.Optim.deTaylor = parXB.deTaylor;end
    end 
    %%% Basis functions in scanner frame
    if ~parXT.exploreMemory && parXB.useSH==2 % Also computed when est(B)==0
        E = updateBasis(E, parXB, NXres);
        if strcmp(rec.Par.Labels.FatShiftDir,'F');E.Dbs.Optim.Mbg=1:floor((1-parXB.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.Dbs.Optim.Mbg=floor(1+parXB.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Dbs.Optim.Mbg=1:NXres(3);end
        if ~isfield(E.Dbs.Optim,'deTaylor'); E.Dbs.Optim.deTaylor = parXB.deTaylor;end
    end 
    %%% B0 Shim
    if ~isfield(parXB, 'modelShim');parXB.modelShim=1;end
    if ~parXT.exploreMemory && isfield(rec.Par.Labels,'Shim') && ( isfield(parXB, 'modelShim') && parXB.modelShim==1 )% Also computed when est(B)==0
        E = updateShim(E, rec.Par.Labels.Shim, NXres, MTRes);
        E.Shim.TE = rec.Par.Labels.TE * 1e-3;%in seconds
    end     
    if parXB.modelShim
        %Move field from object to shim
        z = CNCGUnwrapping(xRes, MSRes ,'MagnitudeGradient4','LSIt'); %Ideally has a magnitude as well
        z = z/2/pi/E.Shim.TE; % in Hz
        B = SH_basis( size(xRes) , 6);%First order SH
        zInt{l} = interpolateBasis(z, B);
        for j = 1:l; E.Shim.B0 = E.Shim.B0 + resampling( zInt{j}, NXres); end
        xRes = xRes .* conj(exp(+1i*2*pi*E.Shim.TE*zInt{l}));
    end
    
    %%% Taylor model
    if ~isfield(rec.Par.Mine,'APhiRecOrig'); rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;end
    if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = []; end
    
    if ~parXT.exploreMemory && any(estB(1:l-1)) && parXB.dephaseCorrection>0;E.Dc.D=resampling(DB,NXres,0,2*ones(1,3));DB=[];end % YB: DB are the atoms in the original dimensions
    if ~parXT.exploreMemory && any(estB) && parXB.useTaylor
        E = updateTaylor (E, xRes, NXres, voxSiz .* ( NX ./ ceil(NX(1:3)./(2.^(-log2(resAni)))) ), parXB, rec.Par.Labels);
        E.Dl.Geom.APhiRecOrig = rec.Par.Mine.APhiRecOrig;
        E.Dl.Geom.permuteHist = rec.Par.Mine.permuteHist;
        E.Dl.Geom.levelScale = ones(size(NX./NXres));%disabled
    end

    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gather(S);rec.M=gather(rec.M);end
    if l==L;S=[];y=[];end
    for n=1:2;yRes=fftGPU(yRes,n)/sqrt(NYres(n));end %YB: Now the measure data is in k-space and sqrt makes it the unitary transform (could have been included in fftGPU)
    %k-space dividing by sqrt(NYres) and PT by NY(res) because buildHarmonicSampling uses the normalised fft and this is internally consistent. This is not true for PT (where we use the fft, so you need to take out the NYres
    if multDimSum( parXT.PT.usePT) > 0
       for n=1:2
           %rec.PT.pSliceRes=fftshiftGPU(rec.PT.pSliceRes,n);
           %rec.PT.pSliceRes=fftGPU(rec.PT.pSliceRes,n)/sqrt(NYres(n));
           rec.PT.pSliceRes=fftGPU(rec.PT.pSliceRes,n)/NYres(n);
           rec.PT.pSliceRes=rec.PT.pSliceRes*sqrt(NY(n));%This is to make signal the same levels as pTimeTest (does not matter since between levels this does not matter)
       end 
    end
    
    %WEIGHTING FUNCTION FOR MOTION ESTIMATION (ROI)
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NXres(3);end
    EH.Mbe=zeros([1 1 NYres(3)],'like',real(xRes));EH.Mbe(E.nF)=1;%To use the extracted region for energy computation
    if rec.Dyn.Debug>=2; plotND([], dynInd(resampling(abs(rec.x),NXres), E.nF,3),[0 .5*1e6],[],0,{[],2},MTRes,[],[],[],3);end %YB to see extracted region

    %HARMONIC SAMPLING
    E.mSt=timeIndexRes;  
    E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0)); % YB: Now mSt a 2D grid with the same states having the value of the state number. Also note that timeIndexRes contains indices of the samples, so biggere than numel(timeIndexRes)
    
    if estT(l) || ~isfield(EH,'We');EH.We=[];end    
    if parXT.fractionOrder~=0;GRes=buildFilter([NXres(1:2) length(E.nF)],'FractionalFiniteDiscreteIso',NX./NXres,gpuIn,parXT.fractionOrder);else GRes=[];end%For fractional finite difference motion estimation           
    fprintf('Building harminic sampling structures.\n');
    [ySt,FSt,GRes,E.iSt,E.nSt,E.mSt,nSa]=buildHarmonicSampling(yRes,GRes,E.mSt,NStates,NXres(1:2));iSt{l}=cat(1,E.iSt{:});yRes=gather(yRes);%Samples as cell arrays
    if any(nSa(1:NStates)==0);nowithin=1;fprintf('Probably not DISORDER, returning without performing corrections\n');return;end   
    fprintf('Number of acquisition sweeps: %d\n',NSweepsAcquisition);
    fprintf('Temporal resolution of acquisition sweep: %.2f sec\n',NSamplesPerAcquisitionSweep*rec.Par.Labels.RepetitionTime(1)/1000);
    fprintf('Number of reconstruction groups: %d\n',NStates);
    fprintf('Minimum/Maximum number of echoes per reconstruction grouping at full resolution: %d/%d\n',min(NSamplesPerState),max(NSamplesPerState));
    fprintf('Minimum/Maximum number of echoes per reconstruction grouping at current resolution: %d/%d\n',([min(nSa),max(nSa)]));
    fprintf('Minimum/Maximum temporal resolution reconstruction grouping at current resolution: %.2f/%.2f sec\n',([min(nSa),max(nSa)])*rec.Par.Labels.RepetitionTime(1)/1000);
    
    %PILOT TONE GROUPING 
     if multDimSum(parXT.PT.usePT)>0 %PTHandling
        %Data type handling - only do after all the fft's/fftshift's are applied
        if parXT.PT.useMagn==1;rec.PT.pSliceRes = abs(rec.PT.pSliceRes); parXT.PT.isComplex = 0;
        elseif parXT.PT.usePhase==1;rec.PT.pSliceRes = angle(rec.PT.pSliceRes); parXT.PT.isComplex = 0;
        else; parXT.PT.isComplex = 1;
        end
        if ~parXT.PT.isComplex; [rec.PT.Ab, rec.PT.Af]=parUnaFun({rec.PT.Ab,rec.PT.Af},@(x) real(x));end
        
        rec.PT.pTimeRes = permute(rec.PT.pSliceRes,[1 2 5 3 4]);%5 is for the multiple repeats/shots
        rec.PT.pTimeRes = reshape(rec.PT.pTimeRes,[prod(NYres([1:2 5])) 1 1 rec.PT.NChaPT]);%rec.PT.pTimeRes = resPop(rec.PT.pTimeRes ,1:3,[],1);

        rec.PT.pTimeRes = dynInd( rec.PT.pTimeRes,dynInd(iSt{l},1:length(find(E.mSt<=NStates)),1),1);% Within each grouping, samples are not sorted in temporal order, so plot might look weird. Does not matter since you group all those samples (mean)
        rec.PT.pTimeRes = permute( rec.PT.pTimeRes , [4 1 2:3]);%Make PT signal as NCha x NStates
        if parXT.PT.usePTWeightsFlag; [rec.PT.WePT, rec.PT.outlWePT]= getPTWeighting(rec.PT.pTimeRes , E.mSt, NStates);else;rec.PT.WePT=ones([1 NStates]);rec.PT.outlWePT=zerosL(rec.PT.WePT);end
            visResiduals(rec.PT.WePT,rec.PT.outlWePT,timeStateoutlWePT,[],[],0, [],[],899);
            plotPTSignals(rec.PT.pTimeRes,179,sprintf('Temporal PT signal at resolution level %d: Resolution %s - Relative grouping %.1f',l, sprintf('%.2fmm ',MSRes), groupSweepsRecon),rec.PT.outlWePT,permute( E.mSt(E.mSt<=NStates),[2 1]));
        rec.PT.pTimeResGr = groupPT( (rec.PT.pTimeRes), permute( E.mSt(E.mSt<=NStates),[2 1]) );
            plotPTSignals(rec.PT.pTimeResGr,180,sprintf('Grouped temporal PT signal at resolution level %d (%s)',l, sprintf('%.2fmm ',MSRes)),rec.PT.outlWePT);              
    end
    assert((NSweeps/NSweepsAcquisition)<=1.25,'NSweeps cannot be smaller than the number of acquired sweeps in order to make the weights to be calculated (see DISORDER paper).')
    
    %if leff<=1;mStSweeps=E.mSt;end%YB: TODO: change because this was originally rigt because no subdividing before leff 2 - now we do with parXT.subDiv
    if groupSweepsRecon >= 1
        mStSweeps=E.mSt;
        stateSampleSweep=stateSample;
        NSweepsForWe=NSweeps;
    else %Need to recalculate buildHarmonicSampling.m with right number of sweeps
        mStTT=timeIndexRes;  
        mStTT(timeIndexRes~=0)=stateSampleSweep(timeIndexRes(timeIndexRes~=0)); 
        NSweeps=NSweepsForWe;
        [~,~,~,~,~,mStSweeps]=buildHarmonicSampling(yRes,GRes,mStTT,NSweepsForWe,NXres(1:2));        
    end

    BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
    E.bS=round(rec.Dyn.GPUbS.^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
    if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
    E.Je=1;%To use joint encoding/decoding
    E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NYres(n);E.Uf{n}.NX=NXres(n);end%Folding
    EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NYres(n);EH.Ub{n}.NX=NXres(n);end%Unfolding
    [E.UAf,EH.UAb]=buildFoldM(NXres(1:3),NYres(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)

    %CONTROL OF TRANSFORM CONVERGENCE / CONSTRAIN THE TRANSFORM / ROI IN
    %READOUT FOR MOTION ESTIMATION
    E.w=w;E.cT=cT;E.winit=parXT.winit;
    if parXT.meanT;CT.mV=5;else CT=[];end
    E.tL=parXT.traLimXT*resIso(l);%Minimum update to keep on doing motion estimation

    cont=0; 
    if isfield(E,'En');E=rmfield(E,'En');end    

    %ITERATIVE SOLVER
    for n=1:rec.Alg.nExtern(l)
        if rec.Dyn.Debug>=2;fprintf('--------------Iteration: %d ---------------\n',n);end
        %RESET CONVERGENCE OF MOTION
        if mod(cont,nReset)==0 || all(E.cT(:)) || parXT.convTransformJoint
            nReset=nReset+1; cont=0;%We increment the number of iterations till next reset
            E.cT(:)=0;%We reset
            E.cT(cT==1)=1;%We fix to 1 the low frequency states in case that for whatever reason cT was set to 1
            if rec.Dyn.Debug>=2;fprintf('Resetting motion convergence\n');end
        end
        if (n==1 || l<L) && any(subTT(:)==2)   %YB: Remember, subTT has same dimensions as NStates whereas subT as NSweeps        
            E.cT(subTT==1)=1;
        end%For within-shot in the first iteration we only estimate those coming from outliers
        if rec.Dyn.Debug>=2;fprintf('Explored motion states: %d of %d\n',NStates-sum(single(E.cT)),NStates);end
        cont=cont+1;

        %TRANSFORM PARAMETERS PER SEGMENT/OUTLIERS/DATA/HARMONICS/SEGMENT INDEXES
        if isfield(parXT,'disableGrouping'); disableGrouping = parXT.disableGrouping;else disableGrouping = 0;end
        if parXB.useSH==1
            [E.Tr,E.Fs,E.nEc, E.Db.cr,E.Db.c_idx ]=compressMotion(rec.T,FSt, resIso(end), parXT, E.Db.c, disableGrouping);%YB: E.nEc are the starting indiced - used in the reconstruction (decode.m)
        elseif parXB.useSH==2
            [E.Tr,E.Fs,E.nEc, E.Dbs.cr,E.Dbs.c_idx ]=compressMotion(rec.T,FSt, resIso(end), parXT, E.Dbs.c, disableGrouping);%YB: E.nEc are the starting indiced - used in the reconstruction (decode.m)
        else
            [E.Tr,E.Fs,E.nEc]=compressMotion(rec.T,FSt,resIso(end),parXT,[],disableGrouping);
        end
        %if parXB.useB1;E.B1m.cr = E.B1m.c;end%CHANGE!
        if rec.Dyn.Debug>=2;fprintf('Number of binned motion states: %d of %d\n',size(E.Tr,5),NT(5));end
        if isfield(E,'Dl') || isfield(E,'Ds');E = transformationConversion(E);end%Make motion parameters for linear model 
        
        %SEGMENTS / MOTION STATES / BLOCK SIZES / TRANSFORM FACTORS 
        if ~isempty(E.Fs);E.NSe=length(E.Fs{1});E.NMs=size(E.Tr,5);else E.NSe=1;E.NMs=1;end%Number of segments (including out of elliptic shutter) / Number of motion states            
        E.NRe=E.NSe-E.NMs; %Number of outern segments
        E.dS(1)=E.NSe;

        %CG SOLVER
        ncx=1:eivaS(1+2*estT(l));
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn;E.Sf=gpuArray(E.Sf);end

        %PRECONDITION
        if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
        else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
        end            

        yX=dynInd(ySt,ncx,4);            
        %LAST ITERATION OUTLIER REJECTION
        EH.We=[]; %YB added - MAYBE CHANGE!
        %Weights based on residuals
        if  (  (leff>0 && estT(l)==0) || (isfield(parXB,'weightedEst') && parXB.weightedEst)  ) ...
            && rec.Alg.AlignedRec>=2 && exist('We','var') %YB: use weihts during all levels + for GDsolverLinear 
            EH.We=ones([size(yX,1) 1],'like',Residuals);
            if length(We)==NSweeps
                for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end %Soft weights    %YB: weights only used for leff>0 and est(T)=0   (this means at full resolution and when at the end of a resolution or when a the level where estT=0       
                visResiduals(We,outlWe,timeSweep,[],[],0, [],[],899);
                fprintf('Robust reconstruction: setting weights from k-space outlier detection.\n');
            else
                warning('Robust reconstruction with weights from k-space outlier detection DISABLED.\n Weights size is not the same as the number of Sweeps.\n Check the rec.parXT.PT.subDiv.');%
            end
        end
        %Weights based on PT signal - priority on the We from residuals and over-written
        if ((leff>0 && estT(l)==0) ||  parXT.PT.usePTWeightsFlag==2) && parXT.PT.usePTWeightsFlag  && rec.Alg.AlignedRec>=2 && isfield(rec.PT,parXT.PT.usePTWeightsMethod) 
            EH.We=ones([size(yX,1) 1],'like',Residuals);
            if strcmp(parXT.PT.usePTWeightsMethod, 'WePT')
                for s=1:max(E.mSt(E.mSt<=NStates));EH.We(E.mSt==s)=rec.PT.WePT(s);end %YB: Now use E.mSt instead of mStSweeps because PT weights defined for finest grouping
            elseif strcmp(parXT.PT.usePTWeightsMethod, 'outlWePT')
               for s=1:max(E.mSt(E.mSt<=NStates));EH.We(E.mSt==s)=1-single(rec.PT.outlWePT(s));end %YB: Now use E.mSt instead of mStSweeps because PT weights defined for finest grouping
            end
            visResiduals(rec.PT.WePT,rec.PT.outlWePT,timeStateoutlWePT,[],[],0, [],[],899);
            fprintf('Robust reconstruction: setting %s from PT outlier detection.\n',parXT.PT.usePTWeightsMethod);
        end
        %EH.We=ones([size(yX,1) 1],'like',Residuals);
        %warning('YB: weighting disabled')
        if estT(l)
            if l==1 && n==1;nX=nIt(2);
            elseif l==1;nX=max(nX-1,1);
            elseif resPyr(l)~=resPyr(l-1) && n==1;nX=max(nIt(2)-1,1);
            elseif resPyr(l)~=resPyr(l-1);nX=max(nX-1,1);
            elseif resPyr(l)==resPyr(l-1) && n==1;nX=0;
            elseif resPyr(l)==resPyr(l-1) && n==2;nX=max(nIt(2)-2,1);
            elseif resPyr(l)==resPyr(l-1);nX=max(nX-(n-2),1);
            end
        else
            nX=nIt(1);%YB: You also arrive here when convergence of previous is set to zero, so estT(l) is set to 0 there and you run one full CG of this level, set subT=2 and go to next level
        end
        if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

        %%% PIlot tone motion prediction %PTHandling
        if multDimSum(parXT.PT.usePT)>0 && parXT.PT.usePT(l)>0 && ~(parXT.PT.usePT(l)==1 && n<parXT.PT.nItActivate(l))
            TrPTForw = gather( pilotTonePrediction (rec.PT.Af , rec.PT.pTimeResGr, 'PT2T','forward', rec.PT.NChaPT, size(rec.T,6)) );
            TrPTBack = gather( pilotTonePrediction (rec.PT.Ab , rec.PT.pTimeResGr, 'PT2T','backward', rec.PT.NChaPT, size(rec.T,6)));
            if strcmp(parXT.PT.calibrationModel, 'forward'); E.Tr = TrPTForw;else; E.Tr=TrPTBack;end
            if parXT.PT.useRASMotion; [~,E.Tr] = transformationConversion([],E.Tr,ssPTConversion,[],[], 0, size(x));end %Inverse conversion: From RAS->IJK
            E = transformationConversion(E);%Update motion parameters in E.Dl
            rec.T=E.Tr; 
            fprintf('<strong>Pilot Tone</strong>: Using PT signal to predict motion states using %s model.\n', parXT.PT.calibrationModel);
            if parXT.PT.usePT(l)==2; estT(l)=0;fprintf('<strong>Pilot Tone</strong>: usePT(%d)==2, so setting estT(%d)=0.\n',l,l);end
        end
        
        if nX~=0 %YB: solving for X and B0
            %%% BASIS FUNCTIONS in head frame
            if estT(l) && estB(l) && parXB.useSH==1 && ( l>1 || n >= parXB.Optim.basisDelay)
                if ~isfield( parXB, 'nGD_B' );  parXB.nGD_B = 4;end
                if E.Db.Optim.deTaylor(1)>1 && length(E.Db.Optim.deTaylor)>1 && (l>1||n>E.Db.Optim.deTaylor(2)); E.Db.Optim.deTaylor=1;end                
                fprintf('** B0 estimation - basis functions in head reference- %d iterations \n', parXB.nGD_B+ 0*single(n==1) );
                EnBefore = computeEnergy(yX,xRes,E,R,EH);
                [E, xRes] = GDsolver(yX,E,EH,P,CX,R,xRes,parXB.nGD_B+ 0*single(n==1),tolType{estT(l)+1},tol(estT(l)+1) , rec.Dyn.Debug);
                EnAfter = computeEnergy(yX,xRes,E,R,EH);
                if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end
                %visBasisCoef(E.Db);
            end
            %%% BASIS FUNCTIONS in scanner reference
            if estT(l) && estB(l) && parXB.useSH==2 && ( l>1 || n >= parXB.Optim.basisDelay)
                if ~isfield( parXB, 'nGD_B' );  parXB.nGD_B = 4;end%if E.Dbs.Optim.deTaylor(1)>1 && length(E.Dbs.Optim.deTaylor)>1 && (l>1||n>E.Dbs.Optim.deTaylor(2)); E.Dbs.Optim.deTaylor=1;end
                fprintf('** B0 estimation - basis functions in scanner reference- %d iterations \n', parXB.nGD_B+ 0*single(n==1) );
                EnBefore = computeEnergy(yX,xRes,E,R,EH);
                [E, xRes] = GDsolverScannerRef(yX,E,EH,P,CX,R,xRes,parXB.nGD_B+ 0*single(n==1),tolType{estT(l)+1},tol(estT(l)+1) , rec.Dyn.Debug);
                EnAfter = computeEnergy(yX,xRes,E,R,EH);
                if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end
                plotND({abs(xRes),1,0},dephaseBasis(E.Dbs.B,dynInd(E.Dbs.cr,1,5),NXres, [], 1) ,[],[],1,{[],2},rec.Par.Mine.APhiRec*Tpermute(rec.Par.Mine.permuteHist{end}),[],[],[],345 );
            end
            %%% B1 MODULATION
            if estT(l) && estB(l) && parXB.useB1 && ( l>1 || n >= parXB.Optim.B1Delay)
                if ~isfield( parXB, 'nGD_B' );  parXB.nGD_B = 4;end                
                fprintf('** B1 estimation - %d iterations \n', parXB.nGD_B+ 0*single(n==1) );
                EnBefore = computeEnergy(yX,xRes,E,R,EH);
                [E, xRes] = GDsolverB1(yX,E,EH,P,CX,R,xRes,parXB.nGD_B+ 0*single(n==1),tolType{estT(l)+1},tol(estT(l)+1) , rec.Dyn.Debug);
                EnAfter = computeEnergy(yX,xRes,E,R,EH);
                if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end
                plotND(abs(xRes),dephaseBasis(E.B1m.B,dynInd(E.B1m.cr,1,5),NXres, [], 1) ,[],[],1,{[],2},rec.Par.Mine.APhiRec*Tpermute(rec.Par.Mine.permuteHist{end}),[],[],[],346 );
            end
            %%% TAYLOR MODEL
            if estT(l) && estB(l) && parXB.useTaylor && any(E.Tr(:)~=0) && ( l>1 || n >= parXB.Optim.taylorDelay) 
                if ~isfield(parXB, 'nGD_Taylor');  parXB.nGD_Taylor = 4;end
                fprintf('** B0 estimation - Taylor voxel model - %d iterations \n', parXB.nGD_Taylor + 0*single(n==parXB.Optim.taylorDelay));
                %Change rotation angles to PT prediction
                if parXT.PT.useForB0
                    fprintf('<strong>Pilot Tone</strong>: Using PT motion prediction for B0 modelling.')
                    TrOrig=E.Tr;
                    E.Tr = TrPTBack;[E] = transformationConversion(E);
                end    
                EnBefore = computeEnergy(yX,xRes,E,R,EH);
                E = GDsolverTaylor(yX,E,EH,P,CX,R,xRes,parXB.nGD_Taylor + 0*single(n==parXB.Optim.taylorDelay),tol(estT(l)+1), rec.Dyn.Debug);
                EnAfter = computeEnergy(yX,xRes,E,R,EH);
                if parXT.PT.useForB0
                    E.Tr = TrOrig;[E] = transformationConversion(E);
                    TrOrig=[];
                end
                if  rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end  
                %Unwrap linear maps
                if E.Dl.C.unwrapTaylor &&  mod(n,8)==1  ; for ii = 1:2; E.Dl.D = dynInd(E.Dl.D, ii,6, unwrapLinear(dynInd(E.Dl.D,ii,6), abs(xRes),0, E.Dl.TE)) ;end; end
                %Susceptibility based denoising
                if E.Dl.C.denoiseTaylor &&  mod(n,3)==1
                    wT = multDimMea(abs(E.Dl.Tr), 5); wT = dynInd(wT,5:6,6); wT = wT/normm(wT);
                    E.Dl.D = TaylorMapsDenoising(E.Dl.D, abs(xRes), wT, MTRes, 40,[],[],1);
                end
                %residualsPerState(yX,xRes,E,EH);
                if isfield(E,'Dl');plotND({abs(xRes),1,0},permute(E.Dl.D,[1:3 6 4 5])*pi/180, [-10 10 0 .5e6],[], 1,{1:3,2},MTRes,[],single(abs(xRes)>.5*multDimMea(abs(xRes))),{3},100 );end
            
            end
            %%% IMAGE ESTIMATION
            fprintf('** Image estimation - %d iterations \n', nX + 4*single(n==1) );
            EnBefore = computeEnergy(yX,xRes,E,R,EH);
            [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,xRes,nX + 8*single(n==1),tolType{estT(l)+1},tol(estT(l)+1));
            EnAfter = computeEnergy(yX,xRes,E,R,EH);
            if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end  
            plotND({angle(xRes),1,0}, abs(xRes),[0 .5e+6 -pi pi],[],0,{[],2},MTRes,[],[],[],123);
            %plotND({abs(xRes),1,0}, angle(exp(+1i*2*pi*E.Shim.TE*E.Shim.B0)),[-pi pi 0 1.4],[],0,{[],2},MTRes,[],[],[],124);

        end
        EnAfter = computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1)/numel(yX);ReSeg =  res_segments(EnAfter, E.NMs, E.nEc);
        %figure(4566); plot(1:E.NMs, ReSeg);
    
        if ~isempty(EnX)%Print and store final energy
            EnX=EnX/numel(yX);
            if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
            EnX=[];
        end

        %ENERGY COMPUTATION
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn; E.Sf=gpuArray(E.Sf);end
        yX=dynInd(ySt,ncx,4);
        [Re{l}(:,n),~,fullResid]=computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1);Re{l}(:,n)=Re{l}(:,n)/numel(yX);
        En{l}(:,n)=normm(yX,[],setdiff(1:numDims(yX),1))/numel(yX);
        if rec.Dyn.Debug>=2;fprintf('Residuals: %0.6g\n',sum(Re{l}(:,n)));end        

        %ENERGY PER STATES
        fullResid=gather(dynInd(fullResid,iSt{l},1,fullResid));
        Residuals=Re{l}(:,n);Residuals(iSt{l})=Residuals;
        Residuals=reshape(Residuals,NYres([1:2 5]));
        
        mStSort=mStSweeps;mStSort(iSt{l})=mStSort;%YB: now in mStSort the samples correspond to the right sweep
        Sweeps=reshape(mStSort,NYres([1:2 5]));%YB: So if we reshape we have the PE plane with the sweep number per TR         

        %OUTLIER DETECTION
        if leff>0;percUse=(0.15:0.05:0.85)*100;
        else percUse=(0.3:0.05:0.7)*100;
        end
        WeO=zeros([NSweeps length(percUse)],'like',Residuals);
        for s=1:NSweeps;WeO(s,:)=prctile(log(Residuals(Sweeps==s)),percUse);end
        WeO=permute(WeO,[1 3 2]);           
        Westd=diff(prctile(WeO,parXT.percRobustShot*100,1),1,1)./diff(norminv(parXT.percRobustShot)); 
        Wmea=prctile(WeO,mean(parXT.percRobustShot)*100,1);
        Wmea=Wmea+Westd*sqrt(2)*erfcinv(2*mean(parXT.percRobustShot));%sqrt(2)*erfcinv(2*mean(parXT.percRobustShot))=norminv(mean(parXT.percRobustShot))
        We=bsxfun(@minus, WeO,Wmea);%YB changed for version R2015
        We=bsxfun(@rdivide, We, Westd);%YB changed for version R2015
        We=mean(We,3);         
        We=min((1-normcdf(We))/((1-parXT.enerRobustShot)/NSweeps),1);     
        outlWe=We<1; %YB: outlWe flags which sweep is an outlier
        %visResiduals(We,outlWe,timeSweep,[],[],0);%YB see every iteration

        if rec.Alg.WriteSnapshots && estT(l)==0;visResiduals(We,outlWe,timeSweep,resampling(Residuals,NY(1:2),1),kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end
        %CS-LIKE/ROBUST RECONSTRUCTION
        if leff>0 && estT(l)==0
            subT(:)=1;
            if rec.Alg.AlignedRec==3 %YB Changed because AlignedRec==2 would also be subdivided
                outlWeV=find(outlWe);
            elseif rec.Alg.AlignedRec==4
                outlWeV=1:length(outlWe);
            else
                outlWeV=[];
            end%YB: AlignedRec==4 explore all, AlignedRec==3, explore only artifacted, AlignedRec==[ 1 2] no subdivision
            for o=1:length(outlWeV);subT(sweepToSample{outlWeV(o)})=2;end

            if all(subT(:)==1) || nowithin;L=l;end%If full resolution and final image recon +no outliers detected --> stop here
        end
        
        if (leff>0 && estT(l)==0) || (rec.Alg.AlignedRec==2 && leff>=0) %YB: WHY? Why only AlignedRec==2 and not if AlignedRec==3?
            EH.We=ones([size(yX,1) 1],'like',Residuals); %THIS COULD BE PROBLEMATIC IN CASE OF SEVERAL REPEATS?
            if ismember(rec.Alg.AlignedRec,[3 5])
                for s=1:NSweeps;EH.We(mStSweeps==s)=1-single(outlWe(s));end % YB: This makes EH.We have the same size(y) but remember that weights are only determined for segments, not for subdivisions!
            else%AlignedRec==4 or AlignedRec==1,AlignedRec ==2: use soft weights 
                for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end
            end
        end
        if leff>0 && (estT(l)==0 || n==rec.Alg.nExtern(l))
            if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xRes,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Di-l=%d',l)));end
            if parXT.computeCSRecon
                fprintf('** Regularised image estimation - %d iterations \n', nX );
                if nX~=0;xCS=CSsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1),rec.Dyn.Debug);else xCS=xRes;end
                if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xCS,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected regularized ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Re-l=%d',l)));end
            end
            EH.We=[];%YB: change so you can also have EH.We for CGsolver for robust (not CS) reconstruction
        end%TODO / WHY: ask Lucilio how we can streamline this

        if parXB.useSH==1 % Restore the binned coefficients
            line = 1:length(E.Db.c_idx);
            for iii = 1:max(E.Db.c_idx(:)); E.Db.c = dynInd(E.Db.c, line(E.Db.c_idx==iii), 5, repmat( dynInd(E.Db.cr,iii,5), [on{4} multDimSum(E.Db.c_idx==iii) ])) ;end; iii=[];
            E.Db.cr = E.Db.c;
            rec.Db.c = E.Db.c;
        elseif parXB.useSH==2
            line = 1:length(E.Dbs.c_idx);
            for iii = 1:max(E.Dbs.c_idx(:)); E.Dbs.c = dynInd(E.Dbs.c, line(E.Dbs.c_idx==iii), 5, repmat( dynInd(E.Dbs.cr,iii,5), [on{4} multDimSum(E.Dbs.c_idx==iii) ])) ;end; iii=[];
            E.Dbs.cr = E.Dbs.c;
            rec.Dbs.c = E.Dbs.c;
        end
        
        %%% PILOT TONE -- PTHandling
        %Set flag to subdivide segments
        if  parXT.PT.subDiv(l)>0 && (n==rec.Alg.nExtern(l) || estT(l)==0) %multDimSum(parXT.PT.usePT)>0 &&
            subT(:)=1;
            outlWeV=1:length(outlWe);%All of them are subdivided. TODO: Make compatible that only those with PT weights are subdivided (will only save time I think)
            for o=1:length(outlWeV);subT(sweepToSample{outlWeV(o)})=2^parXT.PT.subDiv(l);end
            estT(l)=0;
            fprintf('<strong>Pilot Tone</strong>: subdividing activated and estT(%d) set tot 0.\n',l)
        end
        
        %Calibrate PT data %PTHandling
        if multDimSum(parXT.PT.usePT)>0  &&  ( parXT.PT.usePT(l)>0 || ( l<L && multDimSum(parXT.PT.usePT(l+1:end)>0) || n==rec.Alg.nExtern(l) ) ) && n>2
            tStaPTCalib = tic;
            if parXT.PT.useRASMotion; [~,TrPTCalib] =  transformationConversion([],E.Tr,ssPTConversion,[],[], 1, size(x));else;TrPTCalib=E.Tr;end%Forward conversion:IJK->RAS
            %Weighted calibration
            if parXT.PT.weightedCalibrationFlag 
                if strcmp(parXT.PT.weightedCalibrationMethod,'WePT');WePTCalib=rec.PT.WePT;
                elseif strcmp(parXT.PT.weightedCalibrationMethod,'outlWePT');WePTCalib=1-single(rec.PT.outlWePT);
                end
                fprintf('<strong>Pilot Tone</strong>: Using %s for weighted calibration.\n',parXT.PT.weightedCalibrationMethod)
            else
                WePTCalib=[];
            end
            [rec.PT.Af,rec.PT.condAf]  = pilotToneCalibration ( rec.PT.pTimeResGr, TrPTCalib, 'forward', parXT.PT.calibrationOffset, WePTCalib);
            [rec.PT.Ab,rec.PT.condAb] = pilotToneCalibration ( rec.PT.pTimeResGr, TrPTCalib, 'backward', parXT.PT.calibrationOffset, WePTCalib);
            fprintf('<strong>Pilot Tone</strong>: Calibrating PT matrix both forward and backward way.\n')
            fprintf('<strong>Pilot Tone</strong>: Forward calibrating matrix conditioning number: %.1f.\n', rec.PT.condAf)
            fprintf('<strong>Pilot Tone</strong>: Backward calibrating matrix conditioning number: %.1f.\n', rec.PT.condAb)
            %pTimeResGrPredForw = pilotTonePrediction (rec.PT.Af , E.Tr, 'T2PT','forward',rec.PT.NChaPT,size(rec.T,6));
            %pTimeResGrPredBackw = pilotTonePrediction (rec.PT.Ab , E.Tr, 'T2PT','backward',rec.PT.NChaPT,size(rec.T,6));
            %plotPTSignalsToCompare (rec.PT.pTimeResGr, pTimeResGrPredForw, pTimeResGrPredBackw, 765);
            
            TGrPredForw = gather(pilotTonePrediction (rec.PT.Af , rec.PT.pTimeResGr, 'PT2T','forward',rec.PT.NChaPT,size(rec.T,6)));
                if parXT.PT.useRASMotion; [~,TGrPredForw] =  transformationConversion([],TGrPredForw,ssPTConversion,[],[], 0, size(x));end           
            TGrPredBackw = gather(pilotTonePrediction (rec.PT.Ab ,rec.PT.pTimeResGr, 'PT2T','backward',rec.PT.NChaPT,size(rec.T,6)));
                if parXT.PT.useRASMotion; [~,TGrPredBackw] =  transformationConversion([],TGrPredBackw,ssPTConversion,[],[], 0, size(x));end           
            plotPTSignalsToCompare (E.Tr, TGrPredForw, TGrPredBackw, 766);
            
            %plotPTSignals(rec.PT.pTimeResGr, 767, 'Grouped PT signal');plotPTSignals(pTimeResGrPredForw, 768, 'PredictedGrouped  PT signal');
            fprintf('PT calib took %.1f seconds.\n',toc(tStaPTCalib) );
        end
        
        %Implicit calibration
        if parXT.PT.implicitCalibration && leff==0
            if parXT.PT.useRASMotion>0; parXT.PT.useRASMotion=0; warning('Implicit calibration with RAS motion params not implemented. parXT.PT.useRASMotion set to 0.');end
            E.Tr=rec.T;E.Fs=FSt; % YB: Use the original (non-compressed motion) because the compression is for the Xsolver, but not the Tsolver.
            if isfield(E,'Dl') || isfield(E,'Ds');E = transformationConversion(E);end%Make motion parameters for linear model 
            E.NMs=size(E.Tr,5);E.NSe=length(E.Fs{1});
            if parXT.exploreMemory;E.cT(:)=1;end%To test flow without running the main methods       

            %EXTRACT RELEVANT ROI AND COILS AND CALL THE SOLVER                       
            nct=1:eivaS(2); %YB: reminder, the second element gives the #coils to use for motion estimation
            E.Sf=dynInd(SRes,nct,4);
            ySf=dynInd(ySt,nct,4);
            E.dS=[E.NMs eivaS(2)];
            if gpuIn;E.Sf=gpuArray(E.Sf);end
            if parXT.fractionOrder~=0;E.Fd=GRes;end
            [E, rec.PT] = GDsolverPTCalibration(ySf,E,EH,[],CT,[],xRes,5,[],2, rec.PT);
            rec.PT.condAb=cond(rec.PT.Ab);
            fprintf('<strong>Pilot Tone</strong>: Implicit backward calibration matrix conditioning number: %.1f.\n', rec.PT.condAb)
            rec.T=gather(E.Tr);
            if isfield(E,'Dl') || isfield(E,'Ds');[E] = transformationConversion(E);end
            THist=cat(1,THist,gather(rec.T));
            visMotion(rec,[],[],0,[],[],[],[],[],456);
            %if isfield(E,'Dl');visMotion(E.Dl.Tr,[],[],0,[],[],[],[],[],457);end
        end
        
        % MOTION UPDDATE 
        if estT(l) && ~(n==rec.Alg.nExtern(l))  
            %LM SOLVER - note EH.We never used in LMsovler.m since it is a separate Least Squares
            fprintf('** Motion estimation\n');
            %TRANSFORM PARAMETERS PER SEGMENT/HARMONICS/BLOCK SIZES
            E.Tr=rec.T;E.Fs=FSt; % YB: Use the original (non-compressed motion) because the compression is for the Xsolver, but not the Tsolver.
            if isfield(E,'Dl') || isfield(E,'Ds');E = transformationConversion(E);end%Make motion parameters for linear model 
            E.NMs=size(E.Tr,5);E.NSe=length(E.Fs{1});
            if parXT.exploreMemory;E.cT(:)=1;end%To test flow without running the main methods       

            %EXTRACT RELEVANT ROI AND COILS AND CALL THE MOTION SOLVER                       
            nct=1:eivaS(2); %YB: reminder, the second element gives the #coils to use for motion estimation
            E.Sf=dynInd(SRes,nct,4);
            ySf=dynInd(ySt,nct,4);
            E.dS=[E.NMs eivaS(2)];
            if gpuIn;E.Sf=gpuArray(E.Sf);end
            if parXT.fractionOrder~=0;E.Fd=GRes;end
            if ~(parXT.PT.implicitCalibration && leff==0); [E,xRes]=LMsolver(ySf,E,xRes,CT);end%YB:Only when not estimating motion implicitely
            if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end  
            rec.T=E.Tr;
            if isfield(E,'Dl') || isfield(E,'Ds');[E] = transformationConversion(E);end
            THist=cat(1,THist,gather(rec.T));
            visMotion(rec,[],[],0,[],[],[],[],[],456);
            if isfield(E,'Dl');visMotion(E.Dl.Tr,[],[],0,[],[],[],[],[],457);end
        else
            if l==L && parXT.saveFinal
                recW=rec;recW.Alg.SaveRaw=3;recW.Dyn.Batch=1;
                recW.E=E;recW.EH=EH;recW.P=P;recW.CX=CX;recW.R=R;recW.xRes=xRes;recW.yX=yX;recW.fullResid=fullResid;
                recW.timeState=timeState;recW.We=We;recW.Residuals=Residuals;recW.Sweeps=Sweeps;recW.iSt=iSt{l};
                recW.S=[];recW.y=[];
                writeRaw(recW);recW=[];
            end
            break%YB: If est(l)==0 then break;
        end
        E.Sf=[];

        %CHECK CONVERGENCE
        if all(E.cT(:))%Weak convergence
            %if l==L;estT(l)=0;else break;end%If last level we perform another estimation for x           
            estT(l)=0;            
            if rec.Alg.WriteSnapshots;visMotion(rec,voxSiz,timeState,2,strcat(folderSnapshots,filesep,'Motion'),strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end                
        end
    end

    %RECONSTRUCTION TO BE USED TO INITIALIZE THE NEXT RESOLUTION LEVEL
    x=resampling(xRes,NX);
    if ~parXT.exploreMemory && isfield(E,'Ds') && isfield(E.Ds,'chi') && l~=L ;DC=E.Ds.chi;end %YB: save for next level  
    
    if gpuIn;rec.M=gpuArray(rec.M);end
    if discMa;x=rec.M.*x;end
    
    if leff>0
        rec.d=dynInd(rec.d,leff,4,gather(x)); %YB: solution for different levels are stored in 4th dimension
        if parXT.computeCSRecon
            if discMa;xCS=rec.M.*xCS;end
            rec.r=dynInd(rec.r,leff,4,gather(xCS));
        end
    end
    if ~isempty(En{l})
        EnSort=En{l};EnSort(iSt{l},:)=EnSort;
        ReSort=Re{l};ReSort(iSt{l},:)=ReSort;  
        mStSort=E.mSt;mStSort(iSt{l})=mStSort;
        rec.Par.Mine.DISORDER.Residuals{l}=gather(reshape(ReSort,[NYres([1:2 5]) size(ReSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Energy{l}=gather(reshape(EnSort,[NYres([1:2 5]) size(EnSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Motions{l}=THist;
        rec.Par.Mine.DISORDER.States{l}=gather(reshape(mStSort,NYres([1:2 5])));
    end

    if parXT.writeInter && leff>0 || l==L 
        % YB added here before E = [];
        rec.timeState=timeState;
        if (parXB.dephaseCorrection>0 || parXB.useTaylor) && ~parXT.exploreMemory %YB: added  && ~parXT.exploreMemory
                if parXB.dephaseCorrection>0; rec.D=gather(E.Dc.D); else rec.D=gather(E.Dl.D);  end
                rec.D = permute(rec.D, perm);
                %rec.D =extractROI(rec.D,rec.Enc.ROI,0,1); %YB oes not work
                for t=typ2RecF'
                    if ~any(rec.Dyn.Typ2Rec==t);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,t);rec.Dyn.Typ2Wri(t)=1;end
                end
                typ2Rec=rec.Dyn.Typ2Rec;
        end
        rec.E = E;%YB: also export E (delete later as this might require a lot of memory)
        E=[];EH=[];

        drec=rec.d;
        if parXT.computeCSRecon;rrec=rec.r;end
        ROIrec=rec.Enc.ROI;          
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if gpu;rec.(datTyp)=gpuArray(rec.(datTyp));end
            if n~= 12; rr = 1:leff;  else rr = 1;end %YB added because rec.x only 3D
            rec.(datTyp)=dynInd(rec.(datTyp),rr,4);
            rec.(datTyp)=flip(rec.(datTyp),4);
        end
        %if l==L && leff~=1;typ2RecL=setdiff(typ2Rec,[5;12]);elseif l==L && leff==1;typ2RecL=setdiff(typ2Rec,5);else typ2RecL=typ2RecI;end
        if l==L && leff~=1;typ2RecL=setdiff(typ2Rec,[5, 29]);elseif l==L && leff==1;typ2RecL=setdiff(typ2Rec,[5,29]);else typ2RecL=typ2RecI;end%YB change
        %FINAL ADJUSTMENTS
        for n=typ2RecL';datTyp=rec.Plan.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),perm);
            rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0,1);
        end
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if rec.Alg.MargosianFilter;rec.(datTyp)=margosianFilter(rec.(datTyp),rec.Enc);end
            rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);%REMOVE OVERDECODING  
        end
        rec.Enc=rmfield(rec.Enc,'ROI');
        rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file
        typ2RecI=setdiff(typ2RecI,12);
        writeData(rec);
        if l~=L
            rec.Enc.ROI=ROIrec;            
            rec.d=drec;
            if parXT.computeCSRecon;rec.r=rrec;end
        end
    end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time till level %d: %.3f seconds = %.3fminutes\n',l,tend,tend/60);end
    if l==L
        if rec.Dyn.Debug>=2; diary off;end %YB Stop logging
        break;
    end%Early termination if no outliered states detected
    
    
end
