
function rec=solveSensit7T(rec, recB)

%SOLVESENSIT7T Solves for the sensitivities in 7T data by building a virtual body coil.
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and eigenmaps rec.W
%

if nargin < 2; recB=[];end
gpu=gpuDeviceCount;

%HARDCODED GENERAL PARAMETERS
% resolESPIRiT=6;
% resol=max(resolESPIRiT,max(rec.Enc.AcqVoxelSize));%Resolution of sensitivities
resol=max(rec.Enc.AcqVoxelSize);%Resolution of sensitivities % ZN: use the max of reference line resol for espirit
rec.Alg.parE.nItESPIRIT=1;%Number of iterations of ESPIRIT refinements
rec.Alg.parS.nItVARIATI=2;%Number of iterations of variational refinements
ND=3;%Number of dimensions for estimation
rec.Alg.parS.nItAniDiff=0;%Number of iterations of anisotropic diffusion
rec.Dyn.Debug=1;
rec.Dyn.Typ2Rec=[];
rec.Alg.writeSnapshots =1; %Use for reporting ESPIRiT results (see Miki Lustig's code)

NY=size(rec.y);NY(end+1:ND)=1;
rec.Par.Mine.pedsUn=1;%YB:??
if ~isfield(rec,'Enc') || ~isfield(rec.Enc,'AcqVoxelSize');rec.Enc.AcqVoxelSize=ones(1,ND);end
rec.Par.Labels.SliceGaps=0;

%%% RESAMPLE TO REQUIRED RESOLUTION
NResol=round(NY(1:ND).*rec.Enc.AcqVoxelSize(1:ND)/resol);
%rec.y=filtering(rec.y,buildFilter(NY(1:ND),'tukeyIso',1,gpu,0.05));
rec.y=resampling(rec.y,NResol);
[rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec] = mapNIIGeom( rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec, 'resampling',[], NY(1:3), NResol(1:3));
rec.Enc.FOVSize = NResol(1:3);

if gpu;rec.y=gpuArray(rec.y);end
%rec.y=filtering(rec.y,buildFilter(NResol,'tukeyIso',1,gpu,0.05),1);

%%% CONSTRUCT BODY COIL
blockBodyCoilUsage=0;
if ~blockBodyCoilUsage && ~isempty(recB)%Over-write the body coil data from the array receivers
    fprintf('Using body coil from separate body-coil acquisition.\n') 
    assert(size(recB.y,4)==2,'Body coils are not the 2 quadrature modes');
%recB.y = mapVolume(recB.y, rec.y, recB.Par.Mine.APhiRec,rec.Par.Mine.APhiRec,[],[],'spline');    
    B = dynInd(recB.y,1,4) - 1i*dynInd(recB.y,2,4);%Quadrature combination
    phU = CNCGUnwrapping(B,recB.Enc.AcqVoxelSize, 'Magnitude','LS');%Use voxelsize of the original array since mapped to it
    M = abs(B)>.2*multDimMax(abs(B));
    ph = phU;% - multDimMea(phU(M));
    x = abs(B).*exp(1i*ph);
x = mapVolume(x, rec.y, recB.Par.Mine.APhiRec,rec.Par.Mine.APhiRec,[],[],'spline');
B = mapVolume(B, rec.y, recB.Par.Mine.APhiRec,rec.Par.Mine.APhiRec,[],[],'spline');
    rec.x = x;
    rec.B = abs(x);
    rec.Alg.parE.nItESPIRIT=1;%Don't use the iterative ESPIRiT
    rec.Alg.parS.nItVARIATI=1;    
    plotND([], cat(4,RSOS(rec.y),matchHist(RSOS(B),RSOS(rec.y))),[],[],0,[],rec.Par.Mine.APhiRec,{'RSOS array receiver';'RSOS body receiver'},M,{2,[],[],inf},90,'Array vs reveiver comparison');
else
   fprintf('Using virtual body coil from array compression.\n') 
    rec.x=compressCoils(rec.y,1);
    rec.B=RSOS(rec.y);%YB: This is the RootSumOfSquares(RSOS) - does not have phase!
    rec.B=anidiffND(rec.B,rec.Alg.parS.nItAniDiff,[]);
    %rec.x=rec.B.*exp(1i*CNCGUnwrapping(rec.x,rec.Enc.AcqVoxelSize, 'Magnitude','LSIt'));%YB: modified with new function
    %rec.x=rec.B.*exp(1i*CNCGUnwrapping(rec.x,rec.Enc.AcqVoxelSize, 'Magnitude','LS'));%YB: modified with new function
    rec.x=rec.B.*exp(1i*CNCGUnwrappingOld(rec.x,[],rec.Enc.AcqVoxelSize));%YB: assign phase of the compressed coil to the RSOS image

    %Now x is has the phase of the compressed coil and magnitude of the RSOS image. Note the phase can contain singularities and the phase unwrapper can be used to filter it! ('LS' instead of 'LSIt')
    %This is a first guess since phase will be update in the iterative ESRPIRiT
end

%FIRST GUESS ON SENSITIVITIES
if rec.Dyn.Debug && rec.Alg.parS.nItVARIATI>0; fprintf('\n--- Estimating coil sensitivities: Variational method ---\n');end
for l=1:rec.Alg.parS.nItVARIATI
    rec.Alg.parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
    rec.Alg.parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
    rec.Alg.parS.Otsu=[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)    
    rec.Alg.parS.lambda=1e6;%00000;%Last parameter has been recently changed, toddlers have been processed with it set to 1, but it was introducing noise for neonates, 2 may suffice, but 5 gives a bit of margin
    rec.Alg.parS.order=2;%Regularization order for coil estimation-Something around 0.2 is optimal for order 2
    rec.Alg.parS.nErode=0;%2;%Erosion for masking (in mm)
    rec.Alg.parS.nDilate=3;%Dilation for masking (in voxels then in mm)
    rec.Alg.parS.conComp=2;%Whether to get the largest connected component after erosion (1) and to fill holes (2)
    rec.Alg.parS.GibbsRingi=[0 0];%Gibbs ringing for coil profiles
    rec.Alg.parS.ResolRatio=[1 1];%Resolution ratio for coil profiles
    rec.Alg.parS.TolerSolve=1e-5;%S-solver tolerance
    rec.Alg.parS.nIt=300;%S-solver maximum number of iterations
    rec.Alg.parS.softMaskFactor=10;%Factor to normalize the body coil data for soft masking 

    if isfield(rec,'S');rec=rmfield(rec,'S');rec=rmfield(rec,'M');end
    rec.x=rec.B.*sign(rec.x);%YB: .*sign() is adding the phase of the argument
    rec=solveS(rec);
    if isempty(recB);rec.x=SWCC(rec.y,rec.S);end %bsxfun(@times,sum(bsxfun(@times,conj(rec.S),rec.y),4),1./(normm(rec.S,[],4)+1e-6));%Image estimate based on sens estimate    
end

%RECURSIVE ESPIRIT
if rec.Dyn.Debug && rec.Alg.parE.nItESPIRIT>0; fprintf('\n--- Estimating coil sensitivities: ESPIRiT ---\n');end
%rec.B does not change and stays the RSOS: over iterations, the phase is changed!
rec.Alg.parE.NCV=1;%Maximum number of eigenmaps
rec.Alg.parE.NC=resol*ones(1,ND);%[4 4 4];%Resolution (mm) of calibration area to compute compression
rec.Alg.parE.K=100*ones(1,ND);%Resolution (mm) of target coil profiles - kernel size =(1./parE.K)./DeltaK  with  DeltaK=(1./(N(1:ND).*voxsiz(1:ND)));%1/FOV
rec.Alg.parE.subSp=resol*ones(1,ND);%Subsampling in image space to accelerate
rec.Alg.parE.mirr=0*ones(1,ND);%Whether to mirror along a given dimension
rec.Alg.parE.Ksph=round(150^(ND/3));%Number of points for spherical calibration area, unless 0, it overrides other K's
rec.Alg.parE.eigTh=0.03;%0.01;%Threshold for picking singular vectors of the calibration matrix (relative to the largest singular value). YB:If <=0 will be automatically detected
% Zihan % was 0.03

rec.Alg.parE.absPh=0;%Flag to compute the absolute phase using virtual conjugate coils - see paper
rec.Alg.parE.virCo=1;%Flag to use the virtual coil to normalize the maps%YB:What do other values mean?
rec.Alg.parE.eigSc=[0.85 0.3];%0.25;%YB:Not used in solveESPIRIT.m
rec.Alg.parE.dimLoc=[];%Dimensions along which to localize to compute virtual coils YB: The method BerkinBiglic proposed?
rec.Alg.parE.Kmin=6;%Minimum K-value before mirroring%YB:??
rec.Alg.parE.saveGPUMemory=1;
rec.Alg.parE.gibbsRinging=0.5; % Zihan % was 0.5
rec.Alg.parE.showresult = 0; % ZN: 0 - do not show the mid-result or writesnapshots
for s=1:rec.Alg.parE.nItESPIRIT
    rec.x=rec.B.*sign(rec.x);%YB: here you add the phase of the iamge estimate to the RSOS to create body coil
    rec.Alg.parE.nIt = s;
    rec=solveESPIRIT(rec);
    %rec=solveESPIRIT_BART(rec);
    if gpu;rec.S=gpuArray(rec.S);rec.y=gpuArray(rec.y);rec.W=gpuArray(rec.W);rec.x=gpuArray(rec.x);end
    rec.S = dynInd(rec.S,1,6);%Only take first sensitivities - multiple might be computed to validate the choice of one eigenmap
    if isempty(recB);rec.x= SWCC(rec.y,rec.S);end%YB: create new estimate of image with new senisitivities estimate
end

rec = gatherStruct(rec);
% rec.S=gather(rec.S);%Coil sensitivities
% rec.y=gather(rec.y);%Raw coil images
% rec.W=gather(rec.W);%Eigenmaps
% rec.x=gather(rec.x);%Reconstructed image
% rec.B=gather(rec.B);%RSOS that served as virtual body coil

%%% PLOT PARAMETERS FOR LOGGING
rec.Alg.parE
rec.Alg.parS

end