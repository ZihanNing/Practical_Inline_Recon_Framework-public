function rec=reconAlgorithm(rec)

%RECONALGORITHM   States the algorithmic parameters for reconstruction
%   REC=RECONALGORITHM(REC)
%   * REC is a reconstruction structure without the algorithmic parameters. 
%   At this stage it may contain the naming information (.Names), the 
%   status of the reconstruction (.Fail), the .lab information (.Par), the 
%   fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   ** REC is a reconstruction structure with filled algorithmic 
%   information (.Alg)
%

Alg=[];

Alg.JointSplitScans=1;%To estimate split scans together, problem is that the bed may have moved, so the isocenter may not coincide... 
if nargin<1;rec=Alg;return;end
if rec.Dyn.Debug>=1;fprintf('Study detected as: %s\n',rec.Par.Mine.curStud.Stud);end

%GENERIC PARAMETERS
Alg.UseSBSensi=1;%To estimate and use the SB sensitivities, 1 to use only the SB to mask the MB, 2 to re-estimate the sensitivities, 3 to re-estimate using ESPIRIT and outlier robustness, 4 for ESPIRIT only
if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiPha');Alg.UseSBSensi=2;end
%ESPIRIT AND ROBUST (FETAL CASE AND FMRI ONLY)
if ( ~isempty(regexp(rec.Names.Name,'sf','ONCE')) || ~isempty(regexp(rec.Names.Name,'fd','ONCE')) || ~isempty(regexp(rec.Names.Name,'f[5-9]','ONCE')) ) && rec.Par.Mine.Modal==9;Alg.UseSBSensi=3;end%WE ASSUME THIS IS A FETAL CASE, IT MAY NOT
Alg.solverType='CG';%'CG';%'IRWLS';%Solver type. One of the following: 'CG','IRWLS'
Alg.UnfoldFat=0;%To unfold the fat (2) or activate lower eigenmaps (1)
%if rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1 && Alg.UseSBSensi==3 && rec.Plan.Quick<3
if Alg.UseSBSensi==3 && rec.Plan.Quick<3
    Alg.solverType='IRWLS';%'IRWLS'
    if Alg.UseSBSensi==3;Alg.UnfoldFat=1;end
end
if Alg.UseSBSensi==2 && rec.Par.Mine.Modal==9;Alg.UseSBSensi=4;end%ESPIRIT FOR FMRI

Alg.artifRobust=0;%To activate artifact-robust reconstructions
Alg.artifSigma=0.5;%Level of artifact-robustness
Alg.IRWLSPower=-1.5;%-1 for standard Huber, -2 for Huber with full saturation, positive for lp.

Alg.ReestSensi=0;%To reestimate the sensitivities
Alg.NoiseStand=1;%To standardize noise
Alg.DistoSensi=0;%To use distorted sensitivities for EPI

Alg.EstimGFact=[0 1];%To estimate the g-factor (first element, 1 for quick analytic, 0.5 for analytic,  >1 for Monte-Carlo NR-1 repeats)
if ~rec.Plan.Quick
    if rec.Par.Mine.Modal==10;Alg.EstimGFact=[0.5 1];%DWI denoising basic
    elseif rec.Par.Mine.Modal==5;Alg.EstimGFact=[2 1];%T2 reconstruction basic, second element spatial regularization kernel
    end
end
%if rec.Par.Mine.Modal==6;Alg.EstimGFact=1;end
%if rec.Par.Mine.Modal==7;Alg.EstimGFact=1;end
Alg.EstimCFact=1;%To estimate the chi2-factor
if rec.Par.Mine.Modal==7;Alg.EstimCFact=0;end%MEMORY PROBLEMS IN THIS CASE

Alg.WriteSnapshots=1;%Write snapshots
Alg.WriteSensi=0;%To write sensitivity information
Alg.StabilFreq=1;%1 for stabilizing frequency as pre-processing

Alg.RemoveChannels=[];%To remove a series of channels from the reconstruction
%Alg.RemoveChannels=18;%To remove a series of channels from the reconstruction
%Alg.RemoveChannels=2;%[];%To remove a series of channels from the reconstruction

Alg.PhasCorRef=0;%Indicates whether to correct for phase mismatches in multishot DWI using navigator data
Alg.GhosCorRef=1;%Indicates whether to correct for ghosting using calibration data
%ALG.GHOSCORREF TO 0 TO RECONSTRUCT JANA'S DATASETS
if rec.Par.Mine.Modal==9;Alg.GhosCorRef=2;end%For fMRI has shown certain advantages in Giulio's fetal datasets, not much difference for neonates
Alg.GhosCorFilter=0;
if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra');Alg.GhosCorFilter=1;end%Indicates whether to filter the ghosting information in the slice direction, necessary for split scans where the isocenter may change and no reference is available...
Alg.DetectAnomalies=0;%0;%Activate to run the code to detect anomalies in the data.
                      %1-> spikes, sort data-NOT TESTED FOR A WHILE
                      %2-> spikes, invert data, explore first batch detect single dynamic
                      %3-> spikes, invert data, explore all batchs, detect single dynamic
                      %4-> spikes, invert data, explore all batchs, detect all dynamics
                      %5-> detuning, invert data, explore first batch detect single dynamic
                      %6-> detuning, invert data, explore all batchs, detect single dynamic
                      %7-> detuning, invert data, explore all batchs, detect all dynamics
Alg.CheckGain=0;%Activate to run the code to check the gains
Alg.SaveRaw=0;%Activate to save the raw data before reconstruction either non-anonymized (1) or anonymized (2), after reconstruction either non-anonymized (3) or anonymized (4) or after reading the profiles (5)

Alg.UsePrevGhos=1;%To connect ghosting information among sequences
if Alg.GhosCorRef==0;Alg.UsePrevGhos=0;end

if ismember(rec.Par.Mine.Modal,9:10);Alg.PE=1;else Alg.PE=2;end%Phase encode dimension of the reconstruction
Alg.GhosCorRec=0;%To call the Nyquist ghost corrected reconstruction

%if rec.Par.Mine.Modal==10 && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1 && Alg.UseSBSensi~=3;Alg.GhosCorRec=1;end%To call the Nyquist ghost corrected reconstruction
%if rec.Par.Mine.Modal==10 && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1;Alg.GhosCorRec=1;end%NOT A STRONG DIFFERENCE

%HOW DO WE SET THIS??? CHECK WITH MB FMRI
Alg.TolerSolve=1e-3;%Tolerance for solving the reconstruction problem---note it was 1e-5 for maximum performance (in noise estimation)
if ~isempty(rec.Par.Mine.AdHocArray) && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1;Alg.TolerSolve=1e-4;end

Alg.PlugNoise=0;%To plug noise for g-factor tests. 0-> data reconstructions. 1/2-> noise before/after readout transform. 3/4-> noise before/after phase encode transform. 5-> noise in the solveX function. 6-> noise in the processing methods
if rec.Par.Mine.Modal==2;Alg.PlugNoise=0;end%For reference estimation this does not make sense
Alg.SVDRecover=1;%To recover the data using patch-based SVD shrinkage
%Alg.SVDRecover=0;%To recover the data using patch-based SVD shrinkage
Alg.MargosianFilter=1;%To perform Margosian correction
if rec.Par.Mine.Modal~=10;Alg.SVDRecover=0;end%Force not to filter if not DWI

%if ismember(rec.Par.Mine.Modal,9:10);rec.Par.Mine.OverP=2;end
Alg.OverDec=[1 1 1];%To overdecode the reconstructions
if isfield(rec.Par.Mine,'OverP');Alg.OverDec(2)=rec.Par.Mine.OverP;end
if ~isempty(rec.Par.Encoding.KzOversampling);Alg.OverDec(3)=rec.Par.Encoding.KzOversampling(1);end
Alg.OnlyJSON=0;%To only extract and store JSON information
Alg.UseSoftMasking=1;%To use soft instead of hard masking. If 2 we don't mask
%Alg.OverDec(2)=-2;%NOTE THAT WITH DIFFERENT PE's THIS IS ONLY APPROXIMATED

%PARAMETERS FOR VISUALIZATION
Alg.parV.visGhost=0;%To visualize ghost correction term
Alg.parV.visMotion=0;%To plot motion estimation results
Alg.parV.visResiduals=0;%To plot information on the residuals
Alg.parV.visReconstruction=0;%To show an image reconstruction example
Alg.parV.visEmpiricalSamples=0;%To plot the empirical samples versus simulated histograms
Alg.parV.visSensitivities=0;%To plot the sensitivitie
if ismember(rec.Par.Mine.Modal,2);Alg.ParV.visSensitivities=0;end%To only observe if it is data

%PARAMETERS FOR COIL AND MASK ESTIMATION
Alg.parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
Alg.parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
Alg.parS.Otsu=[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)
if ~ismember(rec.Par.Mine.Modal,9:10);Alg.parS.lambda=20;else Alg.parS.lambda=5;end %Last parameter has been recently changed, toddlers have been processed with it set to 1, but it was introducing noise for neonates, 2 may suffice, but 5 gives a bit of margin
Alg.parS.order=2;%Regularization order for coil estimation-Something around 0.2 is optimal for order 2
Alg.parS.nErode=0;%2;%Erosion for masking (in mm)
Alg.parS.nDilate=3;%Dilation for masking (in voxels then in mm)
Alg.parS.conComp=2;%Whether to get the largest connected component after erosion (1) and to fill holes (2)
Alg.parS.GibbsRingi=[0 0];%Gibbs ringing for coil profiles
Alg.parS.ResolRatio=[1 1];%Resolution ratio for coil profiles
Alg.parS.TolerSolve=1e-5;%S-solver tolerance
Alg.parS.nIt=300;%S-solver maximum number of iterations
Alg.parS.softMaskFactor=10;%Factor to normalize the body coil data for soft masking%It was 20, have changed it to 10, not much difference

%PARAMETERS FOR ESPIRIT
if strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiPha');Alg.parE.NCV=2;else Alg.parE.NCV=3;end%Maximum number of eigenmaps, not enough memory on toddlers
if rec.Par.Mine.Modal==10;Alg.parE.NC=[4 4 4];else Alg.parE.NC=[2 2 2];end%Resolution (mm) of calibration area to compute compression
Alg.parE.K=[100 100 100];
Alg.parE.eigTh=0;%0.02%Threshold for picking singular vectors of the calibration matrix (relative to the largest singular value)
Alg.parE.absPh=0;%1%Flag to compute the absolute phase
Alg.parE.virCo=1;%Flag to use the virtual coil to normalize the maps
Alg.parE.eigSc=[0.85 0.3];%0.25;
Alg.parE.subSp=[1 1 1];%[1 1 1];%Subsampling in image space to accelerate
Alg.parE.dimLoc=[];%Dimensions along which to localize to compute virtual coils
Alg.parE.mirr=[8 8 8];%Whether to mirror along a given dimension
Alg.parE.lambda=50;%Regularization factor
Alg.parE.Kmin=6;%Minimum K-value before mirroring
Alg.parE.Ksph=200;%200 before%Number of points for spherical calibration area, unless 0, it overrides other K's

%PARAMETERS FOR STANDARD RECONSTRUCTION
Alg.parX.UseGiRingi=0;%To use Gibbs ringing before (1) or after (2) SENSE unfolding
Alg.parX.refineMask=0;%To refine masking (currently for dynamic scans)
Alg.parX.GibbsRingi=0.2;%Gibbs ringing before/after SENSE unfolding
if rec.Par.Mine.Modal==7
    Alg.parX.perc=0.9;%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved)
    Alg.parX.UseGiRingi=1;
    Alg.parX.GibbsRingi=0.05;
elseif ismember(rec.Par.Mine.Modal,5:6)
    Alg.parX.UseGiRingi=1;
    Alg.parX.GibbsRingi=0.4;%[0.2 0.2];
    Alg.parX.perc=0.99;
else
    Alg.parX.perc=[];
end

%PARAMETERS FOR ALIGNED RECONSTRUCTION (VOLUMETRIC)
if ~rec.Plan.Quick;Alg.AlignedRec=2;else Alg.AlignedRec=0;end%To call the aligned reconstruction (DISORDER) method
if (rec.Par.Mine.Modal~=7 && ~strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra')) || ~ismember(rec.Par.Mine.Modal,5:7) || (rec.Par.Mine.Modal==7 && strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra'));Alg.AlignedRec=0;end%Force not to align if not volumetric and not neonatal or not structural

Alg.parXT.NWend=1;%Number of within shot subdivisions at which to estimate motion, if multiple of 2 we assume tiling is used, otherwise rounded for temporal subdivisions
Alg.parXT.accel=[1 0];%Avoids motion estimation in the given number of levels of the spatio-temporal pyramid (to accelerate) - reconstruction in the given number of levels
%if Alg.AlignedRec==2;Alg.parXT.accel=[0 0];end
if Alg.AlignedRec>2;Alg.parXT.NWend=16;Alg.parXT.accel=[0 0];end%Full within-shot corrections
Alg.parXT.resolMax=4;%2;%4;%Coarser resolution on which to compute motion (in mm)
Alg.parXT.apod=[0.1 0 0];%Apodization factors
Alg.parXT.redFOV=1/3;%Factor to reduce the FOV in the inferior direction for motion estimation
Alg.parXT.perc=[0.9 0.9 0.9];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
Alg.parXT.traLimX=[0.05 0.02];%Dispersion limit of translations/rotations for binning motion states
Alg.parXT.traLimXT=[0.05 0.02];%Respectively translation and rotation limits for XT estimation
Alg.parXT.echoesMin=4;%Minimum required number of samples to compute motion
Alg.parXT.winit=1e-3;%Initial weight of the quasi-Newton update
Alg.parXT.meanT=0;%To constrain the parameters of the transforms solution to sum up to 0
Alg.parXT.UseGiRingi=1;%To use Gibbs ringing filter
Alg.parXT.GibbsRingi=0.05;%Gibbs ringing before SENSE unfolding
Alg.parXT.discardHighRes=0.75;%Those volumetric scans whose resolution is below this value (in mm) are not reconstructed
if Alg.AlignedRec==0;Alg.parXT.discardHighRes=0;end
Alg.parXT.convTransformJoint=0;%If 1 all transforms are estimated at each step, otherwise only those that have not converged
Alg.parXT.tolerSolve=1e-3;%Energy decrease in last iteration
Alg.parXT.percRobustShot=[0.125 0.25];%Percentiles for robust computation of expected inter-shot dispersion
Alg.parXT.enerRobustShot=0.95;%Ratio of the error for acceptance of shots
Alg.parXT.groupSweeps=1;%Factor to group sweeps
Alg.parXT.exploreMemory=0;%To explore convergence without running the main methods
Alg.parXT.writeInter=1;%To write intermediate data
Alg.parXT.computeCSRecon=1;%To compute final CS-like reconstruction if different from zero. It gives the resolution ratio of reconstructions with respect to baseline acquisition
Alg.parXT.fractionOrder=0;%Order for fractional derivative-based motion estimation
Alg.parXT.saveFinal=0;%Save final results for further inspection
Alg.parXT.maximumDynamics=4;%Maximum number of dynamics allowed

%PARAMETERS FOR ALIGNED RECONSTRUCTION (MULTISLICE)
Alg.parXT.toler=[10 10];%5e-7%Tolerance for convergence
if rec.Par.Mine.Modal==5;Alg.parXT.alpha=[0.25 1];else Alg.parXT.alpha=[0.25 1];end%For T2MS/T1MS
%if rec.Par.Mine.Modal==5;Alg.parXT.alpha=[0.25 10];else Alg.parXT.alpha=[0.25 10];end%For T2MS/T1MS
Alg.parXT.correct=1;%1;%1;%1;%0 does not correct for motion / 1 corrects for motion
Alg.parXT.threeD=1;%0 does not use the slice profile / 1 uses the slice profile / 2 uses slice links but not the slice profile
Alg.parXT.outlP=1.1;%1.2;%inf;%Oultier rejection threshold, inf for no outlier rejection
Alg.parXT.thplc=2;%0 does not correct for within plane motion / 1 does not correct for through-plane motion / 2 performs full corrections
Alg.parXT.NzSlab=7;%Size of slab to accelerate
Alg.parXT.refineSoftMask=1;%To refine soft masking
Alg.parXT.fractionReg=2;%2%Fraction order regularizer
Alg.parXT.outputResol=[];%Output resolution in plane, if empty native resolution
Alg.parXT.SVRoutputResol=0.5;%[];%Output resolution of SVR
%Alg.parXT.threeD=0;Alg.parXT.outlP=inf;Alg.parXT.thplc=1;Alg.parXT.alpha=[0 0];Alg.parXT.toler=1e-4;

%PARAMETERS FOR PREPROCESSING (DISTORTION ESTIMATION / MOTION CORRECTION)
%DISTORTION
Alg.parU.writeVideo=0;%To write a video with fMRI information
Alg.parU.useShimBox=1;%Whether to use shim box (1) or brain mask (0) for distortion correction
Alg.parU.Upsampling=1;%Upsampling for unwrapping
Alg.parU.GibbsRingi=1;%Gibbs ringing filter for smoothing the field 
Alg.parU.weightZ=1/10;%Weight of the slice information for unwrapping
Alg.parU.weightT=1/100;%Weight of the temporal information for unwrapping
Alg.parU.LUnwrap=1;%Norm for unwrapping-less than 1 unstable
Alg.parU.UnwrapMeth='CNCG';%'CNCG';%Phase unwrapping method, PUMA is based on graph cuts / LUCK is based on Laplacian equation solver (disabled) / CNCG is based in lp-norms
if rec.Par.Mine.Modal==10;Alg.parU.weightFunc='Magnitude';
else Alg.parU.weightFunc='Magnitude';
end%A combination of the following: 'Magnitude', 'Gradient'

%MOTION
Alg.parU.corrMotion=2;%0;%2;%2;%2;%To correct for motion. 0: no correction. 1: volumetric. 2: per-excitation
Alg.parU.useUndist=1;%0;%1;%1;%If 1->dynamic correction, if 2-> static correction given by the first volume
if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') || strcmp(rec.Par.Mine.curStud.IssuId,'DiscardDWI');Alg.parU.useShimBox=0;end
if rec.Par.Mine.Modal==10 && Alg.parU.corrMotion==2;Alg.parU.corrMotion=1;end%Force only volumetric for DWI
Alg.parU.fractionOrder=0;%0.5%Order for fractional derivative-based motion estimation
if rec.Par.Mine.Modal==10;Alg.parU.fractionOrder=1;end%Force fraction order to be 1 for DWI
Alg.parU.iterMask=1;%To iterate the masking process
Alg.parU.Lambda=1;%To regularize multislice reconstructions
Alg.parU.interpolTime=[0 0];%To interpolate using shearing by temporal order of multiband excitations (distortion / motion)

%DEFAULT PARAMETERS FOR SIGNAL RECOVERY
Alg.parR.illustration=0;%Whether to show information at the different steps of the algorithm
Alg.parR.Gamma=0.2:0.05:0.95;%Random matrix aspect ratio, vector of values interpreted as candidates for patch size estimations
Alg.parR.Subsampling=[2 2 2];%Subsampling factor for patch construction
Alg.parR.DrawFact=3;%Factor to multiply the subsampling factor to sweep the image space to estimate patch sizes
Alg.parR.NoiEstMeth='None';%Noise level estimation method. One of the following: 'None' / 'Exp1' / 'Exp2' / 'Medi'
Alg.parR.ShrinkMeth='Frob';%Shrinkage method. One of the following: 'Hard' / 'Soft' / 'Frob' / 'Oper' / 'Nucl' / 'Exp1' / 'Exp2'
Alg.parR.DirInd=3;%Direction of independent data to accelerate computations of ESD
Alg.parR.HalfScanCorrection=2;%Type of correction for half scan filter. 0-> no correction, 1-> by regridding, 2-> by ESD estimation ramped,3-> by ESD estimation zero-filled,4-> by regridding in pairs of PEs,5-> by zero-filled in singles
Alg.parR.GFactorCorrection=2;%Type of correction for noise profiles after reconstruction: 0-> no correction, 1-> by spatial normalization and local covariance, 2-> by local full covariance
Alg.parR.LinearPhaseCorrection=1;%Not to run correction / Run linear phase correction
Alg.parR.FilteredPhaseCorrection=3;%0->None, 1-> low pass, 2-> shearlet, 3-> full
Alg.parR.WeightAssemb='Gauss';%Type of window weighting for patch assembling. One of the following: 'Gauss' / 'Unifo' / 'Invva'
Alg.parR.PatchSimilar=2;%Similarity metric to build the patches. 2 for Euclidean / 1 for Manhattan
Alg.parR.ExtractROI=1;%To extract ROI in PE direction
Alg.parR.ESDMeth=0;%0->Use simulation / 1-> Use SPECTRODE / 2-> Use MIXANDMIX / 3-> Force computation of simulated spectra (to get AMSE estimates when using noise estimates instead of noise modeling)
Alg.parR.UseComplexData=1;%To use complex data denoising
Alg.parR.ESDTol=1e-4;%SPECTRODE / MIXANDMIX algorithm tolerance
Alg.parR.NR=1;%Number of realizations of random matrix simulations
Alg.parR.useRatioAMSE=1;%To compute the AMSE ratio (1), AMSE (0) or both (2)
Alg.parR.Verbosity=0;%Use 2 to plot some timing information

%ASSIGNMENT
rec.Alg=Alg;

%ERROR CONTROL
if min(rec.Par.Scan.AcqVoxelSize(:))<Alg.parXT.discardHighRes && rec.Par.Mine.Modal==7;fprintf('Voxel resolution is%s and finest supported resolution has been fixed to %.2f. Change corresponding parameter try the reconstruction\n',sprintf(' %.2f',min(rec.Par.Scan.AcqVoxelSize,[],1)),Alg.parXT.discardHighRes);rec.Fail=1;return;end
