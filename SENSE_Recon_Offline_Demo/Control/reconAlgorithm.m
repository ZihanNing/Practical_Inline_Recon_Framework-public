function rec=reconAlgorihm(rec)

%RECONALGORITHM   Parameters for a reconstruction
%   REC=RECONALGORITHM(REC)
%   * REC is a recon structure
%   ** REC is a recon structure
%

%To control which pipeline blocks to run
rec.Pip.estimateSensitivities=1;%To estimate sensitivities
rec.Pip.invertData=1;%To invert data
rec.Pip.readFromFile=[1 1 1];%To read previous reconstructions of body/sensitivities/built-in sensitivities from file

%To control results
rec.Alg.unringingStrength='FilterHard';%Alternatives for unringing, 'None' for preserving resolution, 'FilterModerate' for moderate unringing, 'FilterHard' for hard filtering, 'Subvoxels' for subvoxel-shift method (non-linear)
if strcmp(rec.Alg.unringingStrength,'None');rec.Alg.gibbsRinging=1;
elseif strcmp(rec.Alg.unringingStrength,'FilterModerate');rec.Alg.gibbsRinging=[0.2 0.2];
elseif strcmp(rec.Alg.unringingStrength,'FilterHard');rec.Alg.gibbsRinging=[0.4 0.4];
elseif strcmp(rec.Alg.unringingStrength,'Subvoxels');rec.Alg.gibbsRinging=[];
else error('%d unringing method not contemplated',rec.Alg.unringingStrength)
end 
rec.Alg.writeRaw=0;%To write raw data for KRGDSVR
rec.Alg.correctMotion=1;%To correct for motion

%rec.Alg.supportReadoutS=[0.25 0.75];%To extract a readout area to accelerate sensitivity estimation
%rec.Alg.supportReadout=[0.25 0.75];%To extract a readout area to accelerate reconstructions
rec.Alg.supportReadoutS=[];%To extract a readout area to accelerate sensitivity estimation
rec.Alg.supportReadout=[];%To extract a readout area to accelerate reconstructions

rec.Alg.oversamplingFactors=ones(1,3);%Oversampling factors
rec.Alg.removeOversamplingReadout=1;%To remove oversampling in readout when reading the data
rec.Alg.bodyCoilEstimationMethod='BCC';%One of the following: RSoS (Root Sum of Squares), CC (Coil Compression), BCC (Block Coil Compression), CCPh (Coil Compression Phase)
%rec.Alg.bodyCoilEstimationMethod='RSoS';%One of the following: RSoS (Root Sum of Squares), CC (Coil Compression), BCC (Block Coil Compression), CCPh (Coil Compression Phase)
rec.Alg.sensitivityEstimationMethod='ESPIRIT';%One of the following: Filter (Quick estimate), Standard (as in Philips suite), ESPIRIT
%rec.Alg.sensitivityEstimationMethod='Filter';%One of the following: Filter (Quick estimate), Standard (as in Philips suite), ESPIRIT
rec.Alg.maskThS=0.0001;%Threshold of the body coil intensities for mask extraction, if larger less likely occurrence of background artifacts, if too large, boundaries may be eroded
%rec.Alg.maskThS=0.001;%Threshold of the body coil intensities for mask extraction, if larger less likely occurrence of background artifacts, if too large, boundaries may be eroded
%MAY MAKE SENSE TO MAKE THIS LARGER OR IMPROVE

rec.Alg.parS.maskTh=0.2;%Threshold of the body coil intensities for mask extraction%It was 1
rec.Alg.parS.lambda=1;%This parameter has been constantly changing, toddlers have been processed with it set to 1, but it was introducing noise for neonates, 2 may suffice, but 5 gives a bit of margin%20 previously
rec.Alg.parS.Otsu=[];%[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)

rec.Alg.ThreeDEspirit=1;%To use 3D ESPIRIT
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=2*ones(1,3);else rec.Alg.parE.NC=2*ones(1,2);end%Resolution (mm) of calibration area to compute compression
if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=1.5*ones(1,3);else rec.Alg.parE.NC=1.5*ones(1,2);end%Resolution (mm) of calibration area to compute compression
if rec.Alg.ThreeDEspirit;rec.Alg.parE.eigSc=[0.98 0.3];else rec.Alg.parE.eigSc=[0.9 0.3];end%Possible cut-off for the eigenmaps in soft SENSE, used the first for mask extraction by thresholding the eigenmaps%Default was [0.85 0.3]
if rec.Alg.ThreeDEspirit;rec.Alg.parE.Ksph=200;else rec.Alg.parE.Ksph=50;end%Number of points for spherical calibration area, unless 0, it overrides other K's%Default was 200
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.125;else rec.Alg.parE.factScreePoint=0.25;end%Factor over threshold computed with scree point
if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.5;else rec.Alg.parE.factScreePoint=1;end%Factor over threshold computed with scree point
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=2;else rec.Alg.parE.factScreePoint=4;end%Factor over threshold computed with scree point
rec.Alg.parE.mirr=[0 0 0];%Whether to mirror along a given dimension%Default was [8 8 8]
rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
%rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
rec.Alg.parE.factorBody=1;%Factor for virtual body coil in ESPIRIT (it was 1)
%if strcmp(field,'refscan')
    %rph.Alg.parE.mirr=[0 0 2];
    %rph.Alg.parE.factorBody=10;
%end

rec.Alg.useMasking=0;%To use masking in reconstruction---BUG: does not work well when reading previous sensitivity estimation
%rec.Alg.useBuiltInCalibration=2;%To use built in calibration, 1 uses only phases, 2 uses sensitivities
rec.Alg.useBuiltInCalibration=1;%To use built in calibration, 1 uses only phases, 2 uses sensitivities
rec.Alg.averagePolaritySensitivities=0;%0->not average at all, 1-> average after calculation, 2-> average before calculation

%Exceptions
% if contains(rec.Nam.caseIn,'DISORDER_INVIVO_7T_02');rec.Alg.supportReadoutS=[0.4 0.85];rec.Alg.supportReadout=[0.45 0.8];end
% if contains(rec.Nam.caseIn,'DISORDER_03') || contains(rec.Nam.caseIn,'DISORDER_05');rec.Alg.supportReadoutS=[0.275 0.725];rec.Alg.supportReadout=[0.325 0.675];end
% 
% if contains(rec.Nam.caseIn,'24GiulioPET');rec.Alg.ovesamplingFactors(3)=0.8;end

rec.Alg.maxNumberRepeats=10;

%%HEREHEREHERE---WE NEED TO INTRODUCE THE MASK, SOFT MASKING ETC