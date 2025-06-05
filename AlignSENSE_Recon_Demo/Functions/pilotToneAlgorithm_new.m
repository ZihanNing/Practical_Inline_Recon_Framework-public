
function rec = pilotToneAlgorithm_new(rec)

%PILOTTONEALGORITHM   Sets the default parameters for the Pilot Tone model used in the reconstruction
%   REC=PILOTTONEALGORITHM(REC)
%   * REC is a reconstruction structure without the algorithmic parameters. 
%   At this stage it may contain the naming information (.Names), the 
%   status of the reconstruction (.Fail), the .lab information (.Par), the 
%   fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   ** REC is a reconstruction structure with filled algorithmic information (.Alg)
%
%   Yannick Brackenier - 2023-04-05

%%% PILOT TONE DATA EXTRACTED FROM ACQUIRED K-SPACE
if ~isfield(rec,'PT');rec.PT=[];end
if ~isfield(rec.PT,'yProjRO');rec.PT.yProjRO=[];end %The projection of the PT signal along the readout direction. Used to detect the PT signal.
if ~isfield(rec.PT,'pSliceImage');rec.PT.pSliceImage=[];rec.parXT.usePT=0;end%The sice of PT signal in the image domain. If no PT signal detected, disable its usage.
if ~isfield(rec.PT,'pTimeTest');rec.PT.pTimeTest=[];end%Stored temporal PT signal from data extraction. Only used for comparison.
if ~isfield(rec.PT,'factorFOV');rec.PT.factorFOV=[];end%Factor indicating where the PT signal is situated in the over-sampled FOV. 0 for centre FOV and 1 for the edge of the over-sampled FOV.
if ~isfield(rec.PT,'idxRO');rec.PT.idxRO=[];end%Index in logical units of where the PT signal is situated in the over-sampled FOV.
if ~isfield(rec.PT,'idxMB');rec.PT.idxMB=[];end%Index where the peak PT signal is situated in the extracted Multi-Band PT singal.
if ~isfield(rec.PT,'offsetMBHz');rec.PT.offsetMBHz=[];end%Bandwidth of the extracted multiband information in Hz.
if ~isfield(rec.PT,'offsetMBIdx');rec.PT.offsetMBIdx=[];end%Multi-Band region expressed in logical coordinates.
if ~isfield(rec.PT,'FWHMHz');rec.PT.offsetMBIdx=[];end%FWHM of the multiband information in Hz.
if ~isfield(rec.PT,'isAveragedMB');rec.PT.isAveragedMB=0;end%Whether the MB information is averaged when extracting the data.
if ~isfield(rec.PT,'isRelativePhase');rec.PT.isRelativePhase=0;end%Whether phase of the PT signal is referened to a certain coil.
if ~isfield(rec.PT,'relativePhaseIdx');rec.PT.relativePhaseIdx=[];end%Coil index for the phase referencing.
if ~isfield(rec.PT,'transmitFreq');rec.PT.transmitFreq=[];end%Calculated transmit frequency of the injected PT signal (not used).

%%% ALGORITHM
%Optimisation scheme
if ~isfield(rec.Alg.parXT,'PT');rec.Alg.parXT.PT=[];end
if ~isfield(rec.Alg.parXT.PT,'usePT');rec.Alg.parXT.PT.usePT=0;end %Don't use PT (0), use PT for alternating estimation (1), use PT at end to predict next step (2)
if ~isfield(rec.Alg.parXT.PT,'subDiv');rec.Alg.parXT.PT.subDiv=0;end %Sub-dividing shots after level
if ~isfield(rec.Alg.parXT.PT,'nItActivate');rec.Alg.parXT.PT.nItActivate=10;end

%Part of Pilot Tone signal to use
if ~isfield(rec.Alg.parXT.PT,'signalUsage');rec.Alg.parXT.PT.signalUsage=[];end
if ~isfield(rec.Alg.parXT.PT.signalUsage,'SupportIS');rec.Alg.parXT.PT.signalUsage.SupportIS=[];end%Percentage of elements to use for PT signalin the Inferior-Superior direction
if ~isfield(rec.Alg.parXT.PT.signalUsage,'combineMB');rec.Alg.parXT.PT.signalUsage.combineMB=1;end%Take centre peak (1) / Average (2) if data extraction allow it (see rec.PT.isAveragedMB)
if ~isfield(rec.Alg.parXT.PT.signalUsage,'referencePhase');rec.Alg.parXT.PT.signalUsage.referencePhase=0;end%Which channel to reference the phase of the PT signal to. If zero, not applied.
if ~isfield(rec.Alg.parXT.PT.signalUsage,'orderPreProcessing');rec.Alg.parXT.PT.signalUsage.orderPreProcessing=[];end%(1) de-mean / (2) de-mean phase / (3) normalise / (4) whiten / (5) unit variance / (6) real-imag concatenation / (7) magnitude
if ~isfield(rec.Alg.parXT.PT.signalUsage,'eigTh');rec.Alg.parXT.PT.signalUsage.eigTh=onesL(rec.Alg.resPyr);end%Threshold for SVD thresholding. If negative, the numbers to keep.
if ~isfield(rec.Alg.parXT.PT.signalUsage,'useRealImagSVD');rec.Alg.parXT.PT.signalUsage.useRealImagSVD=1;end%To perform the SVD on the concatenated real/imaginary components.
if ~isfield(rec.Alg.parXT.PT.signalUsage,'includeSVScaling');rec.Alg.parXT.PT.signalUsage.includeSVScaling=1;end%To include the SV scaling in the vectors used as a new basis.

%Signal filtering
if ~isfield(rec.Alg.parXT.PT,'PTFilter');rec.Alg.parXT.PT.PTFilter=[];end
%if ~isfield(rec.Alg.parXT.PT.PTFilter,'Type');rec.Alg.parXT.PT.PTFilter.Type='median';end%Type of filtering of the temporal PT signal.
if ~isfield(rec.Alg.parXT.PT.PTFilter,'medFiltKernelWidthms');rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthms=0;end%Size of the median filter kernelin ms used for removing outliers in the temporal PT signal.
rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthidx = max(1,rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthms/rec.Par.Labels.RepetitionTime);%Converted to logical units.
rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthidx = rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthidx - mod(rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthidx,2)+1;%Enforce it to be odd
if ~isfield(rec.Alg.parXT.PT.PTFilter,'golayFiltKernelWidthms');rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthms=0;end%Size of the Savitzky-Golay filter kernel in ms used for removing outliers in the temporal PT signal.
rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthidx = max(1,rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthms/rec.Par.Labels.RepetitionTime);%Converted to logical units.
rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthidx = rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthidx - mod(rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthidx,2)+1;%Enforce it to be odd

%Outlier detection
if ~isfield(rec.Alg.parXT.PT,'weightedRecon');rec.Alg.parXT.PT.weightedRecon=[];end
if ~isfield(rec.Alg.parXT.PT.weightedRecon,'Flag');rec.Alg.parXT.PT.weightedRecon.Flag=0;end%Whether to use weights from PT instead of residual-based weights in last image reconstruction (1) or in every iteration(2)
if ~isfield(rec.Alg.parXT.PT.weightedRecon,'Method');rec.Alg.parXT.PT.weightedRecon.Method='outlWePT';end%Which weights to use for the above. See getPTWeighting.m.

%Calibration
if ~isfield(rec.Alg.parXT.PT,'Calibration');rec.Alg.parXT.PT.Calibration=[];end
if ~isfield(rec.Alg.parXT.PT.Calibration,'externalFit');rec.Alg.parXT.PT.Calibration.externalFit=[];end
if ~isfield(rec.Alg.parXT.PT.Calibration.externalFit,'Flag');rec.Alg.parXT.PT.Calibration.externalFit.Flag=0;end
if ~isfield(rec.Alg.parXT.PT.Calibration.externalFit,'Ab');rec.Alg.parXT.PT.Calibration.externalFit.Ab=[];end
if ~isfield(rec.Alg.parXT.PT.Calibration.externalFit,'MT');rec.Alg.parXT.PT.Calibration.externalFit.MT=[];end
if ~isfield(rec.Alg.parXT.PT.Calibration.externalFit,'NX');rec.Alg.parXT.PT.Calibration.externalFit.NX=[];end
if ~isfield(rec.Alg.parXT.PT.Calibration.externalFit,'forInit');rec.Alg.parXT.PT.Calibration.externalFit.forInit=0;end
if ~isfield(rec.Alg.parXT.PT.Calibration,'implicitCalibration');rec.Alg.parXT.PT.Calibration.implicitCalibration=1;end %To calibrate the raw k-space (1) or the surrogate motion parameters (0)
if ~isfield(rec.Alg.parXT.PT.Calibration,'implicitCalibrationDelay');rec.Alg.parXT.PT.Calibration.implicitCalibrationDelay=0;end %
if ~isfield(rec.Alg.parXT.PT.Calibration,'calibrationModel');rec.Alg.parXT.PT.Calibration.calibrationModel='backward';end%Whther to use the forward or backward calibration model.
if ~isfield(rec.Alg.parXT.PT.Calibration,'calibrationOffset');rec.Alg.parXT.PT.Calibration.calibrationOffset=1;end%Whether to include an offset in the calibration matrix.
if ~isfield(rec.Alg.parXT.PT.Calibration,'calibrationOrder');rec.Alg.parXT.PT.Calibration.calibrationOrder=1;end%Order of the calibration model.
if ~isfield(rec.Alg.parXT.PT.Calibration,'useRASMotion');rec.Alg.parXT.PT.Calibration.useRASMotion=0;end%To use the transformation expressed in real world coordinates w.r.t. isocentre (1). If set to 2, the logaritmic parameters of the latter are used for calibration.
if ~isfield(rec.Alg.parXT.PT.Calibration,'weightedCalibrationFlag');rec.Alg.parXT.PT.Calibration.weightedCalibrationFlag=0;end %To have a weighted calibration.
if ~isfield(rec.Alg.parXT.PT.Calibration,'weightedCalibrationMethod');rec.Alg.parXT.PT.Calibration.weightedCalibrationMethod='outlWePT';end%Which weights to use for the above. See getPTWeighting.m.

%Clustering
if ~isfield(rec.Alg.parXT.PT,'profileClustering');rec.Alg.parXT.PT.profileClustering=[];end
if ~isfield(rec.Alg.parXT.PT.profileClustering,'Flag');rec.Alg.parXT.PT.profileClustering.Flag=zerosL(rec.Alg.resPyr);end %Whether to use the PT-based binning.
if ~isfield(rec.Alg.parXT.PT.profileClustering,'numClusters');rec.Alg.parXT.PT.profileClustering.numClusters=300;end %The number of clusters to use for profile binning
if ~isfield(rec.Alg.parXT.PT.profileClustering,'traLimX');rec.Alg.parXT.PT.profileClustering.traLimX=[];end %Threshold in within-cluster-distance for determining numClusters (over-writes the parameters numClusters)
if ~isfield(rec.Alg.parXT.PT.profileClustering,'traLimXPerc'); rec.Alg.parXT.PT.profileClustering.traLimXPerc = 1;end %Percentage of TR's that need to fullfil traLimX criteria to find appropriate number of clusters.
if ~isfield(rec.Alg.parXT.PT.profileClustering,'clusterMotion');rec.Alg.parXT.PT.profileClustering.clusterMotion=1;end %Whether to cluster motion parameters (1) or the raw PT signal (0)
if ~isfield(rec.Alg.parXT.PT.profileClustering,'equalBinSize');rec.Alg.parXT.PT.profileClustering.equalBinSize=1;end %Whether to enforce clusters to have the same number of profiles

            
%%% CHECKS
assert(multDimSum(single([any(rec.Alg.parXT.PT.profileClustering.Flag)  any(rec.Alg.parXT.PT.subDiv)]))<=1,'Cannot do subdivision and clustering in same reconstruction.');
