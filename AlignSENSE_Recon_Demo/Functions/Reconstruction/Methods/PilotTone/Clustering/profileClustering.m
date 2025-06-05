

function [stateSample, PT, parPT, TClustered, dClusters, dSamples, dSamplesOwnCluster, dTranClust, dRotClust] = profileClustering(p,PT,parPT,MS,stateSampleInit,Tinit,deb)

%PROFILECLUSTERING groups readout profiles into clusters that can be used for a final reconstruction. Clustering is based on pilot tone signal
%   provided per readout.
%
%   [STATESAMPLE, PT, PARPT, TCLUSTERED, DCLUSTERS, DSAMLES, DSAMPLESOWNCLUSTER] = PROFILECLUSTERING(P,PT,PARPT,MS,DEB)
%   * P is the Pilot Tone signal with dimensions NChannelsxNSamples.
%   * PT is the structure with all the PT data stored. See pilotToneAlgorithm.m for an overview.
%   * {PARPT} is the structure with all the PT parameters stored that inform clustering. See pilotToneAlgorithm.m for an overview.
%   * {MS} is the resolution of the native resolution (independent of the resolution level this is called from).
%   * {DEB} a flag to report intermediate information.
%   ** STATESAMPLE returns the cluster index for every sample.
%   ** PT is the structure with all the PT data stored.
%   ** PARPT is the structure with all the PT parameters stored. 
%   ** TCLUSTERED are the motion parameters for the clusters.
%   ** DCLUSTERS is the averaged distance within each cluster. Distance is not averaged across motion parameters.
%   ** DSAMPLES are the distances for every sample w.r.t every cluster. Distance is not averaged across motion parameters.
%   ** DSAMPLESOWNCLUSTER are the distances for every sample w.r.t its assigned cluster. Distance is not averaged across motion parameters.
%   ** DTRANCLUST are the averaged translation distances per cluster (averaged over all 3 translation dimensions).
%   ** DROTCLUST are the averaged rotation distances per cluster (averaged over all 3 rotation dimensions).
%

if nargin<2 || isempty(PT); error('profileClustering:: Not enough input arguments.');end
if nargin<3 || isempty(parPT); parPT=[];parPT.profileClustering = []; end
if nargin<4 || isempty(MS); MS=ones([1 3]);end
if nargin<7 || isempty(deb); deb=0;end

%%% SET PARAMETERS
if ~isfield(parPT.profileClustering,'Method'); parPT.profileClustering.Method = 'KMeans';end
if ~isfield(parPT.profileClustering,'numClusters'); parPT.profileClustering.numClusters = 100;end
if ~isfield(parPT.profileClustering,'clusterMotion'); parPT.profileClustering.clusterMotion = 1;end %To cluster motion parameters vs. orignal pilot tone signal
if ~isfield(parPT.profileClustering,'traLimX'); parPT.profileClustering.traLimX = [];end 
if ~isfield(parPT.profileClustering,'traLimXPerc'); parPT.profileClustering.traLimXPerc = 1;end %Percentage of TR's that need to fullfil traLimX criteria to find appropriate number of clusters.
if length(parPT.profileClustering.traLimX)==1; parPT.profileClustering.traLimX=parPT.profileClustering.traLimX*[1 1];end

dimS=5;%Dimension where motion states are stored - hard-coded
dimM=6;%Dimension where motion parameters are stored - hard-coded

%%% PREDICT MOTION ON A TR LEVEL
if strcmp(parPT.Calibration.calibrationModel,'forward')
    T = pilotTonePrediction (PT.Af, p, 'PT2T','forward', PT.NChaPTRec, 6);
elseif strcmp(parPT.Calibration.calibrationModel,'backward')
    T = pilotTonePrediction (PT.Ab, p, 'PT2T','backward', PT.NChaPTRec, 6);
else
    error('profileClustering:: One calibration model must be implemented to call this function.');
end

if ~parPT.profileClustering.clusterMotion
   warning('profileClustering:: Clustering based on raw PT signal. Only returned stateSample and distances are meaningfull.'); 
end

%%% PREPARE DATA
if parPT.profileClustering.clusterMotion
    perm = [dimS dimM setdiff(1:numDims(T),[dimS, dimM])];
    X = permute(T, perm);%Re-order data so that rows are the data points (TR) and columns the variables (motion parameters)
    X(:,1:3) = X(:,1:3).*MS;%In units of mm
    X(:,4:6) = convertRotation(X(:,4:6),'rad','deg');%In units of deg
    
    XClusteredInit = permute(Tinit, perm);%Re-order data so that rows are the data points (TR) and columns the variables (motion parameters)
    XClusteredInit(:,1:3) = XClusteredInit(:,1:3).*MS;%In units of mm
    XClusteredInit(:,4:6) = convertRotation(XClusteredInit(:,4:6),'rad','deg');%In units of deg
else
    X=p.';
end

%%% REPORT INITIAL GROUPING CHARACTERISTICS
if parPT.profileClustering.clusterMotion
    [dSamp,dTranSamp,dRotSamp,dTranClust,dRotClust] = computeDispersion(X,XClusteredInit,stateSampleInit);
    if deb
        fprintf('Initial dispersion from DISORDER sample grouping: \n');
        fprintf('   Mean translation dispersion across clusters = %.2f mm\n',multDimMea(dTranClust) );
        fprintf('   Mean rotation dispersion across cluster = %.2f deg\n',multDimMea(dRotClust) );
    end       
end

%%% APPLY CLUSTERING
K = parPT.profileClustering.numClusters;
kInit = 100; KUp = 50; KMax = 200;%Hard-coded parameters to grid-search for # of cluster
options = statset('UseParallel',1);%,'Streams',stream);
distMetric = 'cityblock';%L1 distance. Options: 'sqeuclidean' (default) | 'cityblock' | 'cosine' | 'correlation' | 'hamming'
numIt=100000;
nReplicas = 10;
warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
    
tSta = tic;
if ~isempty(parPT.profileClustering.traLimX) && parPT.profileClustering.clusterMotion%Over-write the number of clusters and see look for a K that satifies the required within distance
    
    for KTemp = kInit:KUp:KMax
        %Run KMeans
        [stateSample, XClustered, dClusters, dSamples] = kmeans(X, KTemp,'Options', options,'Distance', distMetric,'MaxIter',numIt,'Replicates',nReplicas);
        %Compute dispersion within clusters   
        [dSamp,dTranSamp,dRotSamp,dTranClust,dRotClust] = computeDispersion(X,XClustered,stateSample);
        %Report
        if deb
            fprintf('Trying clustering with K = %d \n',KTemp);
            fprintf('   Mean translation dispersion across clusters = %.2f mm\n',multDimMea(dTranClust) );
            fprintf('   Mean rotation dispersion across cluster = %.2f deg\n',multDimMea(dRotClust) );
        end
        %Check criteria
        if multDimSum( dTranClust < parPT.profileClustering.traLimX(1) & dRotClust < parPT.profileClustering.traLimX(2) ) >= round(parPT.profileClustering.traLimXPerc*length(dTranClust))
            K = KTemp; 
            break;
        elseif KMax-KTemp < KUp
            fprintf('No K found that fullfils traLimX criteria. K set to maximum value = %d.\n',KTemp);
            K = KTemp;
        end
    end
else
    [stateSample, XClustered, dClusters, dSamples] = kmeans(X, K,'Options',options,'Distance',distMetric,'MaxIter',numIt,'Replicates',nReplicas);%Run KMeans
end
parPT.profileClustering.numClusters=K;%Assign the number of clusters used in the end
tEnd = toc(tSta);
warning('on','MATLAB:handle_graphics:exceptions:SceneNode');

%%% COMPUTE FINAL DISPERSIONS AND REPORT
if parPT.profileClustering.clusterMotion
    [dSamp,dTranSamp,dRotSamp,dTranClust,dRotClust] = computeDispersion(X,XClustered,stateSample);
    if deb
        fprintf('Final dispersion from optimised sample clustering: \n');
        fprintf('   Mean translation dispersion across clusters = %.2f mm\n',multDimMea(dTranClust) );
        fprintf('   Mean rotation dispersion across cluster = %.2f deg\n',multDimMea(dRotClust) );
    end       
end    

dClusters = dClusters./ accumarray(stateSample, onesL(stateSample));%Make it mean distance for every sample
dSamplesOwnCluster = zeros( [size(dSamples,1) ,1], 'like', dSamples);
for k=1:K; idx = stateSample==k; dSamplesOwnCluster(idx) = dSamples(idx,k);end

%%%ARRANGE DATA IN FINAL FORM FOR RECONSTRUCTION
if parPT.profileClustering.clusterMotion
    TClustered = XClustered;
    %Restore right units
    TClustered(:,1:3) = TClustered(:,1:3)./MS;
    TClustered(:,4:6) = convertRotation(TClustered(:,4:6),'deg','rad');
    
    TClustered = ipermute(single(TClustered),perm);
else
    if strcmp(parPT.Calibration.calibrationModel,'forward')
        TClustered = pilotTonePrediction (PT.Af, XClustered.', 'PT2T','forward', PT.NChaPTRec, 6);
    elseif strcmp(parPT.Calibration.calibrationModel,'backward')
        TClustered = pilotTonePrediction (PT.Ab, XClustered.', 'PT2T','backward', PT.NChaPTRec, 6);
    end

    dSamp = T; for k=1:K;idx=stateSample==k; dSamp=dynInd(dSamp,idx,dimS, dynInd(dSamp,idx,dimS)- dynInd(TClustered,k,dimS) );end
    dTranSamp = multDimMea(abs(dynInd(dSamp,1:3,dimM)),dimM);%Over dims - in mm
    dRotSamp = multDimMea(abs(dynInd(dSamp,4:6,dimM)),dimM);%Over dims - in rad

    dTranClust = zerosL(dClusters); dRotClust = zerosL(dClusters);
    for k=1:K; idx = stateSample==k; dTranClust(k) = multDimMea( dynInd(dTranSamp,idx,dimS));dRotClust(k) = multDimMea( dynInd(dRotSamp,idx,dimS));end
end      
        
TClustered = gather(TClustered);
stateSample = single(stateSample.');

%%% REPORT 
if deb
    fprintf('Grouping TR''s into clusters:\n');
    fprintf('   Clustering based on motion parameters = %d.\n',parPT.profileClustering.clusterMotion);
    fprintf('   Number of samples = %d.\n',size(X,1));
    fprintf('   Number of clusters = %d.\n',K);
    fprintf('   Time computing = %.1f seconds.\n',tEnd);
    fprintf('   Averaged within-cluster distance per parameter = %.2f deg-mm.\n',multDimMea(dClusters/size(X,2)));
end

end


%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dSamp,dTranSamp,dRotSamp,dTranClust,dRotClust] = computeDispersion(X,XClustered,stateSample)
    %Get the number of clusters
    K = size(XClustered,1);
    %Compute the distance per sample
    dSamp = X(:,1:6); for k=1:K;idx=stateSample==k; dSamp=dynInd(dSamp,idx,1, dynInd(dSamp,idx,1)- XClustered(k,1:6));end
    %Compute distance for translations and rotations separately
    dTranSamp = multDimMea(abs(dSamp(:,1:3)),2);%Over dims - in mm
    dRotSamp = multDimMea(abs(dSamp(:,4:6)),2);%Over dims - in deg
    %Look at distance within cluster
    dTranClust = zeros([1 K],'like',dSamp);
    dRotClust = zeros([1 K],'like',dSamp);
    for k=1:K
        idx = stateSample==k;
        dTranClust(k) = multDimMea( dynInd(dTranSamp,idx,1));
        dRotClust(k) = multDimMea( dynInd(dRotSamp,idx,1));
    end
end