

function [ WePT, outlWePT, MADShot, stdShot ] = getPTWeighting( pTime , mSt , NStates, percMAD, ZScoreOutl, useSVD, nComponents, lossFunctionInfo)

%GETPTWEIGHTING
%   * PTIME the pilot tone signal over time with the channels in the rows and time along the columns.
%   * MST the indices defining the grouping.
%   * {NSTATES} the motion states. Note this is only used when there is an elliptical shutter (only relevant for DISORDER).
%   * {PERCMAD} the percentile for relating the MAD estimator to the standard deviation of a Normal distribution.
%   * {ZSCOREOUTL} Z-score above which shots are detected as outliers.
%   * {USESVD} whether to use the Singular vectors to detect the outliers on.
%   * {NCOMPONENTS} the number of Singular vectors to remain for the analysis.
%   * {LOSSFUNCTIONINFO} a structure containing information about the loss  function for outlier detection
%   ** WEPT the soft weights from PT
%   ** OUTLWEPT indicating which state is an outlier detected from PT
%

%%% Coeffients to tweak
%   percMAD: bigger creates a smaller sigma for the Normal distribution
%   ZScoreOutl: making bigger will result in less outliers
%   lossFunctionInfo.delta: bigger will increase weights

if nargin < 2 || isempty(mSt); warning('getPTWeighting:: mSt not provided. Assumed to be a single shot.');mSt=ones([1 pTime]);end %List containing the motion states of every sample
if nargin < 3 || isempty(NStates); NStates=max(mSt);end %Number of states
if nargin < 4 || isempty(percMAD); percMAD=.5;end %Percentile for estimating sigma based on MAD estimator
if nargin < 5 || isempty(ZScoreOutl); ZScoreOutl=1.5;end %ZScore above which to categorise shots as outloers based on their standard deviation
if nargin < 6 || isempty(useSVD); useSVD=1;end %Whether to use the singular components to do the outlier detection or the raw coil signals.
if nargin < 7 || isempty(nComponents); nComponents=1;end %Number of singular components to keep
if nargin < 8 || isempty(lossFunctionInfo) %Loss function to compute weights
    lossFunctionInfo=[];
    lossFunctionInfo.Type='FairPotential';%Good approximation for L1 which is differentiable: https://web.eecs.umich.edu/~fessler/papers/files/talk/16/um-sjtu-ji.pdf
    lossFunctionInfo.delta=15;%Might need to fine-tune
end 
assert(nComponents==1,'getPTWeighting:: Not implemented with nComponents>1.');

N = size(pTime);
NCha = N(1);
NSamples = N(2);

%% Take SVD and extract first nComponents
if useSVD %Use singular vectors
    [~, ~, v] = svd(pTime,"econ");
    pTime = v(:,1:min(nComponents,NSamples)).';
else %Use coil with biggest coil amplitude
    tt= multDimMea( abs(pTime),2);
    [ ~, idxMax ] = max (tt,[],1);
    pTime = dynInd(pTime,idxMax, 1);
end

%Take absolute value
pTime=abs(pTime); %Now only 1 (virtual) channel

%% Shot handling
%Group the shots based on stateSample indices       
stateSample= permute( mSt(mSt<=NStates),[2 1]);%Exclude elliptical shutter
samplesShot = regionprops(stateSample,pTime,'PixelValues');
samplesShot={samplesShot.PixelValues};

%Compute MAD and standard deviation for every shot
MADShot = cellfun(@(x) mad(x(:),1), samplesShot);
stdShot = cellfun(@(x) std(x(:)), samplesShot);

%Get an global MAD estimator for the motion trace
MADGlobal = median(MADShot); %Can also take mean. Lucilio anticipated that it won't matter that much.

%Use global MAD estimator to infer the standard deviation for a distribution of pilot tone sgnal (that does not contain variations due to motion)
k = 1/(sqrt(2)*erfinv(percMAD));%Assume a normal distribution: https://en.wikipedia.org/wiki/Median_absolute_deviation
sigmaDistribution = k * MADGlobal ;

%% Normalse standard deviation of every shot w.r.t. the distribution
ZScore = stdShot./sigmaDistribution;

%% Copmute weights bad on the the method Lucilio used in the IRLS for robust regression
if strcmp(lossFunctionInfo.Type,'FairPotential') 
    fDeriv = @(x) ( x / (1+abs(x/lossFunctionInfo.delta)) );
else
   error('getPTWeighting:: Loss function %s not implemented yet.',lossFunctionInfo.Type) 
end
WePT = fDeriv(ZScore)./ZScore;
WePT = min(WePT,1);%Set supremum to 1

%% Detect outliers 
outlWePT = ZScore>ZScoreOutl;






