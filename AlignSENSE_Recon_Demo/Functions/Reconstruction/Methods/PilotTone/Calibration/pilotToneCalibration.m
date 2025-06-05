
function [A, condA] = pilotToneCalibration(p, T, calibrationModel, calibrationOffset, We)

if nargin<3 || isempty(calibrationModel); calibrationModel='backward';end %Modelling motion as a linaer combination of PT signal T = Ab*p (backward) or PT as p=A*T (forward = physical setup)
if nargin<4 || isempty(calibrationOffset); calibrationOffset=1;end 
if nargin<5 || isempty(We); We=[];end %Weights for a weighted calibration

calibToDegrees=1;%Test for rot-ran imbalance

%%% CHECK DATA CONSISTENCY
NP = size(p);
NCha=NP(1);
NT = size(T);
NStates = NP(2);
nD = length(NT);
assert(NT(nD-1)==NP(2) , 'pilotToneCalibration:: Number of states in PT signal and motion traces is not the same.');

isComplex = ~isreal(p);%Whether to model PT signal and calibration matrix as complex valued.
T=permute(T, [6 5 1:4]);%Now NMotparam x NSamples
if calibToDegrees;T = dynInd(T,4:6,1,convertRotation( dynInd(T,4:6,1),'rad','deg') );end%In degrees

if ~ isempty(We)%Weighted Least Squares
    We = reshape(We, [ 1 NStates ]);
    We = single(We);%In case We are logical values
    p = bsxfun(@times,p, sqrt(We) );
    T = bsxfun(@times,T, sqrt(We) ); 
end
%%% CALIBRATE USING LEAST SQUARES
if strcmp(calibrationModel, 'forward') %p = A*T
    if calibrationOffset; T = cat(1, T, onesL(dynInd(T,1,1)));end%Offset in motion parameters
    
    A = p*pinv(T);% Same as A = p / T or A = p*inv(T);

elseif strcmp(calibrationModel, 'backward') %T = A*p
    if calibrationOffset; p = cat(1, p, onesL(dynInd(p,1,1)) );end%Offset in PT signal (along channels)
   
    A = T*pinv(p); % Same as A = T / p or A = p*inv(T);      
else
    error('pilotTonePrediction:: %s not recognised as a calibrationModel.',calibrationModel)
end

if ~isComplex; A=real(A);end

%%% REPORT CALIBRATION CHARACTERISTICS
condA=cond(A);