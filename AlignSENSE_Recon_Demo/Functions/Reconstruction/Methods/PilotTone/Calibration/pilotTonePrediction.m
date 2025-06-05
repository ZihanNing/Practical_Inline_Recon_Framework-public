

function xOut = pilotTonePrediction (A, xIn, typeConversion, calibrationModel, NCha, NMotParam)
% Predict PT as a NCha*NStates signal

if nargin<3 || isempty(typeConversion); typeConversion='T2PT';end%Motion parameters to PT signal (T2PT) or from PT signal to motion parameters (PT2T)
if nargin<4 || isempty(calibrationModel); calibrationModel='backward';end %Modelling motion as a linaer combination of PT signal T = Ab*p (backward) or PT as p=A*T (forward = physical setup)
if nargin<5 || isempty(NCha); NCha=[];end %Just a check to see if combination of A, typeConversion and calibrationModel is consistent
if nargin<6 || isempty(NMotParam); NMotParam=[];end %Just a check to see if combination of A, typeConversion and calibrationModel is consistent

calibToDegrees=1;%Test for rot-ran imbalance

%%% CHECK DATA CONSISTENCY
calibrationOffset=0;
if strcmp(calibrationModel,'forward') %p = A*T
    if ~isempty(NCha); assert(size(A,1)==NCha,'pilotTonePrediction:: A, typeConversion or calibrationModel not internally consistent');end
    if size(A,2)==(NMotParam+1); calibrationOffset=1;end
elseif strcmp(calibrationModel,'backward') %T = A*p
     if ~isempty(NCha); assert(size(A,2)==NCha || size(A,2)==(NCha+1),'pilotTonePrediction:: A, typeConversion or calibrationModel not internally consistent');end
     if size(A,2)==(NCha+1); calibrationOffset=1;end
end

%%% CONVERT
if strcmp(typeConversion,'T2PT') %From motion parameters to PT signal   
    T = xIn;
    N = size(T);
    nD = ndims(T);
    T = permute(T,[nD, nD-1, 1:nD-2]);% Permute to have as 6 * NStates
    if calibToDegrees;T = dynInd(T,4:6,1,convertRotation( dynInd(T,4:6,1),'rad','deg') );end%In degrees
    
    if strcmp(calibrationModel, 'forward') %p = A*T
         if calibrationOffset; T = cat(1, T, onesL(dynInd(T,1,1)));end%Offset in motion parameters
         p = A * T;
        
    elseif strcmp(calibrationModel, 'backward') %T = A*p
        if isequal( A, zerosL(A)); error('pilotTonePrediction:: Cannot estimate motion parameters with calbration matrix of zeros.');end
        p = pinv(A)*T;        
    else
        error('pilotTonePrediction:: %s not recognised as a calibrationModel.',calibrationModel)
    end
    p = dynInd(p, 1:NCha,1);%extract DC term
    xOut = p;
    
elseif strcmp(typeConversion, 'PT2T') %From PT signal to motion parameters 
    
    p = xIn;
    isComplex = ~isreal(p) || ~isreal(A);
        
    if strcmp(calibrationModel, 'forward') %p = A*T
        if isequal( A, zerosL(A)); warning('pilotTonePrediction:: Cannot estimate motion parameters with calbration matrix of zeros.');end
        T = pinv(A)*p;
        
    elseif strcmp(calibrationModel, 'backward') %T = A*p
        if calibrationOffset; p = cat(1, p, onesL(dynInd(p,1,1)) );end%Offset in PT signal (along channels)
        T=A*p;        
    else
        error('pilotTonePrediction:: %s not recognised as a calibrationModel.',calibrationModel)
    end
    
    if calibToDegrees;T = dynInd(T,4:6,1,convertRotation( dynInd(T,4:6,1),'deg','rad') );end%In degrees
    
    if isComplex; T = real(T);end
    T = ipermute(T, [6 5 1:4]);    
    T = dynInd(T, 1:NMotParam,6);%extract offset
    
    xOut=gather(T);
    
else
    error('pilotTonePrediction:: %s not recognised as a typeConversion.',typeConversion)
    
end