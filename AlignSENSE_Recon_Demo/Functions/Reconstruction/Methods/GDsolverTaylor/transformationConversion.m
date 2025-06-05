

function [E, TrConverted] = transformationConversion(E, Tr, Dstruct, DstructName, block, di, N)

%TRANSFORMATIONCONVERSION Converts the transformation parameters from the logical coordinates (IJK) (units of voxels and radians around the centre), to
%   transformation parameters in the Right-Anterior-Superior (RAS) frame w.r.t. the scanner's iso-centre.
%   [E,TRCONVERTED]=TRANSFORMATIONCONVERSION(E,TR,{DSTRUCT},{DSTRUCTNAME},{BLOCK},{DI},{N})
%   * E is the encoding structure as used in the reconstructions (see solveXTB.m).
%   * {TR} the motion parameters with the last 2 dimensions containing the #states and the #motion parameters (6 for rigid transformation). Only rigid motion implemented.
%   * {DSTRUCT} a substructure of E used for B0 modelling. Can be created independently when not doing reconstructions involving the E structure.
%             Look at the code to see what needs to be provided in DSTRUCT.
%   * {DSTRUCTNAME} the name of the D structure so that it can be assigned to E.
%   * {BLOCK} a flag to block the conversion.
%   * {DI} a flag to indicate the direction. 
%           Forward conversion(1=default) is from IJK --> RAS.
%           Inverse conversion (o) is from RAS --> IJK.
%   * {N} the size of Dstruct.Dl. This allows you to run this function without having to store a large array for array size computation only. 
%   ** E the updated encoding structure.
%   ** TRCONVERTED the converted motion parameters.
%

if nargin < 1 || isempty(E); E = []; end
if nargin < 2 || isempty(Tr); assert(isfield(E,'Tr'),'transformationConversion:: Motion parameters must be provided explicitely or in E structure');Tr = E.Tr; end
if nargin < 3 || isempty(Dstruct); Dstruct = []; end
if nargin < 4 || isempty(DstructName); DstructName = 'Dl'; end %Taylor model. Function also used for Ds (susceptibility modelling).
if nargin < 6 || isempty(di); di = 1; end
if nargin < 7 || isempty(N); N = []; end

%%% CHECK INPUT DATA
if ~isempty(Dstruct) 
    E.(DstructName) = Dstruct;
elseif isempty(Dstruct) && ~isfield(E,DstructName)
    warning('transformationConversion:: No %s structure found in either E or as input argument.\n Return input motion parameters.',DstructName);
    TrConverted=Tr;
    return;
end

if nargin < 5 || isempty(block); if isfield(E.(DstructName).Geom, 'block'); block=E.(DstructName).Geom.block; else; block = 0; end;end
if isfield( E.(DstructName).Geom, 'useLog');useLog = E.(DstructName).Geom.useLog;else; useLog=0;end
if useLog; fprintf('transformationConversion:: Logarithm of rotation matrix used in conversion. Should only be use for PT studies.\n');end

%%% IF BLOCKED, FILL IN ORIGINAL TR AND STORE IN E
if block
    TrConverted = Tr;
    if isfield(E,DstructName); E.(DstructName).Tr = Tr; end
    return; 
end

%%% INITIALISE
TrConverted = zerosL(Tr);
if isempty(N)
    if strcmp(DstructName,'Dl'); N = size(E.(DstructName).D);
    elseif strcmp( DstructName,'Ds');N = size(E.(DstructName).chi);
    else; error('transformationConversion:: Not defined for other structure');
    end
end
N = N(1:3);
dimM = numDims(Tr);
dimS=dimM-1;
NStates = size(Tr,dimS);

%%% CREATE COMBINED TRANFORMATION MATRIX
Tcomb = E.(DstructName).Geom.APhiRecOrig;%s_form = Active meaning: from RAS -> space where motion is applied

if isfield(E.(DstructName).Geom,'permuteHist')
    for q = 1:length(E.(DstructName).Geom.permuteHist)%Permutations
        perm = E.(DstructName).Geom.permuteHist{q};
        P = Tpermute(perm);
        Tcomb = Tcomb * P;%Permute in moving frame, so right multiplication
    end
end

%%% SCALING DUE TO LOWER RESOLUTION - deprecated since T parameters always defined on full resolution.
% if isfield(E.(DstructName).Geom,'permuteHist')
%     S = diag([E.(DstructName).Geom.levelScale 1]);%Scaling due to resampling to (isotropic) volume in between levels but not changing motion parameters accordingly)
%     Tcomb = Tcomb*S;
% end

%%% COMPUTE ROTATION IN RAS FRAME
for i = 1:NStates
    %Compute R from T
    Ti = dynInd(Tr,i,dimS);
    R = convertT(Ti,'T2R'); %transformation matrix for sincRigidTransform.m
    
    %Take into account FOV and optionally use logarithmic representation
    centreFOV = ceil((N+1)/2);
    
    if di==1%From ijk-RAS
        R = Rt(centreFOV) * R * Rt(-centreFOV) ; %Correct for fact that DISORDER transformation convention around center FOV.    %R = R(1:3,1:3);%For now only focus on rotation component (so inverses on 3x3 matrix instead of 4x4)
        RNew = (Tcomb * R)/Tcomb; %same as Tcomb * R * inv(Tcomb);

    elseif di==0 %From RAS-ijk
        if useLog 
            tRot = logarithmExtraction(Ti,0);
            R(1:3,1:3)=tRot;
        end
        RNew = Tcomb\(R * Tcomb); %same as inv(Tcomb) * R * Tcomb;
        RNew = Rt(-centreFOV) * RNew * Rt(centreFOV) ; %Correct for fact that DISORDER transformation convention around center FOV
    end
    
    %Convert to Euler Angles variable again
    TrConverted = dynInd(TrConverted, i,dimS, convertT(RNew,'R2T'));
    
    if useLog && di==1
        logParams = logarithmExtraction(RNew);
        TrConverted = dynInd (TrConverted, {i, 4:6}, dimS:dimM, logParams);
    end
    
end

%%% STORE IN E STRUCTURE
if di==1 %Only store in D structure if it's the forward conversion. Inverse conversion currently only applied without E.
    E.(DstructName).Tr = TrConverted;
end
end


%%% HELPER FUNCTIONS
function xou = logarithmExtraction(xin, di)
if nargin <2 || isempty(di); di=1;end
    if di
        tRot=xin; %Input is rotation matrix -> convert to log parameters
        logtRot = logm(tRot(1:3,1:3));
        logParams = cat(2, logtRot(1,2), logtRot(1,3) , logtRot(2,3) );
        xou = logParams;
    else
        Ti = xin; %Input is motion parameter with logaritmic weights stored --> convert to rotation
        logParams = dynInd(Ti,4:6,6); logParams = reshape(logParams,[1 3]);
        logtRot = zeros([3 3]);
        logtRot(1,2)= logParams(1); logtRot(2,1)=-logtRot(1,2);
        logtRot(1,3)= logParams(2); logtRot(3,1)=-logtRot(1,3);
        logtRot(2,3)= logParams(3); logtRot(3,2)=-logtRot(2,3);
        tRot = expm(logtRot);
        xou = tRot;
    end
end