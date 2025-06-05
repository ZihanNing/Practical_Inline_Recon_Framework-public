
function [E] = updateShim(E, shimPar, NX, MT)

%%%PARAMETERS
addpath(genpath('/home/ybr19/Projects/B0Shimming'))
parameters_ShimCalibration;
debug=1;

%%%CREATE NEW BASIS
[~, E.Shim.coefNamePaper,E.Shim.coefNameTerra] = shimBasis(NX, MT);%Basis not stored to save memmory
E.Shim.NX=NX;
E.Shim.MT=MT;

%%% LOAD PRECALIBRATED DATA
correctionFactorsFileName = strcat(saveDir,filesep, '/correctionFactors.mat');
if exist(correctionFactorsFileName,'file')==2
    ss = load(correctionFactorsFileName); 
    correctFact_mA2au = ss.correctFact_mA2au;%Multiplicative correction
else
    warning('updateShim:: No shim parameters found from calibration.')
end

crossTermFileName = strcat(saveDir,filesep, 'crossTerms.mat');
if exist(crossTermFileName,'file')==2
    ss = load(crossTermFileName); 
    crossTerms = ss.crossTerms;
else
        warning('updateShim:: No shim parameters found from calibration.')
end
ss = load(strcat( saveDir, '/tuneUpSettings.mat')  );

%%%CREATE SHIM B0
shimPar_TuneUp = [];
shimPar_TuneUp.shimCurrents_au  = ss.cTuneUp;%Cotains >13 coefficients since it models the Tune-Up field
shimPar_TuneUp.shimOrder = shimOrder;
shimPar_TuneUp.gamma = ss.gamma;

E.Shim.B0 = getShim(NX, shimPar_TuneUp, MT, 'tune-up', debug);%B0 from Tune-Up
%E.Shim.B0_TuneUp = getShim(NX, shimPar_TuneUp, MT, 'tune-up', debug);

shimPar.shimCurrents_au  = shimPar.shimCurrents_au - ss.shimCurrents_au;%w.r.t. tuneup. Contains 13 coefficients. using crossTerms this will create >13 coeficients
shimPar.correctFact_mA2au = correctFact_mA2au;
shimPar.shimOrder = shimOrder;
shimPar.crossTerms=crossTerms;
E.Shim.B0 = E.Shim.B0 + getShim(NX, shimPar, MT, 'active', debug);%Tune-up includes passive shim + tune-up currents
%E.Shim.B0_Rel = E.Shim.B0 + getShim(NX, shimPar, MT, 'active', debug);%Tune-up includes passive shim + tune-up currents        
