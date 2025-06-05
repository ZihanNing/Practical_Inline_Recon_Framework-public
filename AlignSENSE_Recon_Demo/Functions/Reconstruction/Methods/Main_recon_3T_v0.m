
clc; close all;
clear all
cd /home/ybr19/Reconstruction/
addpath(genpath('.'))

studies_20210329_DISORDER_MotionFree_3T_Reconstructions;

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
rmpath(genpath('/home/ybr19/Data'))
addpath(genpath(pathData{1}))

%fileName = '/home/ybr19/Data/3Ttest/ha_17082019_1007414_3_2_mv4mpragesenseV4001.mat';
%fileName = '/home/ybr19/Data/3Ttest/gr_16082019_1359078_7_2_mv43dbrainviewflairshcsenseV4001.mat';
fileName = strcat(pathData{1}, filesep, fileIn{1}{1}, '001.mat');%Motion corrupted data
%fileName = strcat(pathData{1}, filesep, fileIn{3}{1}, '001.mat');%Hybrid data
load(fileName);

coilTests = 1:8;

for i=1:length(coilTests)  

%%% ASSIGN SENSITIVITIES
fileNameRef = strcat(pathData{1}, filesep, fileIn{1}{coilTests(i)}, '001.mat');%centre-centre
ss =load(fileNameRef);
rec.S = ss.rec.S;
ss = [];

%%% RUN
rec.Par.Labels.H0 = 3;
rec.Alg.parXB.useSH = 0;
rec.Alg.parXB.BasisType = 'SH';% 'SHOW' / 'BSplines' - not yet implemented ar
rec.Alg.parXB.SHorder = 3; % disable by setting < 0
rec.Alg.parXB.nGD_B = 1;
rec.Alg.parXB.deTaylor = [2 inf];
rec.Alg.parXB.Optim.basisDelay=4;  
     
rec.Alg.parXT.meanT=0;       
rec.Alg.parXT.groupSweeps = 2;
rec.Alg.parXT.resolMax = 5;
rec.Alg.parXT.saveFinal=0;
            rec.Names.Name=sprintf('scan1_coills_%i',coilTests(i));
            rec.Names.Name=sprintf('testOrient');
            rec.Names.pathOu=pathData{1};

rec.Alg.parXB.useTaylor = 0;
rec.Alg.parXB.orderTaylor = 1;
rec.Plan.Suff='_DephCorr';
rec.Plan.SuffOu= sprintf('_SH%dTaylor%d', rec.Alg.parXB.SHorder,rec.Alg.parXB.orderTaylor);
            
rec.Alg.parXB.Reg.lambda_sparse = 0e+5;%1e+5; % Ad hoc - need to refine
rec.Alg.parXB.Reg.lambda_smooth = 0e+5; % 3e+5 Ad hoc - need to refine
rec.Alg.parXB.Reg.lambda_susc = 0e+3;
rec.Alg.parXB.nGD_Taylor = 1;
rec.Alg.parXB.filterD.sp = 2;%In mm
rec.Alg.parXB.filterD.gibbsRinging = 0.4;
rec.Alg.parXB.deSH = 0;
rec.Alg.parXB.Optim.FISTA=0;
rec.Alg.parXB.Optim.RMSProp=0;
rec.Alg.parXB.Optim.beta=0.9;
rec.Alg.parXB.Optim.DupLim=inf;
rec.Alg.parXB.Optim.taylorDelay=15;
rec.Alg.nExtern = 2;

rec.Alg.parXB.weightedEst = 0;
rec.Alg.parXB.redFOV = 0;
rec.Alg.parXB.intraVoxDeph = 0;
rec.Alg.parXT.percRobustShot=[0.25 0.75];%Percentiles for robust computation of expected inter-shot dispersion
rec.Alg.parXT.enerRobustShot=0.85;%0.95;%Ratio of the error for acceptance of shots
rec.Alg.parXT.corrFact = 1;
rec.Alg.parXT.disableGrouping = 1;
rec.Alg.parXB.useSusc = 0;
            rec.Alg.parXB.Optim.alphaList=1e+10 * 10.^(-1* [-3:0.7:6]);
rec.Alg.parXB.useB1 = 0;

rec.Alg.parXT.convTransformJoint = 1;%all motion parameters estimated at every iteration
rec.parXT.writeInter=0;
            rec.Alg.parXB.redFOV = 1/(2.7);
            rec.Alg.parXB.intraVoxDeph = 0;

solveXTB_ext(rec);

end