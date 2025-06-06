
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a script to perform MoCo on volumetric GRE-based sequences.
%Code provided by:
%Yannick Brackenier and 
%Lucilio Cordero-Grande.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET PATHS
addpath(genpath(pwd));

%% SPECIFY THE ONES YOU WANT TO RECONSTRUCT
idFile = 1; %List of all the acquisitions to reconstruct
idRef = 1; %Should not be touched - only change if you acquired different ref scans and want to change the one to be used
idRefB = 0; % added by ZN, if there's a refB need to be used

%% RUN
seqType = 0;%seqType=0 for GRE/MPRAGE, 1 for multi-echo GRE and 2 for MP2RAGE
tic
runRecon(idFile, idRef,idRefB, seqType)
Compute_Time = toc/60
