
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a script to perform MoCo on volumetric GRE-based sequences.
%Code provided by:
%Zihan Ning,
%Yannick Brackenier, and
%Lucilio Cordero-Grande.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET PATHS
addpath(genpath(pwd));

%% SPECIFY THE ONES YOU WANT TO RECONSTRUCT
idFile = 1; %List of all the acquisitions to reconstruct
idRef = 1; %Should not be touched - only change if you acquired different ref scans and want to change the one to be used
idRefB = 0; % added by ZN, if there's a refB need to be used

%% RUN
seq_type = 'mege';% mprage, mege supported for this version of the code
tic
runRecon(idFile, idRef,idRefB, seq_type)
Compute_Time = toc/60
