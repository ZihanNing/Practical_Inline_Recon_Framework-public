clear all
close all
clc

%% This code is for recovering the data and prepare a structure for retro-recon
%%% prepare main_to_sub.mat
patientID = '33488490';
Acq_date = '08-Aug-2024';
ori_path = '/home/zn23/gadgetron_test/DATA_DUMP/';
nameRef = strcat(ori_path,Acq_date,'-P',patientID,'/RAW_DATA_HIGHRES.mat');
save_path = strcat(ori_path,Acq_date,'-P',patientID);
isSWI = 0;
fileName = 'T1_MPRAGE_DISORDER_sag_Gad';
save('main_to_sub.mat','patientID','nameRef','save_path','isSWI','fileName');

%% Call subcall bash for a recomputation
handle_connection_disorder_subcall_bash;
