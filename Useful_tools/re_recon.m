function re_recon()
%%
clc
clear all

sorting_info.fileName = 'T1_MPRAGE_DISORDER_sag_Gad_17-25-11';
sorting_info.nameRef = '/home/zn23/data/20-Mar-2025-Pxxxxxxxxxxxxxxxxxxxxxxxx/RAW_DATA_HIGHRES_T1_MPRAGE_DISORDER_sag_Gad_17-25-11.mat';
sorting_info.patient_ID = 'xxxxxxxxxxxxxxxxxxxxxxxx';
sorting_info.save_path = '/home/zn23/data/20-Mar-2025-Pxxxxxxxxxxxxxxxxxxxxxxxx';
sorting_info.save_time = '2025-03-20-10-26-34';

%% the SWI DIS (bucket)
sorting_info.fileName = 'SWI_DISORDER_R2_TS58_still_13-04-45';
sorting_info.nameRef = '/home/gadgetron/data/28-Feb-2025-Pxxxxxxxxxxxxxx/RAW_DATA_HIGHRES_SWI_DISORDER_R2_TS58_still_13-04-45.mat';
sorting_info.patient_ID = 'xxxxxxxxxxxxxx';
sorting_info.save_path = '/home/zn23/data/28-Feb-2025-Pxxxxxxxxxxxxxx';
sorting_info.save_time = '2025-02-28-10-26-34';

%% the SWI DIS (buffer)
sorting_info.fileName = 'SWI_DISORDER_R2_TS58_still_17-15-59';
sorting_info.nameRef = '/home/zn23/data/24-Feb-2025-Pxxxxxxxxxxxxxx/RAW_DATA_HIGHRES_SWI_DISORDER_R2_TS58_still_17-15-59.mat';
sorting_info.patient_ID = 'xxxxxxxxxxxxxx';
sorting_info.save_path = '/home/zn23/data/24-Feb-2025-Pxxxxxxxxxxxxxx';
sorting_info.save_time = '2025-02-24-10-26-34';

%% the SPACE (buffer)
sorting_info.fileName = 'FLAIR_sag_DISORDER_2Ave_PFoff_RR22caipi_still_14-44-39';
sorting_info.nameRef = '/home/zn23/data/11-Mar-2025-Pxxxxxxxx/RAW_DATA_HIGHRES_FLAIR_sag_DISORDER_2Ave_PFoff_RR22caipi_still_14-44-39.mat';
sorting_info.patient_ID = 'xxxxxxxx';
sorting_info.save_path = '/home/zn23/data/11-Mar-2025-Pxxxxxxxx';
sorting_info.save_time = '2025-03-11-10-26-34';

%% General pipeline (SWI)
sorting_info.fileName = 'SWI_DISORDER_R2_TS58_still_16-21-40';
sorting_info.nameRef = '/home/zn23/data/19-Mar-2025-Pxxxxxxxxxxxxxx/RAW_DATA_HIGHRES_SWI_DISORDER_R2_TS58_still_16-21-40.mat';
sorting_info.patient_ID = 'xxxxxxxxxxxxxx';
sorting_info.save_path = '/home/zn23/data/19-Mar-2025-Pxxxxxxxxxxxxxx';
sorting_info.save_time = '2025-03-19-16-38-21';

%% MPRAGE
save('main_to_sub_MPRAGE.mat')
handle_connection_disorder_MPRAGE_subcall_bash

%% SWI (DISORDER - buffer)
save('main_to_sub_SWI.mat')
handle_connection_disorder_MEGE_subcall_bash

%% SWI (DISORDER - bucket)
save('main_to_sub_SWI.mat')
handle_connection_disorder_bucketybuffer_MEGE_subcall_bash

%% SWI (SENSE)
save('main_to_sub_SWI.mat')
handle_connection_disorder_MEGE_SENSE_subcall_bash

%% FLAIR
save('main_to_sub_FLAIR.mat')
handle_connection_disorder_SPACE_subcall_bash

%% General
save('main_to_sub.mat')
handle_connection_disorder_general_subcall_bash

