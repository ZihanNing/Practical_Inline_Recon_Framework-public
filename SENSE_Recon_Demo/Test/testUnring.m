addpath(genpath('/home/lcg13/Work/DISORDER7TR02'))

fileIn='/home/lcg13/Mounts/specialist/Data/Siemens/fMRI_2023_07_11/20230711_FID41023_MS_wrist_fMRI_dpg_on_08mmiso_GRAPPA2_20mm_AX_150meas_x.nii.gz';
fileOu='/home/lcg13/Mounts/specialist/Data/Siemens/fMRI_2023_07_11/20230711_FID41023_MS_wrist_fMRI_dpg_on_08mmiso_GRAPPA2_20mm_AX_150meas_x_Unring.nii.gz';

[x,~,MT]=readnii(fileIn);
x=unring(x);
writenii(fileOu,x,[],[],MT);
