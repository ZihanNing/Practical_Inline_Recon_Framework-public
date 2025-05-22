clc
clear
close all

%% read in raw data
home_path = getenv('HOME');
addpath(genpath([home_path,'/matlab/Reconstructionv4/mapVBVD']));
twix = mapVBVD('meas_MID00750_FID57554_SWI_DISORDER_R2_TS58_still.dat');
twix = twix{2};
twix.image.dataSize(1:11)

%% arrange the header
TW_hdr = twix.hdr;
TW_image.Lin = twix.image.Lin;
TW_image.Par = twix.image.Par;
TW_image.NCol = twix.image.NCol;
TW_image.centerCol = twix.image.centerCol;
TW_image.centerLin = twix.image.centerLin;
TW_image.centerPar = twix.image.centerPar;
TW_image.sqzSize = twix.image.sqzSize;
TW_image.dataSize = twix.image.dataSize;
TW_image.slicePos = twix.image.slicePos;
TW_image.Ave = twix.image.Ave;

%% save
save('twix_SWI_UKBB','TW_hdr','TW_image')

