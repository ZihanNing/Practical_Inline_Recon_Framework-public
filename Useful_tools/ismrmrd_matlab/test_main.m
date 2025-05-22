clearvars; close all; clc;
addpath(genpath(pwd))

% define the filename
filename = 'input_highres_7T_MP2RAGE.h5';

% dataset definition
dset = ismrmrd.Dataset(filename, 'dataset');

% Header retrieving
header = ismrmrd.xml.deserialize(dset.readxml);

% Data retrieving
data_struct = dset.readAcquisition();
clearvars dset;

disp('Simple ismrmrd commands allows to retrieve the main Header and the Raw data section of the MRD dataset');
disp('Press Enter to continue !');
disp('.');
pause();