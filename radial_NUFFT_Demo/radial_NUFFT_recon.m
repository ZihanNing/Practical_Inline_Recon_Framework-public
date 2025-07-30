%===========================================================
% Main execution script for MERINA sodium MRI k-space to im-
% age reconstruction. Handles various settings and options, 
% reconstructs images with optional trajectory shift corr.,
% (not supported in demo!), combines channel images (sos 
% only in demo), saves nifti images.
%-----------------------------------------------------------
% Sam Rot (UCL)
%-----------------------------------------------------------
% Notes: reduced functionality for demo
%===========================================================

%% DEPENDENCIES 
addpath(genpath('radial_NUFFT_Demo'));
addpath(genpath('/home/zn23/matlab/finufft')); % ZN: load the finufft path; finufft toolbox need to be installed first

%% CONFIGS FOR INPUTS
home_path = getenv('HOME')
raw_path = [home_path,'/Gadgetron_Parallel_Framework/raw'];

settings.out_dir = 'out2';
settings.rawarrname1 = 'meas_MID00165_FID10766_KCL_miniflash_JC_merina1000proj_tro2_spiral_37e_tept6_t245v_fa9.dat';
settings.Tro = 2;
settings.check_shift = 1;
settings.shift_img = 'Vol';
settings.shift_noTE = 12;

% PARSING INPUT %
out_dir = settings.out_dir;
array_raw1_fname = settings.rawarrname1;
Tro = settings.Tro;
skip_shift = settings.check_shift;
shift_noTE = settings.shift_noTE;
shift_img = settings.shift_img; 

%% DIRECTORY MANAGEMENT %% 
if not(isfolder([raw_path,'/',out_dir]))
    mkdir([raw_path,'/',out_dir])
    disp("Created image output directory")
else
    error("Output directory <" + out_dir + "> already exists. Either remove it, or enter a different directory name")
end

%% KSPACE TO IMAGE %%
% recon of array image(s)
doshift = 'n'; % trajectory shifting not supported in this demo

disp("Commenced recon of array image")
[img_array, res_array, ~, ~, ~] = methods_recon([raw_path,'/',array_raw1_fname],'',Tro,"n",doshift,shift_noTE,[],[],shift_img);
disp("Completed recon of array image")

%% SOS CHANNEL COMBINATION %% 
disp("Commenced sos channel combination")
img_array_sos = squeeze(sqrt(sum(abs(img_array).^2,4)));
disp("Completed sos channel combination")

%% SAVING IMAGES %% 
res_array = res_array * 1e3; 
output_path = [raw_path, '/', out_dir];
disp(['Writing final images as NIFTI file into output path ',output_path]);

tmp_nii = make_nii(fliplr(abs(squeeze(img_array_sos))),[res_array,res_array,res_array]);%,[],[],mat2str(TE*1e3));
save_nii(tmp_nii,[output_path,'/arrayimg1_sos','.nii.gz']);
