function [arrayimg1_sos] = Inline_NUFFT_radial(PIPE_folder,twix_like)

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
% already handled in the handle_connection_background_nufft_radial.m

%% UI FOR INPUTS
% the inline implemented pipeline should be completely automatic 
% therefore the UI has been replaced by a config

% load the config
load([PIPE_folder,'/Inline/Config_NUFFT_radial_MEGE.mat']);

%% PARSING INPUT %%
out_dir = twix_like.hdr.save_path; % save the intermediate result in the same folder as raw
meth_script = settings.meth_script;
meth_nufft = settings.meth_nufft;
finufft_dir = settings.finufft_dir;
no_array_imgs = settings.no_array_images;
array_raw1_fname = settings.rawarrname1;
if no_array_imgs == 2
    array_raw2_fname = settings.rawarrname2;
end
meth_avg = settings.meth_avg;
Tro = settings.Tro;
no_TE = settings.no_TE; 
meth_sense = settings.meth_sense;
sense_noPointscalib = settings.sense_noPointscalib; 
path_fsl = settings.path_fsl; 

skip_vol_recon = settings.check_volume;
vol_raw_fname = settings.rawvolname;
Tro_vol = settings.Tro;
no_TE_vol = settings.no_TE; 
sense_GaussKernel = settings.sense_GaussKernel;

skip_adc = settings.check_adc;

skip_b1_recon = settings.check_b1;
b1_raw_fname = settings.rawb1name;
Tro_b1 = settings.Tro;
no_TE_b1 = settings.no_TE; 

skip_shift = settings.check_shift;
shift_noTE = settings.shift_noTE;
shift_img = settings.shift_img; 

%% DIRECTORY MANAGEMENT %% 
% has been already handled


%% KSPACE TO IMAGE %%
% recon of array image(s)
doshift = 'n'; % trajectory shifting not supported in this demo

disp("Commenced recon of array image")
[img_array, res_array, ~, ~, ~] = methods_recon_inline(twix_like,Tro,"n","n",[],doshift,shift_noTE,shift,shift_TE1,shift_img,"n",meth_nufft, meth_script);
disp("Completed recon of array image")

%% SOS CHANNEL COMBINATION %% 
disp("Commenced sos channel combination")
img_array_sos = squeeze(sqrt(sum(abs(img_array).^2,4)));
disp("Completed sos channel combination")

%% SAVING IMAGES %% 
res_array = res_array * 1e3; 
output_path = twix_like.hdr.save_path; % uniformly saved the images all together
disp(['Writing final images as NIFTI file into output path ',output_path]);

arrayimg1_sos = fliplr(abs(squeeze(img_array_sos)));
tmp_nii = make_nii(arrayimg1_sos,[res_array,res_array,res_array]);%,[],[],mat2str(TE*1e3));
save_nii(tmp_nii,[output_path,'/arrayimg1_sos','.nii.gz']);

end