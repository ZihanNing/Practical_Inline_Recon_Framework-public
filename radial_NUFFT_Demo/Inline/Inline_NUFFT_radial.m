function [arrayimg1_sos] = Inline_NUFFT_radial(PIPE_folder,twix_like)

%===========================================================
% Main execution script for MERINA sodium MRI k-space to im-
% age reconstruction. Handles various settings and options, 
% reconstructs images with optional trajectory shift corr.,
% combines channel images, saves images and settings. If 
% applicable, aligns multiple repeats in image space. 
%-----------------------------------------------------------
% Sam Rot (UCL)
%-----------------------------------------------------------
% Notes: currently only has Finufft recon supported, but pl-
% anning to include legacy reconstruction for completeness. 
% Some extensions to come: B1 mapping and correction, B0 map-
% ping and correction, improved trajectory realignment. 
%-----------------------------------------------------------
% This is an inline implementation version PIPE
% implemented by Zihan Ning (KCL)
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

%%% recon of volume image %%% 
if skip_vol_recon == false
    if shift_img == "Vol" && skip_shift == false
        doshift = 'y'; 
    else
        doshift = 'n';
    end
    disp("Commenced recon of volume image")
    [img_vol, res_vol, ~, shift, shift_TE1] = methods_recon([raw_path,'/',vol_raw_fname],'',Tro_vol,"n","n",[],doshift,shift_noTE,[],[],shift_img,"n",meth_nufft, meth_script);
    disp("Completed recon of volume image")
end

%%% recon of array image(s) %%%
if contains(shift_img,"Arr") && skip_shift == false
    doshift = 'y'; 
else
    doshift = 'n';
end

if ~(exist('shift', 'var')) % if the shift wasn't calculated from volume images
    shift = [];
    shift_TE1 = [];
end

if no_array_imgs == 2
    if meth_avg == "Image"
        disp("Commenced recon of array images")
        [img_array, res_array, ~, ~, ~] = methods_recon([raw_path,'/',array_raw1_fname],'',Tro,"n","n",[],doshift,shift_noTE,shift,shift_TE1,shift_img,"n",meth_nufft, meth_script);
        disp("Reconstructed array image 1/2")
        img_array2 = methods_recon([raw_path,'/',array_raw2_fname],'',Tro,"n","n",[],doshift,shift_noTE,shift,shift_TE1,shift_img,"n",meth_nufft, meth_script);
        disp("Completed recon of array images")
    else
        disp("Commenced recon of array images")
        [img_array, res_array, ~, ~, ~] = methods_recon([raw_path,'/',array_raw1_fname],[raw_path,'/',array_raw2_fname],Tro,"n","n",[],doshift,shift_noTE,shift,shift_TE1,shift_img,"n",meth_nufft, meth_script);
    end
else
    disp("Commenced recon of array image")
    [img_array, res_array, ~, ~, ~] = methods_recon_inline(twix_like,Tro,"n","n",[],doshift,shift_noTE,shift,shift_TE1,shift_img,"n",meth_nufft, meth_script);
    disp("Completed recon of array image")
end

%%% sos combination
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