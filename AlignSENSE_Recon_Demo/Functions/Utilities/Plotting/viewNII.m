

function [] = viewNII(fileName, rr)

%% load nii
nii = load_nii(fileName);
% nii.hdr.hist.srow_x=[0 1 0 0];
% nii.hdr.hist.srow_y=[0 1 0 0];
% nii.hdr.hist.srow_z=[0 0 1 0];
% nii.hdr.hist.sform_code=1;
% nii.original.hdr.hist.srow_x=[0 1 0 0];
% nii.original.hdr.hist.srow_y=[0 1 0 0];
% nii.original.hdr.hist.srow_z=[0 0 1 0];
if nargin < 2 ; rr =  [];end


%% Options and display
option = [];
% option.command = 'init'
% option.command = 'update'
% option.command = 'clearnii'
% option.command = 'updatenii'
% option.command = 'updateimg' (nii is nii.img here)
 option.usecolorbar = 1;%0 | [1]
% option.usepanel = 0 | [1]
% option.usecrosshair = 0 | [1]
% option.usestretch = 0 | [1]
% option.useimagesc = 0 | [1]
% option.useinterp = [0] | 1
% 
% option.setarea = [x y w h] | [0.05 0.05 0.9 0.9]
% option.setunit = ['vox'] | 'mm'
% option.setviewpoint = [x y z] | [origin]
% option.setscanid = [t] | [1]
% option.setcrosshaircolor = [r g b] | [1 0 0]
 option.setcolorindex = 3;%From 1 to 9 (default is 2 or 3)
% option.setcolormap = colormapNII(1000); %(Mx3 matrix, 0 <= val <= 1)
% option.setcolorlevel = No more than 256 (default 256)
% option.sethighcolor = []
% option.setcbarminmax = rr;
% option.setvalue = []
 option.glblocminmax = rr;%[]
% option.setbuttondown = ''
% option.setcomplex = [0] | 1 | 2

view_nii(nii, option);


function cmap= colormapNII(N)

N1 = round(N/2);
N2 = N-N1;

red= [1 0 0 ];
blue = [0 0 1];
white= [1 1 1];
black = [0 0 0];

cmap=linspaceND(black,white,N);
cmap = permute(cmap,[3 2 1]);
