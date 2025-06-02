
clc
cd ('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Projects\GeomTesting')
addpath(genpath('.'));
addpath(genpath('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Software\DISORDER\DefinitiveImplementationRelease07'));
addpath(genpath('C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Software\Utilities'));

%% Load data for specified acquisition
scan = 4;
if scan ==1 
    fileName = 'DataSet1/202103021142_AWP79013_4_GRE_MIDRES_DISO_16_16_LIN_LIP/s004a1001.nii';%Default
    fileNameRaw = 'DataSet1\20210302_FID18142_P_GRE_MIDRES_DISO_16_16_LIN_LIP.dat';
    dimFlip = [];
elseif scan ==2
    fileName = 'DataSet1/202103021142_AWP79013_20_GRE_DISO_PARot10v2/s020a1001.nii';%PA10. Attention: The actual rotation is AP, not PA
    fileNameRaw = 'DataSet1\20210302_FID18151_P_GRE_DISO_PARot10v2.dat';
    dimFlip = [3];
elseif scan ==3
    fileName = 'DataSet1/202103021142_AWP79013_16_GRE_DISO_HFRot10v2/s016a1001.nii';%HF10
    fileNameRaw = 'DataSet1\20210302_FID18149_P_GRE_DISO_HFRot10v2.dat';
    dimFlip = [3];
elseif scan ==4
    fileName = 'DataSet1/202103021142_AWP79013_18_GRE_DISO_LRRot10v2/s018a1001.nii';%LR10 
    fileNameRaw = 'DataSet1\20210302_FID18150_P_GRE_DISO_LRRot10v2.dat';
    dimFlip = [2];
end

%% Load image
[x,MS,MT] = readNIIExt(fileName(1:end-4),{''});
x{1} = single(x{1});

nameSave = 'C:\Users\ybr19\OneDrive - King''s College London\2020 - 2023 PhD\Projects\GeomTesting\Testing\';
writeNII (strcat(nameSave,'X') ,{'orig'},x,MS,MT);

%% Resample
N = size(x{1});
resol = MS{1}*3;
Nnew = round(N(1:3).*MS{1}./resol);
Nnew = N;
Nnew(1) = round(N(1)/.2);
Nnew(2) = round(N(2)/.2);
Nnew(3) = round(N(3)/.2);


xNew{1} = resampling(x{1}, Nnew);
xNew{1} = imresize3(x{1}, Nnew,'cubic');
[MSOu{1}, MTOu{1}] = mapNIIGeom(MS{1}, MT{1}, 'resampling', [], N, Nnew );
writeNII (strcat(nameSave,'X') ,{'resampled'},xNew,MSOu,MTOu);

%% Flip
% flipDim = [2];
% 
% xNew{1} = x{1};
% for fl = flipDim 
%     xNew{1} = flip(xNew{1},fl);
% end
% 
% [MTOu{1}, MSOu{1}] = mapNIIGeom(MT{1}, MS{1}, 'flip', flipDim, size(x{1}) , size(xNew{1}) );
% writeNII (strcat(nameSave,'X') ,{'flipped'},xNew,MSOu,MTOu);

%% Permute
% perm = [2 1 3];
% 
% xNew{1} = permute( x{1}, perm);
% 
% [MSOu{1}, MTOu{1}] = mapNIIGeom(MS{1}, MT{1}, 'permute', perm, size(x{1}) , size(xNew{1}) );
% writeNII (strcat(nameSave,'X') ,{'permute'},xNew,MSOu,MTOu);

%% dynInd
% idx = {23:60, ':' , ':' };
% 
% xNew{1} = dynInd( x{1}, idx, 1:3);
% [MTOu{1}, MSOu{1}] = mapNIIGeom(MT{1}, MS{1}, 'dynInd', idx, size(x{1}) , size(xNew{1}) );
% writeNII (strcat(nameSave,'X') ,{'dynInd'},xNew,MSOu,MTOu);

%% padding
% padDim = [20 30 60];
% 
% xNew{1} = padarray( x{1}, padDim, 0);
% [MTOu{1}, MSOu{1}] = mapNIIGeom(MT{1}, MS{1}, 'pad', padDim, size(x{1}) , size(xNew{1}) );
% writeNII (strcat(nameSave,'X') ,{'padding'},xNew,MSOu,MTOu);
% 
%% translate
% tran = [-10 4 7];
% 
% xNew{1} = circshift( x{1}, tran);
% [MSOu{1}, MTOu{1}] = mapNIIGeom(MS{1}, MT{1}, 'translate', tran, size(x{1}) , size(xNew{1}) );
% writeNII (strcat(nameSave,'X') ,{'translate'},xNew,MSOu,MTOu);


%% rotate

theta = [15 -18 30]/180*pi;
T = [ 0, 0, 0, theta(3), theta(2), theta(1) ];% [ t_x t_y t_z theta_z theta_y theta_x] in voxels and radians
R = T2R(T,1);

[~,kGrid,rkGrid,~,~] = generateTransformGrids(size(x{1}));  
[et] = precomputeFactorsSincRigidTransform(kGrid,rkGrid,T, 1); % Forward

xNew{1} = real(sincRigidTransform(x{1}, et, 1, [],[],0));
[MSOu{1}, MTOu{1}] = mapNIIGeom(MS{1}, MT{1}, 'rotate', R, size(x{1}) , size(xNew{1}) );
writeNII (strcat(nameSave,'X') ,{'rotate'},xNew,MSOu,MTOu);

xTest{1} = mapVolume(xNew{1},x{1},MTOu{1},MT{1});
writeNII (strcat(nameSave,'X') ,{'rotateBackMapped'},xTest,MS,MT);
