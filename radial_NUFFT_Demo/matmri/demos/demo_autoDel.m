%% demo for automatic determination of trajectory delay
%
%  Cite:
%
%  Dubovan PI, Baron CA. Model-based determination of the synchronization
%  delay between MRI and trajectory data. Magn Reson Med. 2022 Sep 26. doi:
%  10.1002/mrm.29460.  PMID: 36161333.
%
%  (c) Corey Baron, 2021

% To make things simple, we'll start with actual acquired data for truth image. 
load data_images.mat

% Pick a slice
nsl = 2;
X = X(:,:,nsl); 
Y = Y(:,:,nsl); 
Z = Z(:,:,nsl); 
im0 = im0(:,:,nsl);
b0map = b0map(:,:,nsl);
Crcvr = Crcvr(:,:,nsl,:);

% The receivers are actually virtual coils from coil compression, so we can
% only use a subset to save time for the purposes of this demo.
Crcvr = Crcvr(:,:,:,1:16);

% Load a trajectory. Can choose between spiral or epi. 
% The trajectory data was acquired at a dwell time of 1 us.
% The datatime variable represents the sample times for a typical MRI
% acquisition (dwell time = 2.5 us). 
% dataStartTime is the time that data acquisition would start for each
% trajectory (i.e., the field probe system began acquiring early)
docase = 0;
switch docase
    case 0
        load data_traj_epi_R3.mat
        dataStartTime = 365;
    case 1
        load data_traj_spiral_R4.mat
        dataStartTime = 50;
end
datatime = datatime + dataStartTime;

% Set spatial grid
grid.x = X; 
grid.y = Y;
grid.z = Z;

% Build sampling object for delayed trajectory
delApplied = 2; % microseconds
[phs_spha_del,phs_conc_del] = interpTrajTime(phs_spha,phs_conc,tdwelltraj,delApplied,datatime);
S = sampHighOrder(b0map,datatime(:),phs_spha_del',phs_conc_del',grid);

% Build receiver operator
R = rcvrOp(Crcvr,0);

% Sample image
data = S*(R*im0);

% Add some noise
ns_std = 1e-4;
data = data + ns_std*(randn(size(data)) + 1i*randn(size(data)));

%% Normal highorder recon using conjugate gradient method

% Build sampling object based on trajectory with no delay
fprintf('Performing recon with delay error...\n')
[phs_spha_del0,phs_conc_del0] = interpTrajTime(phs_spha,phs_conc,tdwelltraj,0,datatime);
S = sampHighOrder(b0map,datatime(:),phs_spha_del0',phs_conc_del0',grid);

opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
maxIt = 20;
opt.stopOnResInc = 0; 
[im, resvec, mse, xnorm, xdiff] = cgne(opFunc,data,[],maxIt); 

% Find delay automatically
fprintf('Finding delay...\n')
del0 = 0; % starting guess 
maxNit_cgne = length(resvec)-1;
delJumpFact = 3;
numCoarseSearch = 0;
[delFound, delSk_perIt] = findDelAuto(del0,data,tdwelltraj,phs_spha,phs_conc,datatime,grid,b0map,Crcvr,maxNit_cgne,delJumpFact,numCoarseSearch);
fprintf('True delay = %.2f us. Found delay = %.2f us\n', delApplied, delFound);

% Determine image with optimal delay
fprintf('Performing recon with found delay...\n')
[phs_spha_dela,phs_conc_dela] = interpTrajTime(phs_spha,phs_conc,tdwelltraj,delFound,datatime);
S = sampHighOrder(b0map,datatime(:),phs_spha_dela',phs_conc_dela',grid);
opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
imFoundDel = cgne(opFunc,data,[],maxIt); 

% Determine image with true delay
fprintf('Performing recon with ground truth delay...\n')
[phs_spha_del,phs_conc_del] = interpTrajTime(phs_spha,phs_conc,tdwelltraj,delApplied,datatime);
S = sampHighOrder(b0map,datatime(:),phs_spha_del',phs_conc_del',grid);
opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
imAppliedDel = cgne(opFunc,data,[],maxIt); 


%% Plotting
figure;
subplot(2,2,1); imagesc(permute(abs(im0),[2,1]),[0 0.006]); title('Ground truth'); axis('image'); axis('off'); colormap('gray');
subplot(2,2,2); imagesc(permute(abs(im),[2,1]),[0 0.006]); title('Delay error'); axis('image'); axis('off'); colormap('gray');
subplot(2,2,3); imagesc(permute(abs(imAppliedDel),[2,1]),[0 0.006]); title('True delay'); axis('image'); axis('off'); colormap('gray');
subplot(2,2,4); imagesc(permute(abs(imFoundDel),[2,1]),[0 0.006]); title('Found delay'); axis('image'); axis('off'); colormap('gray');









