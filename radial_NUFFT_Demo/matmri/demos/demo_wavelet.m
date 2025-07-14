%% sampHighOrder with wavelet regularization demo
%
%  (c) Corey Baron, 2020

% To make things simple, we'll start with actual acquired data for 
% the B0 map, receiver sensitivity, and a ground truth image
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
docase = 1;
switch docase
    case 0
        load data_traj_epi_R3.mat
        dataStartTime = 365;
    case 1
        load data_traj_spiral_R4.mat
        dataStartTime = 50;
end

% Set spatial grid
grid.x = X; 
grid.y = Y;
grid.z = Z;

% Interpolate trajectory samples to typical MRI sampling rate
[phs_spha,phs_conc] = interpTrajTime(phs_spha,phs_conc,tdwelltraj,dataStartTime,datatime);
S = sampHighOrder(b0map,datatime(:),phs_spha',phs_conc',grid);

% Build receiver operator
R = rcvrOp(Crcvr,0);

% Sample image
data = S*(R*im0);

% Add some noise
ns_std = 1e-4;
data = data + ns_std*(randn(size(data)) + 1i*randn(size(data)));


%% Normal highorder recon using conjugate gradient method
% Note that by default, cgne stops on the first iteration 
opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
maxIt = 25;
opt.stopOnResInc = 0; 
[im, resvec, mse, xnorm, xdiff] = cgne(opFunc,data,[],maxIt); 


%% Wavelet regularization
% Solve argmin(||Ax-b||^2_2 + lambda*||Wx||_1)
% Build discrete wavelet transform object
levelsPerDim = [2 2];
isDec = 0; % Non-decimated to avoid blocky edges
W = dwt(levelsPerDim,size(im),isDec);
wim = W*im;
figure; imagesc(wim); title('Wavelet transform of cgne solution')

% Automatic lambda selection, similar to Varela-Mattatall G, Baron CA,
% Menon RS. Magn. Reson. Med. 2021;86:1403–1419. Here we'll base the lambda
% selection per level on histograms of wavelet xform of "zero filled"
% image, which is more generally just the adjoint of the encoding matrix A.
% We just use the location where the histogram of the magnitude image
% is maximized, which corresponds to the std of the noise if we assume that
% the distribution is dominated by the {noise + aliasing} that we want
% to suppresss with wavelet regularization (i.e., it's roughly a Rayleigh
% distribution). We divide by sqrt(2) b/c we're dealing with magnitude
% while the bfista algorithm is using complex values for softthresholding.
% This assumption is most valid for the 1st level high freq transform, so
% we use this lambda for all the other levels (and set the lambda for the
% low frequencies to 0, since it is not sparse).
im_zf = opFunc(data,'transp');
wim_zf = W*im_zf;
figure; 
subplot(1,2,1);
imagesc(wim_zf); title('Wavelet transform of zero-filled');
subplot(1,2,2)
h = histogram(gather(abs(wim_zf.high{1}(:))),round(numel(wim_zf.high{1}(:))/100));
title('Histogram of magnitude of 1st level of wavelet transform')
xlim([0 6e-4]);
vals = h.Values;
vals(1) = 0; % The first bin can be inflated if there was masking in im0
[~,threshInd] = max(vals);
lamVal = 0.5*(h.BinEdges(threshInd) + h.BinEdges(threshInd+1))/sqrt(2);
lambda = waveletObj(lamVal,wim_zf); % This replicates the scalar lamVal into all wavelet fields using the template wim_zf
lambda.low = 0;

% Do a few cgne iterations to get a better starting guess
x0 = cgne(opFunc,data,[],2); 

% Run the recon. 
NitMax = 200;
opt.resThresh = 1e-4; % use fractional change in (||Ax-b||^2_2 + lambda*||Wx||_1) from last iteration to stop iterations
opt_bfista.gtruth = im0; % For these sims, we know the ground truth and can thus track mean square error
[imCS, resSqAll, RxAll, mseAll] = bfista(opFunc,data,W,lambda,x0,NitMax,opt_bfista);

%% Plotting. Note that imCS has much less noise than im.
figure;
subplot(2,3,1); imagesc(permute(abs(im0),[2,1]),[0 0.006]); title('Ground truth'); axis('image'); axis('off'); colormap('gray');
subplot(2,3,2); imagesc(permute(abs(im),[2,1]),[0 0.006]); title('Early stopping regularization'); axis('image'); axis('off'); colormap('gray');
subplot(2,3,3); imagesc(permute(abs(imCS),[2,1]),[0 0.006]); title('Wavelet regularization'); axis('image'); axis('off'); colormap('gray');
subplot(2,3,4); plot(log10(resSqAll(:,1))); ylabel('Data consistency'); xlabel('iterations');
subplot(2,3,5); plot(log10(RxAll(:,1))); ylabel('Wavelet l1 norm'); xlabel('iterations');
subplot(2,3,6); plot(log10(mseAll)); ylabel('Mean squared error'); xlabel('iterations');








