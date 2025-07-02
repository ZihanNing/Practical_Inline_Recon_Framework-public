%% Expanded encoding model using inversion of dense sampling matrix
%
%  Modern GPU's have enough memory to store the full dense encoding 
%  matrix for a single 2D slice. This example requires ~10GB GPU memory
%
%  (c) Corey Baron, 2022

% To make things simple, we'll start with actual acquired data for 
% the B0 map, receiver sensitivity, and a ground truth image
load data_images.mat

% Pick a slice
nsl = 1;
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
nse =  ns_std*(randn(size(data)) + 1i*randn(size(data)));
data = data + nse;


%% Normal highorder recon using conjugate gradient method
% Stop several iterations after the residual decreases below the noise norm
% (i.e., discrepancy principle). Stopping right at the point where the
% residual is equal to noise tends to over-regularize, so we do several 
% extra iterations
tic1 = tic;
opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
Atnse = opFunc(nse, 'transp');
maxIt = 50;
opt.noiseVar = 2*ns_std^2; % Factor of 2 for complex data 
opt.nseExtraIt = 6; % Extra iterations after reaching residual expected from noise 
[im1, resvec, mse, xnorm, xdiff] = cgne(opFunc,data,[],maxIt,opt); 
fprintf('Iterative method took %d seconds\n', round(toc(tic1)));


%% Direct method
% Clear S and R to save memory
clear S R

tic2 = tic;
% Create dense matrices
[AtA,Aty] = mrSampFuncMat(data,Crcvr,b0map,datatime(:),phs_spha',phs_conc',grid);
fprintf('Took %d seconds to find dense matrix\n', round(toc(tic2)));

% Tikhonov regularization  
% We find the regularization weighting that satisfies the discrepancy
% principle, where the residual is close to what is expected from the known
% noise variance. This is somewhat slow, but it only needs to be done for
% one sample slice.
lam2Vals = 10.^(-4:1:0);
resVals = zeros(length(lam2Vals),1);
regVals = zeros(length(lam2Vals),1);
im2All = zeros([size(im0), length(lam2Vals)]);
di = size(AtA,1);
for nlam = 1:length(lam2Vals)
    lam2 = lam2Vals(nlam);
    AtA(1:di+1:numel(AtA)) = AtA(1:di+1:numel(AtA)) + lam2; % The square of the lambda value just gets added to the diagonal.
    %clear S R 
    im2 = AtA\Aty(:);
    AtA(1:di+1:numel(AtA)) = AtA(1:di+1:numel(AtA)) - lam2;
    tmp = AtA*im2 - Aty(:);
    resVals(nlam) = norm(gather(tmp));
    regVals(nlam) = norm(gather(im2(:)));
    im2All(:,:,nlam) = double(gather(reshape(im2,size(im0))));    
end
figure; 
subplot(2,2,1); plot(log10(resVals.^2), log10(regVals.^2), '-o')
xlabel('log10(||residual||^2)'); ylabel('log10(||x||^2)'); title('L-curve')
subplot(2,2,3); plot(log10(lam2Vals), log10(resVals.^2), '-o')
xlabel('log10(lambda^2)');  ylabel('log10(||residual||^2)');
subplot(2,2,4); plot(log10(lam2Vals), log10(regVals.^2), '-o')
xlabel('log10(lambda^2)');  ylabel('log10(||x||^2)');

% Find regularization param with linear interpolation based on discrepancy
% principle (L-curve is also an option, but in my experience the corner is
% hard to find algorithmically, and it under-regularizes for these recons).
refVal = norm(Atnse(:));
ind_l = find((resVals - refVal) < 0, 1, 'last');
% y=mx+b on log of values
m = (log10(resVals(ind_l+1)) - log10(resVals(ind_l)))/...
    (log10(lam2Vals(ind_l+1)) - log10(lam2Vals(ind_l)));
b = log10(resVals(ind_l)) - m*log10(lam2Vals(ind_l));
lam2_log = (log10(refVal) - b)/m;
lam2 = (10^lam2_log);

% Discrepancy principle tends to over-regularize, so scale lambda by a
% factor. This should only weakly depend on traj and noise level.
lam2 = 0.5*lam2;

fprintf('Took %d seconds to find regularization param\n', round(toc(tic2)));
fprintf('  log10(residual) from discrepancy principle is %g, which yields lam2 ~ %g \n',...
    log10(norm(Atnse(:)).^2), lam2)

tic3 = tic;
AtA(1:di+1:numel(AtA)) = AtA(1:di+1:numel(AtA)) + lam2;
im2 = AtA\Aty(:);
im2 = gather(reshape(im2,size(im0)));
fprintf('Direct method took %d seconds\n', round(toc(tic3)));
clear AtA Aty

figure; 
subplot(2,2,1); imagesc(permute(im0, [2 1])); colormap('gray'); axis('image'); axis('off'); title('Ground truth')
subplot(2,2,3); imagesc(permute(abs(im1), [2 1])); colormap('gray'); axis('image'); axis('off'); title('Iterative')
subplot(2,2,4); imagesc(permute(abs(im2), [2 1])); colormap('gray'); axis('image'); axis('off'); title('Direct')


