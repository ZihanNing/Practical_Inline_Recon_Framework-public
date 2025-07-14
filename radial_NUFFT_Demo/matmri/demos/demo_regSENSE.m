%% Regularized SENSE Demo
%
%  Try the regridding demo before this one. Here, we include receiver
%  sensitivities in the recon to implicitely perform SENSE. Cartesian and
%  non-Cartesian cases are handled similarly.
%
%  (c) Corey Baron, 2020

%% Generate Simulation k-space data 
% Create receiver sensitivities for simulation
N = 128; % matrix size
Nrcvr = 8; % number of receivers
R = simRcvrSens([N,N],Nrcvr);

% Cartesian: define sampling mask. Have a fully sampled region
NF_C = 28;
UR_C = 4;
mask = ones([N,N]);
mask(1:sqrt(UR_C):end,:) = 0;
mask(:,1:sqrt(UR_C):end) = 0;
mask(end/2-NF_C/2+1:end/2+NF_C/2,end/2-NF_C/2+1:end/2+NF_C/2) = 1;

% Non-Cartesian: Create projection trajectory undersampled by a factor of 4
UR_NC = 4;
[kloc,dcf] = projection(round(0.5*pi*N/UR_NC),N,2);
% Create nuFFT operator
Nop = nufftOp([N,N], kloc);

% Create image
im = phantom(N);

% Sample image. Use nufft operator for non-Cartesian, and FFT for cartesian
data_NC = Nop*(R.*im);
data_C  = mask.*ifftnc(R.*im,2);

% Add some noise
SNR = 40;
data_NC = data_NC + 1/SNR*(randn(size(data_NC)) + 1i*randn(size(data_NC)));
data_C = data_C + mask/SNR.*(randn(size(data_C)) + 1i*randn(size(data_C)));


%% Determine receiver sensitivities for Cartesian
im_rc_C = padcrop(padcrop(data_C,[NF_C NF_C]), [N N]);
im_rc_C = fftnc(im_rc_C,2);
Rop_C = rcvrOp(im_rc_C,true,NF_C*[1 1]);


%% Determine receiver sensitivities for non-Cartesian
% Projection traj is fully sampled up to radius/undersampling rate, so we
% can use that to calibrate the SENSE maps. 
% We'll solve argmin(x) || Nop*x - data ||_2^2 for each receiver channel
% Get fully sampled subset of data
inds = sqrt(sum(kloc.^2,2))<=0.5/UR_NC;
kloc_rc = kloc(inds,:);
dcf_rc = dcf(inds);
data_rc = data_NC(inds,:);

% Get starting point for iterations using regridding
NFFT = nufftOp([N,N], kloc_rc, dcf_rc);
x0 = NFFT'*data_rc;

% Solve argmin(x) || Nop*x - data ||_2^2
NFFT.dcf = [];
clear opt
opt.daNFull = size(data_rc);
opt.imNFull = size(x0);
opFunc = @(x,transp) mrSampFunc(x,transp,NFFT,[],opt);
Niterations = 3;
im_rc_NC = lsqr(opFunc,data_rc(:),[],Niterations,[],[],x0(:));
im_rc_NC = reshape(im_rc_NC,opt.imNFull);

% Compute receiver sensitivities
Nfull = ceil(N/UR_NC);
Rop_NC = rcvrOp(im_rc_NC,true,Nfull*[1 1]);


%% Solve SENSE for Cartesian, using early stopping for regularization
% Solve argmin(x) || mask*iFFT*R*x - data ||_2^2
clear opt
opt.daNFull = size(data_C);
opt.imNFull = [N N];
opFunc = @(x,transp) mrSampFunc(x,transp,mask,Rop_C,opt);
Niterations = 5;
im_SENSE_earlyStop_C = lsqr(opFunc,data_C(:),[],Niterations);
im_SENSE_earlyStop_C = reshape(im_SENSE_earlyStop_C,opt.imNFull);

%% Solve SENSE for non-Cartesian, using early stopping for regularization
% Get starting point for iterations using regridding
NFFT = nufftOp([N,N], kloc, dcf);
x0 = Rop_NC'*(NFFT'*data_NC);
% Solve argmin(x) || nuFFT*R*x - data ||_2^2
clear opt
opt.daNFull = size(data_NC);
opt.imNFull = [N N];
NFFT.dcf = [];
opFunc = @(x,transp) mrSampFunc(x,transp,NFFT,Rop_NC,opt);
Niterations = 15;
im_SENSE_earlyStop_NC = lsqr(opFunc,data_NC(:),[],Niterations,[],[],x0(:));
im_SENSE_earlyStop_NC = reshape(im_SENSE_earlyStop_NC,opt.imNFull);

%% Solve SENSE for Cartesian, using Tikhonov regularization
% Solve argmin(x) || mask*iFFT*R*x - data ||_2^2
x0 = mrSampFunc(data_C,'transp',mask,Rop_C);
clear opt
opt.daNFull = size(data_C);
opt.imNFull = [N N];
opt.tikReg = 0.02*sqrt(norm(x0(:),2)); % Depends on SNR and traj...
opFunc = @(x,transp) mrSampFunc(x,transp,mask,Rop_C,opt);
Niterations = 100;
b = [data_C(:); zeros(numel(x0),1)];
im_SENSE_tik_C = lsqr(opFunc,b,[],Niterations,[],[],x0(:));
im_SENSE_tik_C = reshape(im_SENSE_tik_C,opt.imNFull);


%% Solve SENSE for non-Cartesian, using early stopping for regularization
% Get starting point for iterations using regridding
NFFT = nufftOp([N,N], kloc, dcf);
x0 = Rop_NC'*(NFFT'*data_NC);
% Solve argmin(x) || nuFFT*R*x - data ||_2^2
clear opt
opt.daNFull = size(data_NC);
opt.imNFull = [N N];
opt.tikReg = 0.1*sqrt(norm(x0(:),2)); % Depends on SNR and traj...
NFFT.dcf = [];
opFunc = @(x,transp) mrSampFunc(x,transp,NFFT,Rop_NC,opt);
Niterations = 100;
b = [data_NC(:); zeros(numel(x0),1)];
im_SENSE_tik_NC = lsqr(opFunc,b,[],Niterations,[],[],x0(:));
im_SENSE_tik_NC = reshape(im_SENSE_tik_NC,opt.imNFull);

% Use optimal dcf from https://doi.org/10.1002/mrm.26928 to speed up convergence
NFFT.dcf = dcf.^0.25;
opFunc = @(x,transp) mrSampFunc(x,transp,NFFT,Rop_NC,opt);
b = [reshape(NFFT.dcf.*data_NC,[],1); zeros(numel(x0),1)];
im_SENSE_tik_NC_optdcf = lsqr(opFunc,b,[],Niterations,[],[],x0(:));
im_SENSE_tik_NC_optdcf = reshape(im_SENSE_tik_NC_optdcf,opt.imNFull);


%% Plotting
figure; 
subplot(2,4,[1,5]); imagesc(im); axis('image'); axis('off'); colormap('gray'); title('Ground truth'); ylabel('Cartesian')
subplot(2,4,2); imagesc(mask); axis('image'); axis('off'); colormap('gray'); title('Sampling Pattern');
subplot(2,4,3); imagesc(abs(im_SENSE_earlyStop_C)); axis('image'); axis('off'); colormap('gray'); title('regSENSE, Early Stop');
subplot(2,4,4); imagesc(abs(im_SENSE_tik_C)); axis('image'); axis('off'); colormap('gray'); title('regSENSE, Tikhonov');

subplot(2,4,6); plot(kloc(:,1),kloc(:,2),'.'); xlabel('kx'); ylabel('ky'); axis('square'); axis('off'); 
subplot(2,4,7); imagesc(abs(im_SENSE_earlyStop_NC)); axis('image'); axis('off'); colormap('gray'); 
subplot(2,4,8); imagesc(abs(im_SENSE_tik_NC_optdcf)); axis('image'); axis('off'); colormap('gray'); 



