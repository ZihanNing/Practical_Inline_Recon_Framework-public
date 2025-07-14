%% sampHighOrder demo
%
%  Try the regSENSE demo before this one. Here, we include spatially
%  varying phase from eddy currents and B0 inhomogeniety via sampHighOrder.
%
%  sampHighOrder is a very powerful function to execute the forward model
%  of MRI sampling (or its adjoint) with minimal assumptions, using a brute
%  force execution of the summations in the MRI signal equation.
%  Execution times are kept reasonable using precomputation and GPU (at the
%  expense of requiring significant memory).
%
%  (c) Corey Baron, 2020

%% Generate Simulation k-space data 
% Create receiver sensitivities for simulation
UR = 2; % Undersampling rate. Instructive to try different values!
N = 96; % matrix size
Nrcvr = 8; % number of receivers
R = simRcvrSens([N,N],Nrcvr);

% Create ideal EPI trajectory undersampled by a factor of UR
kx = repmat(linspace(-1,1,N+1), [N/UR, 1]);
kx = kx(:,1:end-1);
kx(2:2:end,:) = flip(kx(2:2:end,:), 2); % Alternating direction for EPI
ky = repmat(linspace(-1,1,N/UR+1)', [1, N]);
ky = ky(1:end-1,:);
% Convert kloc to units of rad/m. Assume res of 2 mm isotropic
res = 0.002;
kmax = pi/res;
kx = kmax*kx;
ky = kmax*ky;

% Define sampling times, assuming a 50 ms readout duration
sampTimes = reshape(linspace(0,0.05,numel(kx)), size(kx'))';

% Define image locations. Use an axial slice 3 cm from isocenter.
z0 = 0.03;
[grid.x, grid.y, grid.z] = meshgrid((-N/2:N/2-1)*res, (-N/2:N/2-1)*res, z0);

% Define some exponential eddy currents using random seeds. Ensure that the
% max phase accumulation at edge of FOV is 4*pi.
E0 = 4*pi;
maxVals = 1./[1 N*res/2*ones(1,3) (N*res/2)^2*ones(1,5) (N*res/2)^3*ones(1,7)]; % Spherical harmonic basis function value at edge of FOV
phs_spha = zeros([16,size(kx)]);
allCoeffs = zeros(numel(kx),16);
for n=1:16
    tau = (0.25 + rand)*100/1000; % between 25 and 125 ms time constant 
    eddy = 1 - exp(-sampTimes/tau);
    eddy = eddy - eddy(1);
    eddy = E0*eddy/max(eddy(:))*2*(rand-0.5) * maxVals(n);
    phs_spha(n,:) = eddy(:);
    % Get all normalized time varying coefficients for plotting later
    allCoeffs(:,n) = reshape(eddy',[],1)/maxVals(n); 
end

% Add in nominal trajectory
phs_spha(2,:) = phs_spha(2,:) + kx(:)';
phs_spha(3,:) = phs_spha(3,:) + ky(:)';

% Ignore concomitant terms for this simulation
phs_conc = zeros([4,size(kx)]);

% Define magnetic field inhomogeneity that varies along y linearly. Units
% rad/s. Assume 8pi phase accumulation for a voxel at the edge of the FOV
b0 = linspace(0,1,N)'*8*pi/max(sampTimes(:));
b0 = repmat(b0, [1 N]);

% Build sampling object
% Note: here we use a 2D array for sampTimes to keep the data format
% analogous to the Cartesian EPI trajectory. sampTimes ca also be 1D
% (with phs_spha and phs_conc 2D instead of 3D), which would be the typical
% scenario for non-Cartesian acquisitions. One can make it 1D using
% sampTimes(:) here too, and it should work fine (the data output would
% become 1D for each receiver instead of being 2D as in this example)
S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,grid);

% Create image. 
im = phantom(N);

% Sample image
data = S*(R.*im);

% Add some noise
SNR = 60;
data = data + 1/SNR*(randn(size(data)) + 1i*randn(size(data)));


%% Perform a naive recon with fft then RSOS
im_fftsos = data;
im_fftsos(2:2:end,:) = flip(im_fftsos(2:2:end,:),2);
im_fftsos = fftnc(im_fftsos);
im_fftsos = calcRSOS(im_fftsos,3);


%% Perform regularized SENSE, assuming that there's no eddy currents or field inhomogeneity
% Create receiver operator
Rop = rcvrOp(R,false);
% Zero fill data
data_rs = zeros([N,N,Nrcvr]);
data_rs(1:UR:end,:,:) = data;
data_rs((1+UR):(2*UR):end,:,:) = flip(data_rs((1+UR):(2*UR):end,:,:),2);
mask = zeros([N,N]);
mask(1:UR:end,:) = 1;
% Solve argmin(x) || mask*iFFT*R*x - data ||_2^2
clear opt
opt.daNFull = size(data_rs);
opt.imNFull = [N N];
opFunc = @(x,transp) mrSampFunc(x,transp,mask,Rop,opt);
Niterations = 5;
im_SENSE = lsqr(opFunc,data_rs(:),[],Niterations);
im_SENSE = reshape(im_SENSE,opt.imNFull);


%% Perform model based recon with full model
% Create receiver operator
Rop = rcvrOp(R,false);
% Solve argmin(x) || S*R*x - data ||_2^2
clear opt
opt.daNFull = size(data);
opt.imNFull = [N N];
opFunc = @(x,transp) mrSampFunc(x,transp,S,Rop,opt);
Niterations = 5;
im_modelBased = lsqr(opFunc,data(:),[],Niterations);
im_modelBased = reshape(im_modelBased,opt.imNFull);


%% Plotting
figure;
subplot(2,4,[1,5]); imagesc(im); axis('image'); axis('off'); colormap('gray'); title('Ground truth'); ylabel('Cartesian')
subplot(2,4,2); plot(reshape(kx',[],1), reshape(ky',[],1)); title('Trajectory'); xlabel('kx [rad/m]'); ylabel('ky [rad/m]');
subplot(2,4,3); plot(reshape(sampTimes'*1000,[],1),allCoeffs); title('Eddy current phase at edge of FOV for each spherical harmonic');  xlabel('time [ms]') ; ylabel('phase [rad]')
subplot(2,4,4); imagesc(b0/1000); title('B0 map'); axis('image'); axis('off');  c=colorbar; c.Label.String = 'rad/ms';

subplot(2,4,6); imagesc(abs(im_fftsos)); axis('image'); axis('off'); colormap('gray'); title('FFT then RSOS');
subplot(2,4,7); imagesc(abs(im_SENSE)); axis('image'); axis('off'); colormap('gray'); title('Normal regularized SENSE');
subplot(2,4,8); imagesc(abs(im_modelBased)); axis('image'); axis('off'); colormap('gray'); title('Full model based recon');






