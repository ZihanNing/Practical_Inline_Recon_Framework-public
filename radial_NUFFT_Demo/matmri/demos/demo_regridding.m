%% Regridding Demo
%
%  Regridding is performed by first defining a type-2 nuFFT operator,
%  updating it's density compensation function (dcf), and then applying the
%  adjoint of the nuFFT.
%
%  More accurate reconstructions are obtained by finishing with a few
%  conjugate gradient iterations (without a dcf)
%
%  (c) Corey Baron, 2020

%% Generate Simulation k-space data using nuFFT
% Create fully sampled projection trajectory
N = 128; % matrix size
[kloc,dcf] = projection(round(0.5*pi*N),N,2);

% Create nuFFT operator
Nop = nufftOp([N,N], kloc);

% Create image
im = phantom(N);

% Sample image
data = Nop*im;

% Add some noise
SNR = 20;
data = data + 1/SNR*(randn(size(data)) + 1i*randn(size(data)));


%% Perform regridding
% Set density compensation. 
% For more complicated traj, a numerical approach can be done via 
%   Nop = Nop.findDcf; % (slow, but works well)
Nop.dcf = dcf;

% Perform regridding. Note the aliasing in the corners is just a
% consequence of a projection trajectory
im_r = Nop'*data;


%% Use iterative method that is more accurate and does not require dcf
% Solves argmin(x) || Nop x - data ||_2^2
% Requires few iterations for fully sampled with a good initialization
Nop.dcf = [];
opt.imNFull = size(im); % This causes vectorization of image output that is required for lsqr
opFunc_a = @(x,transp) mrSampFunc(x,transp,Nop,[],opt); % Function needed for lsqr
Niterations = 3;
im_i = lsqr(opFunc_a,data,[],Niterations,[],[],im_r(:)); % Note the initial guess of im_r is not necessary, but it speeds up convergence
im_i = reshape(im_i,[N,N]);


%% Plotting
figure; 
subplot(2,2,1); plot(kloc(:,1),kloc(:,2),'.'); xlabel('kx'); ylabel('ky');
subplot(2,2,2); imagesc(im); axis('image'); colormap('gray');   title('Ground truth');
subplot(2,2,3); imagesc(abs(im_r)); axis('image'); colormap('gray'); title('Regridded');
subplot(2,2,4); imagesc(abs(im_i)); axis('image'); colormap('gray'); title('Iterative');



