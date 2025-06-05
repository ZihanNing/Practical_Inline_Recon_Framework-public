%%%%%%%%%%%%%%%%%%%%% test phase unwrap
% parameters
N = 512;
ampPhase = 20;
noise = 0.5;

[x, y] = meshgrid(linspace(-1,1,N));

%%%%% (1) unweighted case
% original unwrapped phase
phi = exp(-(x.*x+y.*y)/2/0.2^2) * ampPhase + (x + y) * ampPhase/2;

% wrapped phase
psi = wrapToPi(phi + randn(size(phi))*noise);

sub=[1 1];

psiS=psi(1:sub(1):end,1:sub(2):end);
phiS=phi(1:sub(1):end,1:sub(2):end);

% unweighted case
abc = tic;
phi2 = phase_unwrap(psiS);
disp(sprintf('Unweighted phase unwrap of a %dx%d image takes %f seconds', N, N, toc(abc)));

ld=2;
we=[1 1];%sub;%[1 1];%sqrt(sub);%[1 1];%sqrt(sub);%[1 1];

phiw=exp(1i*psiS);
%phiw=phiw(:,1);
%we=we(1);
phi2o=CNCGUnwrapping(phiw,ld,we);

%%%%%HEREHEREHERE---DISPLAY WITH RANGES, TURN BACK LD TO 1, CHECK WHY
%%%%%DYNAMIC RANGE IS CHANGING WITH FOURIER TRANSFORM

% show the images
%close all;
figure;
subplot(2,2,1);
imC=imagesc(phiS); title('Original phase');
rng=[min(phiS(:)) max(phiS(:))];

subplot(2,2,2);
imagesc(psiS,rng); title('Wrapped phase with noise');

%subplot(2,2,3);
%imagesc(ones(N)); title('Weight');

subplot(2,2,3);
imagesc(phi2o,rng); title('Unwrapped phase ours');

subplot(2,2,4);
imagesc(phi2,rng); title('Unwrapped phase reference');

%return

%%%%% (2) now test the weighted case
weight = ones(N);
xregion = floor(N/4):floor(N/2);
yregion = floor(N/4):floor(N/2);
weight(yregion, xregion) = 0;

% change the zero-weighted region to noise only
psi3 = psi;
psi3(yregion, xregion) = randn([length(yregion), length(xregion)]);

psi3S=psi3(1:sub(1):end,1:sub(2):end);
weightS=weight(1:sub(1):end,1:sub(2):end);

phiw=exp(1i*psi3S);
phi3o=CNCGUnwrapping(phiw,ld,we);

% now unwrap
bac = tic;
phi3 = phase_unwrap(psi3S, weightS);
%phi3 = phase_unwrap(psi3);
disp(sprintf('Weighted phase unwrap of a %dx%d image takes %f seconds', N, N, toc(bac)));

% show the images
figure;
subplot(2,2,1);
imagesc(phiS); title('Original phase');
rng=[min(phiS(:)) max(phiS(:))];


subplot(2,2,2);
imagesc(psi3S,rng); title('Wrapped phase with noise');

%subplot(2,2,3);
%imagesc(weight); title('Weight');

subplot(2,2,3);
imagesc(phi3o,rng); title('Unwrapped phase ours');

subplot(2,2,4);
imagesc(phi3,rng); title('Unwrapped phase');

%%%%% (3) test the weighted case (with noise in the border)
weight4 = zeros(N)+0.01;
xregion = floor(N/5):floor(4*N/5);
yregion = floor(N/5):floor(4*N/5);
weight4(yregion, xregion) = 1;

% change the zero-weighted region to noise only
psi4 = randn(size(psi));
psi4(yregion, xregion) = psi(yregion, xregion);

psi4S=psi4(1:sub(1):end,1:sub(2):end);
weight4S=weight4(1:sub(1):end,1:sub(2):end);

phiw=exp(1i*psi4S);
phi4o=CNCGUnwrapping(phiw,ld,we);

% now unwrap
acb = tic;
phi4 = phase_unwrap(psi4S, weight4S);
%phi4 = phase_unwrap(psi4);
disp(sprintf('Weighted phase unwrap of a %dx%d image takes %f seconds', N, N, toc(acb)));

% show the images
figure;
subplot(2,2,1);
imagesc(phi); title('Original phase');
rng=[min(phiS(:)) max(phiS(:))];

subplot(2,2,2);
imagesc(psi4S,rng); title('Wrapped phase with noise');

subplot(2,2,3);
imagesc(phi4o,rng); title('Unwrapped phase ours');

subplot(2,2,4);
imagesc(phi4,rng); title('Unwrapped phase');
