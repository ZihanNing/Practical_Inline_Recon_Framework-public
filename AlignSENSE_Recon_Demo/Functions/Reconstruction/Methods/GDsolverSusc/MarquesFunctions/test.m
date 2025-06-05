

clc; close all; clear all;
addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath('.'))

% Add DISORDER / data repository
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'))
addpath(genpath('/home/ybr19/Software/B0Modelling'))

%% Play with these parameters
rotPar = [ 45 0 0 ];%[theta_x, theta_y, theta_z] In degrees
filterKernel = 1; 

%% Kernel creation
%%%K-space
createInKspace=1;
D1 = create_dipole_kernel([0 0 1], 3*[1 1 1], [150 150 150], createInKspace, rotPar);
if filterKernel; D1 = bsxfun( @times, D1 , (buildFilter([150 150 150],'tukeyIso',.9,[],.5) ) );end
D1 = fftshift(ifftn(D1));

%%% Image space
createInKspace=0;
D2 = create_dipole_kernel([0 0 1], 3*[1 1 1], [150 150 150], createInKspace, rotPar);
if filterKernel; D2 = bsxfun( @times, D2 , (buildFilter([150 150 150],'tukeyIso',.9,[],.5) ) );end
D2 = fftshift(ifftn(D2));

%%% Plot
plotND([],real( cat(4,D1,D2) ),[-1 1]*1e-5,[],1,[],[],{'Dipole created in k-space';'Dipole created in  image space'});
if fileterKernel; title('With filtering k-space kernel');else; title('Without filtering k-space kernel');end

%plotND([],fftshift(buildFilter([150 150 150],'tukeyIso',1,[],.3)),[],[],0);

