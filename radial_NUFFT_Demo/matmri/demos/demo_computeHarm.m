%% demo for computation of trajectory from raw field probe data
%
%  This demo provides an example of computation of spherical harmonic
%  representation of phase over time based on measurements from field
%  probes. 
%
%  A standard least squares solution that also considers concomitant
%  gradient fields can be performed (see Vannesjo et al,
%  https://doi.org/10.1002/mrm.25841), or an algorithm that also mitigates
%  errors from probes placed far from isocenter can be used (manuscript in
%  preparation)
%
%  (c) Corey Baron, 2022
%

% Load field probe data collected on the 7 Tesla MRI at the CFMM for a
% diffusion weighted scan.
% Description of variables:
%   probePositions: position of each of the probes [m]
%   dataRaw:        time-series raw data measured by each field probe. The
%                   spherical harmonics are computed from its phase. 
%   dt:             time spacing between each raw data point [s]
%   fieldOffsets:   local field offset at each probe [T]. Measured using a
%                   calibration scan before the scan of interest
%   gammaProbes:    gyromagnetic ratio for probes. These ones used fluorine
%   gammaMRI:       gyromagnetic ratio for sample (typically proton)
%   B0:             magnetic field [T]
load data_raw_probes_R2spiral.mat

% Perform a ordinary least squares solution. 
opt.fitOrder = 2;          % Fit to 2nd order spherical harmonics
opt.fitDoSteps = false;    % Stepwise correction
opt.fitDoWeights = false;  % Weighted least square correction (based on probe distance from isocenter)
opt.fitDoResAdj = false;   % Weighted residual correction
[phs_spha_ls, phs_conc_ls] = harmonicsFromRaw(probePositions,dataRaw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt);

% Perform a solution that mitigates errors from distant probes
opt.fitDoSteps = true;    % Stepwise correction
opt.fitDoWeights = true;  % Weighted least square correction (based on probe distance from isocenter)
opt.fitDoResAdj = true;   % Weighted residual correction
[phs_spha_co, phs_conc_co] = harmonicsFromRaw(probePositions,dataRaw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt);

% Plotting
figure;
NP = size(phs_spha_ls,2);
nr = floor(sqrt(NP));
nc = ceil(NP/nr);
t = dt*(1:size(phs_spha_ls,1))*1000; % [ms]
ylabelsAll = {'rad', 'rad/m', 'rad/m^2', 'rad/m^3'};
for np = 1:NP
    subplot(nr,nc,np);
    plot(t,phs_spha_ls(:,np));
    hold('all')
    plot(t,phs_spha_co(:,np));
    xlabel('time (ms)')
    if np > 9
        ylabel(ylabelsAll{4});
    elseif np > 4
        ylabel(ylabelsAll{3});
    elseif np > 1
        ylabel(ylabelsAll{2});
    else 
        ylabel(ylabelsAll{1});
    end
end
legend('least squares','corrected')





