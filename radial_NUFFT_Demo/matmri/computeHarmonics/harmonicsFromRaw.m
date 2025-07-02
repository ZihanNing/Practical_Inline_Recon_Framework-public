function [phs_spha, phs_conc] = harmonicsFromRaw(probe_positions,probe_raw,dt,fieldOffsets,gammaProbes,gammaMRI,B0,opt)
% Compute trajectory for high-order encoding model using raw data from
% field probes.
%
% Inputs:
%   probePositions: position of each of the probes [m]
%   dataRaw:        time-series raw data measured by each field probe. The
%                   spherical harmonics are computed from its phase 
%   dt:             time spacing between each raw data point [s]
%   fieldOffsets:   local field offset at each probe [T]. Measured using a
%                   calibration scan before the scan of interest
%   gammaProbes:    gyromagnetic ratio for probes. These ones used fluorine
%   gammaMRI:       gyromagnetic ratio for sample (typically proton)
%   B0:             magnetic field [T]
%
% Outputs:
%   phs_spha:       spherical harmonic coefficients that correspond to
%                   basis functions defined in basisFuncHarm.m
%   phs_conc:       concomitant field coefficients that correspond to basis
%                   functions defined in basisFuncConc.m
%
% (c) Corey Baron and Paul Dubovan, 2021
%

if nargin<8 || isempty(opt) || ~isfield(opt,'fitOrder')
    % Determine below, based on number of probes
    opt.fitOrder = [];
end
if ~isfield(opt,'fitDoSteps') || isempty(opt.fitDoSteps)
    % Whether to fit spatial orders in step-wise fashion
    opt.fitDoSteps = true;
end
if ~isfield(opt,'fitDoWeights') || isempty(opt.fitDoWeights)
    % Whether to use weighted least squares with weights based on distance
    % of probes from isocenter
    opt.fitDoWeights = true;
end
if ~isfield(opt,'fitDoResAdj') || isempty(opt.fitDoResAdj)
    % Whether to adjust phase of probes that are far from isocenter, based
    % on residual of initial fit
    opt.fitDoResAdj = true;
end
if ~isfield(opt,'fitResAdjNit') || isempty(opt.fitResAdjNit)
    % Number of iterations for fitDoResAdj option. Heuristically
    % determined. Should be 2 or 3.
    opt.fitResAdjNit = 3;
end
if ~isfield(opt,'coilParams') || isempty(opt.coilParams) || isempty(opt.coilParams.alph) 
    % Assume a symmetric gradient coil
    opt.coilParams.alph = 0.5;
end
if ~isfield(opt,'fitInds') 
    % Can be used to specify exactly which basis functions are fit by
    % index. Overrides "fitOrder" option.
    opt.fitInds = [];
end
if ~isfield(opt,'noConc') || isempty(opt.noConc) 
    % Ignore concomitant gradients. Should only be used for debugging
    opt.noConc = false;
end
coilParams = opt.coilParams;
fitOrder = opt.fitOrder;


% Correct for probe field offsets
times = (1:size(probe_raw,1))';
phsCor = 2*pi*(gammaProbes*dt)*fieldOffsets(:)'.*times;
probe_raw = probe_raw .* exp(-1i*phsCor);

% Remove a probe if raw probe amplitude falls below threshold (default 0.1%) 
% Probe magnitudes normalized based on the initial point
probemag_thresh = 0.001; %default value based on testing, change if needed
probe_raw_norm = abs(probe_raw)./abs(probe_raw(1,:,:,:,:));
probe_index = 1:size(probe_positions,1);
for nprobe = 1:size(probe_positions,1)
    if any(probe_raw_norm(:,nprobe,1,1) < probemag_thresh)
        remove_ind = nprobe;
        probe_index(probe_index == remove_ind) = [];
        warning('Probe %d removed from fit due to raw probe amplitude falling below threshold\n', nprobe)
    end            
end    
probe_positions = probe_positions(probe_index,:);
probe_raw = probe_raw(:,probe_index,:,:);
    
% Determine number of probes
Nprobe = size(probe_positions,1);

% Set defaults OR verify the designated fit order has the minimum number of probes needed.
% If below the minimum probe requirement due to probe removal, set to appropriate fit
% order based on number of available probes.
if isempty(opt.fitInds)
    if isempty(fitOrder)
        if Nprobe >= 16
            fitOrder = 3;
        elseif Nprobe >= 9
            fitOrder = 2;
        else
            fitOrder = 1;
        end
    else
        if fitOrder == 3
            if Nprobe < 16 & Nprobe >= 9
                fitOrder = 2;
                warning('Fit order set to %d since not enough probes for desired order selection', fitOrder)
            elseif Nprobe < 9
                fitOrder = 1;
                warning('Fit order set to %d since not enough probes for desired order selection', fitOrder)
            end
        elseif fitOrder == 2
            if Nprobe < 9
                fitOrder = 1;
                warning('Fit order set to %d since not enough probes for desired order selection', fitOrder)
            end
        end  
    end  
    switch fitOrder
        case 1
            opt.fitInds = 1:4;
        case 2
            opt.fitInds = 1:9;
        case 3
            opt.fitInds = 1:16;
        otherwise
            error('unknown value for opt.fitOrder')
    end
else
    if max(opt.fitInds) >= 16
        fitOrder = 3;
    elseif max(opt.fitInds) >= 9
        fitOrder = 2;
    else
        fitOrder = 1;
    end
    opt.fitInds = sort(opt.fitInds);
end
Nbasis = length(opt.fitInds);

% Parameter checking
if (Nbasis >= Nprobe) && opt.fitDoResAdj
    warning('opt.fitDoResAdj not recommended when number of probes does not exceed number of fitted basis functions.')
end

% Prep. 
phsRaw = unwrap(angle(probe_raw));
magRaw = abs(probe_raw); 

% Have the phase start at 0 for all probes
phsRaw = phsRaw - phsRaw(1,:,:,:);

% Compute squared distance from isocenter, to be used for weighting
W = [];
if opt.fitDoWeights
    % Using net distance from isocenter
    W = sum(probe_positions.^2,2);
        % Other metrics to potentially try instead:
        % Using sum of all distances from isocenter
        %W = sum(abs(probe_positions),2).^2;
        % Using max z-distance from isocenter
        %W = probe_positions(:,3).^2;

    W = 1./W; % Closer probes should have larger weights

    % Normalize weights
    W = W/max(W);
end

% Create matrix for maxwell terms basis functions
Nc = 4;
B = zeros(Nprobe, Nc, 'like', phsRaw);
for l=1:Nc
    B(:,l) = basisFuncConc(probe_positions(:,1),probe_positions(:,2),probe_positions(:,3),l);
end

% Iterate between computing SphHarm and phase from conc grads
Nit = 3;
if ~opt.fitDoResAdj
    opt.fitResAdjNit = 1;
end
phs_spha = zeros(size(phsRaw,1),Nbasis,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
phs_conc = zeros(size(phsRaw,1),Nc,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
for nv = 1:size(phsRaw(:,:,:),3)
    phsRaw_a = phsRaw(:,:,nv);
    if ~opt.fitDoSteps
        % Fit all the orders of spherical harmonics simultaneously
        norderAll = fitOrder;
    else
        % Incrementally increase order, because solving for all terms
        % simultanously can cause error in the lower order terms.
        % NB: shouldn't start with order 0, because if the probes are net
        % offcenter it will mess up the B0 estimation.
        norderAll = 1:fitOrder;
    end
    % Loop through orders (required for "opt.fitDoSteps" option)
    for ncorr = 1:opt.fitResAdjNit
        phsRaw_conc = 0;
        % Initialize phase that is subtracted from raw phase after fitting
        % lower orders 
        phs_pre = 0; 
        for norder_ind = 1:length(norderAll)
            norder = norderAll(norder_ind);
            % Set matrix equation for normal spherical harmonic expansion of trajectory
            if norder == 3
                Nla = find(opt.fitInds>=10,1,'first');
                Nlb = find(opt.fitInds<=16,1,'last');
            elseif norder == 2
                Nla = find(opt.fitInds>=5,1,'first');
                Nlb = find(opt.fitInds<=9,1,'last');
            elseif norder == 1
                Nla = 1;
                Nlb = find(opt.fitInds<=4,1,'last');
            else
                error('Must choose order between 1 and 3')
            end
            if ~opt.fitDoSteps
                % Here we need to do all orders at once, so we start from the
                % first one (i.e., Nla=1) and go up to the one specified via
                % norder above (i.e., Nl)
                Nla = 1;
            end
            A = zeros(Nprobe, Nlb-Nla+1, 'like', phsRaw);
            for l=Nla:Nlb
                basisIndex = opt.fitInds(l);
                A(:,l-Nla+1) = basisFuncHarm(probe_positions(:,1),probe_positions(:,2),probe_positions(:,3),basisIndex);
            end
            phs_in = phsRaw_a.' - phs_pre;
            if norder_ind == 1
                % Also fit concomitant gradient terms
                [phs_spha_a, phs_conc_a, phsRaw_conc] = doFit(A,B,phs_in,W,0,Nit,dt,gammaProbes,coilParams,B0,opt.noConc);
            else
                phs_spha_a = doFit(A,B,phs_in,W,phsRaw_conc);
            end
            phs_pre = phs_pre + A*(phs_spha_a);
            resid_ord = phsRaw_a.' - phsRaw_conc - phs_pre;
            % Save result into output
            phs_spha(:,Nla:Nlb,nv) = phs_spha_a.';
        end
        if opt.fitDoResAdj
            % Remove residual phase (large residuals are created by
            % step-wise approach). These residuals represent errors from
            % non-linearity / inhomogeneity of fields (gradient or, B0, or
            % other eddy current modes)
            % Perform a larger correction for probes at a further distance
            W2 = 1./W;
            W2 = W2/max(W2);
            phsRaw_a = phsRaw_a - W2'.*(resid_ord');
        end
    end
    if ~isempty(phs_conc_a)
        phs_conc(:,:,nv) = phs_conc_a;
    end
end
%figure; plot(resid)


% Format output and scale phase by gammas
phs_spha = phs_spha * gammaMRI / gammaProbes;
phs_conc = phs_conc * gammaMRI / gammaProbes;


end

function [phs_spha_a, phs_conc_a, phsRaw_conc] = doFit(A,B,phs_in,W,phsRaw_conc,Nit,dt,gammaProbes,coilParams,B0,noConc)
    phs_conc_a = [];
    if (nargin < 7) || noConc
        doConc = 0;
        Nit = 1;
    else
        doConc = 1;
        phsRaw_conc = 0;
    end
    if ~isempty(W)
        pA = pinv(diag(W)*A);
    else
        pA = pinv(A);
    end
    for n=1:Nit
        if ~isempty(W)
            phs_spha_a = pA*(W(:).*(phs_in - phsRaw_conc));
        else
            phs_spha_a = pA*(phs_in - phsRaw_conc);
        end
        if doConc
            linInds = 2:4;
            phs_conc_a = compMaxwellPhase(phs_spha_a(linInds,:).',dt,gammaProbes,coilParams,B0);
            phsRaw_conc = B*(phs_conc_a.');
        end
    end
end


