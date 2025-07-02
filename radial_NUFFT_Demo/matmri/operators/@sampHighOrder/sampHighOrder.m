classdef sampHighOrder
% Perform MRI sampling using direct summation of complex exponentials.
% Supports B0 map, time-varying spherical harmonic distributions of phase,
% and time-varying distributions of phase that are typical for concomitant
% gradients. Supports 1D, 2D, or 3D.
%
% Usage: 
%   Define sampHighOrder object using 
%       S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,imMask,useGPU,useSingle,useInterp,svdThresh,subFact)
%   Sample image using data = S*image. x can have more dimensions than
%       b0, but summations are only performed over first two dims (i.e., 2D
%       imaging is assumed). Fairly high memory overhead due to 
%       precomputations, but quite fast.
%   Perform adjoint of sampling using S'*data. Required for iterative
%       solvers like lsqr.
%
% Inputs: 
%   b0: B0 map in units of rad/s. Caution: make sure orientation
%       of b0 is consistent with phs_grid.x and phs_grid.y.
%       Can provide empty matrix to ignore B0 effects (not recommended).
%   sampTimes: sampling times in units of s. Can be multi-dimensional.
%   phs_spha:  [Ncoeff_spha x size(sampTimes)] array of coefficients for
%       spherical harmonic phase. Units are rad/<spatial>, where <spatial>
%       can be unitless (DC), m (normal k-space units), m^2, etc, depending
%       on the coefficient.
%   phs_conc:  [Ncoeff_conc x size(sampTimes)] array of coefficients for
%       concomitant grad phase. Units are rad/<spatial>, where <spatial>
%       can be unitless (DC), m (normal k-space units), m^2, etc, depending
%       on the coefficient.
%   phs_grid: struct with x, y, and z positions for each voxel in b0. 
%       Should be created using ndgrid or meshgrid. If created with
%       meshgrid, phs_grid.x will have positions varying along 2nd dim
%       (with ndgrid, variation along 1st dim), and vice versa for
%       phs_grid.y. Size of b0 map must be same size as each field of phs_grid.
%           phs_grid.x: value of x at each position in units of m. 
%           phs_grid.y: value of y at each position in units of m. 
%           phs_grid.z: value of z at each position in units of m. 
%   imMask: mask denoted expected location of signal. b0map and 2nd
%       order+ spherical harmonic phase terms are set to 0 outside the
%       mask. This improves the performance of interpolation (useInterp
%       option), because it removes regions with quickly varying phase that
%       are irrevelant to image recon (since there's no signal there).
%   useGPU: whether to use GPU. Default = true.
%   useSingle: if true, use single instead of double precision. Default =
%       false. 
%   useInterp: uses an interpolated approach. Always faster than direct
%       approach on CPU (i.e., when useGPU = 1). Not much benefit for GPU,
%       unless not enough memory for direct approach. 
%		See Wilm et al DOI: 10.1109/TMI.2012.2190991
%   svdThresh (default = 0.05): trades off accuracy with speed when
%       useInterp=1. Decrease to improve accuracy.
%   subFact (default = 5): another trade-off for accuracy and speed for 
%       useInterp=1. Decrease for accuracy. 5 or less should have
%       negligible error. Min val = 1.
%
%
%   Within the class, basis functions for each index of phs_spha and
%   phs_conc are computed using basisFuncHarm.m and basisFuncConc.m,
%   respectively. 
%
%   (c) Corey Baron 2020

	properties
        useSingle = 0;
		useGPU = 1;
		useInterp = 0; 
        useSegmented = 0;    % Whether to use segmented encoding matrix. May be necessary for large data size
                          %   1: use segmented, but still precompute array and save in RAM
                          %   2: use segmented with no precomputation at all. Slowest option, but least memory.
        svdThresh = 0.05; % Threshold for svd used in interpolated method
        subFact = 5;      % Factor to subsample in time by for interpolated approach
        subFactSpc = 1;   % Factor to subsample in space by for interpolated approach
	end

	properties (SetAccess = protected)
		adjoint = 0;
        NDim = [];  % Number of object domain dimensions
        NDimk = []; % Number of k-space dimensions
		b0 = []; % rad/sec
        b0mask = [];
        phs_spha = []; % 0th, 1st, 2nd and 3rd order spherical harmonic terms for phase over time. Must have dimensions 16 x size(sampTimes). 
        phs_conc = []; % conc grad terms for phase over time. Must have dimensions 4 x size(sampTimes). 
        phs_grid = []; % Struct with fields phs_grid.x, phs_grid.x, phs_grid.x. Each must have dimensions equivalent to b0. MUST be in magnet frame.
		sampTimes = []; % sec
		kbase = []; % spatial part of spherical harmonics that is multiplied with phs_spha or phs_conc terms
        phiDiv = []; % temporal derivative of kbase. Can be used to find global delays. See DOI: 10.1002/mrm.29460
		traj = [];		  % Precomputed values for interpolated approach
        svdSpace = [];    % Precomputed values for interpolated approach
        svdTime = [];     % Precomputed values for interpolated approach
        phsShft = [];
		kSize = [];
		imSize = [];
        trajFromRaw = [];
        useSegmentedNdiv = [];
        ksphaDiv = [];
        kconcDiv = [];
        SegNpnts = [];     % Seg* variables are for segmented encoding matrix with precomputation
        SegNpntsAdj = [];
        SegPhs = [];
	end

	methods

		function obj = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,b0mask,useGPU,useSingle,useInterp,svdThresh,subFacts,useSegmented)
			if nargin == 0
				obj.tests;
				return;
			end
            if nargin>5
				obj.b0mask = b0mask;
            end
            if nargin<6 || isempty(useGPU)
                % Determine options based on memory available
                % Precomputations using GPU requires about 3.3 times numel(b0)*numel(sampTimes)*8 bytes of memory 
                % After this calc, sampHighOrder object requires numel(b0)*numel(sampTimes)*8 bytes
                % If the calc is successful, there is enough working memory to run
                % multiplications with sampHighOrder object.
                % To be safe, we allow for 4 times numel(b0)*numel(sampTimes)*8
                if ~(gpuDeviceCount>0)
                    obj.useGPU = 0;
                    obj.useInterp = 1;
                else
                    G = gpuDevice;
                    avMem = G.AvailableMemory;
                    numbytes = numel(phs_grid.x)*numel(sampTimes)*8;
                    if numbytes < 0.2*avMem
                        obj.useGPU = 1;
                        obj.useSingle = 0;
                        obj.useInterp = 0;
                    elseif numbytes < 0.4*avMem
                        obj.useGPU = 1;
                        obj.useSingle = 1;
                        obj.useInterp = 0;
                    else
                        obj.useGPU = 1;
                        obj.useSingle = 1;
                        obj.useInterp = 0;
                        obj.useSegmented = 1;
                        obj.useSegmentedNdiv = ceil(numbytes / (0.2*avMem)); % 0.2 factor was determine heuristically
                        if obj.useSegmentedNdiv > max(numel(b0),numel(sampTimes))/4
                            % Unlikely, but possible. Still use single,
                            % because matrix size must be massive.
                            obj.useGPU = 0;
                            obj.useSingle = 1;
                            obj.useInterp = 1;
                            obj.useSegmented = 0;
                        end
                    end
                end
            else
                if nargin>6 && ~isempty(useGPU)
                    obj.useGPU = useGPU;
                end
                if nargin>7 && ~isempty(useSingle)
                    obj.useSingle = useSingle;
                end
                if nargin>8 && ~isempty(useInterp)
                    obj.useInterp = useInterp;
                end
                if nargin>9 && ~isempty(svdThresh)
                    obj.svdThresh = svdThresh;
                end
                if nargin>10 && ~isempty(subFacts)
                    obj.subFact = subFacts(1);
                    if numel(subFacts)>1
                        obj.subFactSpc = subFacts(2);
                    end
                end
                if nargin>11 && ~isempty(useSegmented)
                    obj.useSegmented = useSegmented;
                    obj.useSegmentedNdiv = 2;
                end
            end
            if obj.subFact < 1
                obj.subFact = 1;
            end
            % Check for GPU support
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected or GPU not supported. Using CPU.')
                obj.useGPU = 0;
            end  
            obj.phs_grid = phs_grid;
			obj.imSize = size(phs_grid.x);  
            if isempty(b0)
                obj.b0 = zeros(size(phs_grid.x), 'like', phs_grid.x);
            elseif ~all(size(b0) == obj.imSize)
                error('Size mismatch between b0 and phs_grid.x');
            else
                obj.b0 = b0;
            end
			obj.sampTimes = sampTimes;
            obj.kSize = size(obj.sampTimes);
            if (length(obj.imSize) == 2) && (obj.imSize(2) == 1)
                obj.imSize = obj.imSize(1);
            end
            if (length(obj.kSize) == 2) && (obj.kSize(2) == 1)
                obj.kSize = obj.kSize(1);
            end
            obj.NDim  = length(obj.imSize);
            obj.NDimk = length(obj.kSize);
            if (obj.NDim > 3) || (obj.NDimk > 3)
                error('Only up to 3D in object- or sampling-domain allowed')
            end
			obj.phs_spha = phs_spha;
			sz_a = size(phs_spha);
			if any(sz_a(2:end) ~= obj.kSize)
				error('Dimension mismatch between phs_spha and sampTimes')
			end
            % Check if gpu is possible
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected. Using CPU. To disable this warning, set option useGPU to 0')
                obj.useGPU = 0;
            end
			obj.phs_conc = phs_conc;
			% Use single precision if requested
			if obj.useSingle
                obj.sampTimes = single(obj.sampTimes);
				obj.phs_spha = single(obj.phs_spha);
				obj.phs_conc = single(obj.phs_conc);
				obj.b0 = single(obj.b0);
				obj.phs_grid.x = single(obj.phs_grid.x);
				obj.phs_grid.y = single(obj.phs_grid.y);
				obj.phs_grid.z = single(obj.phs_grid.z);
                obj.b0mask = single(obj.b0mask);
			end
            % Move variables to GPU
			if obj.useGPU 
				obj.sampTimes = gpuArray(obj.sampTimes);
				obj.phs_spha = gpuArray(obj.phs_spha);
				obj.phs_conc = gpuArray(obj.phs_conc);
				obj.b0 = gpuArray(obj.b0);
				obj.phs_grid.x = gpuArray(obj.phs_grid.x);
				obj.phs_grid.y = gpuArray(obj.phs_grid.y);
				obj.phs_grid.z = gpuArray(obj.phs_grid.z);
                obj.b0mask = gpuArray(obj.b0mask);
			end
			if obj.useInterp
                [obj.svdTime,obj.svdSpace,obj.traj,obj.phsShft] = prepForInterp(obj);
                if obj.useGPU
                    obj.svdTime = gpuArray(obj.svdTime);
                    obj.svdSpace = gpuArray(obj.svdSpace);
                    obj.phsShft = gpuArray(obj.phsShft);
                end
            elseif ~obj.useSegmented
				obj.kbase = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes);
            elseif obj.useSegmented 
                obj.SegNpnts = ceil(prod(obj.imSize)/obj.useSegmentedNdiv);
                obj.SegNpntsAdj = ceil(prod(obj.kSize)/obj.useSegmentedNdiv);
                if obj.useSegmented == 1
                    obj.SegPhs = prepForSegmented(obj);
                end
			end
		end

		function y = mtimes(obj,x)			
            if obj.useSingle
                if (isa(x,'gpuArray') && ~isaUnderlying(x,'single')) || ~isa(x,'single')
                    %warning('useSingle specified, but input is not single. Forcing to be single.')
                    x = single(x);
                end
            end

            % Choose method (all should give the same result, but with
            % different trade-offs in terms of speed, memory usage, and
            % accuracy (only obj.useInterp can reduce accuracy)
            if obj.useInterp
                y = useInterpWorker(obj,x);
            elseif obj.useSegmented
                y = useSegmentedWorker(obj,x);
            else
                y = useDirectWorker(obj,x);
            end
        end
		
		function y = useInterpWorker(obj,x)
            % Use nufft's and interpolation
            szx = [size(x), 1, 1];
            if obj.useGPU && ~isa(x, 'gpuArray')
                x_a = gpuArray(x);
            else
                x_a = x;
            end
            nbins = size(obj.svdTime,2);
            if obj.adjoint 
                NRep = numel(x_a)/prod(obj.kSize);
                x_a = reshape(x_a, [obj.kSize, NRep]);
                y = x_a.*conj(reshape(obj.svdTime, [obj.kSize,1,nbins])); 
                y = conj(obj.phsShft).*y;
                y = obj.traj'*y;
                y = y.*conj(reshape(obj.svdSpace, [obj.imSize,1,nbins]));
                y = sum(y,ndims(y));
                y = reshape(y, [obj.imSize, szx(obj.NDimk+1:end)]);
            else
                if ~isempty(obj.ksphaDiv)
                    error('phiDiv not yet implemented for interpolated approach')
                end
                NRep = numel(x_a)/prod(obj.imSize);
                x_a = reshape(x_a, [obj.imSize, NRep]);  
                y = x_a.*reshape(obj.svdSpace, [obj.imSize,1,nbins]);
                y = obj.traj*y;
                y = obj.phsShft.*reshape(y, [obj.kSize,NRep,nbins]);
                y = y.*reshape(obj.svdTime, [obj.kSize,1,nbins]);
                y = sum(y,ndims(y));
                y = reshape(y, [obj.kSize, szx(obj.NDim+1:end)]);
            end
            if ~isa(x,'gpuArray')
                y = gather(y);
            end
        end
        
        function y = useDirectWorker(obj,x)
            % Use direct model
            szx = [size(x), 1, 1];
            if obj.useGPU && ~isa(x, 'gpuArray')
                y = gpuArray(x);
            else
                y = x;
            end
            if obj.adjoint
                y = reshape(x, prod(obj.kSize), []);
                y = exp(-1i*obj.kbase)*y;
                y = reshape(y, [obj.imSize, szx(obj.NDimk+1:end)]);
            else
                y = reshape(y, prod(obj.imSize), []);
                if isempty(obj.ksphaDiv)
                    y = exp(1i*obj.kbase.')*y;
                else
                    tmp = 1i*exp(1i*obj.kbase).*prepForDirect(obj,obj.ksphaDiv,obj.kconcDiv,1);
                    tmp = tmp.';
                    y = tmp * y;
                    clear tmp
                end
                y = reshape(y, [obj.kSize, szx(obj.NDim+1:end)]);
            end
            % Normalization so that a cartesian Fourier transform would have
                    % the adjoint equal to the inverse
            y = y/sqrt(prod(obj.imSize));
            if ~isa(x,'gpuArray')
                y = gather(y);
            end
		end

        function y = useSegmentedWorker(obj,x)
            % Use direct model in multiple segments (requires less memory,
            % but is slower)
            szx = [size(x), 1, 1];
            if obj.useGPU && ~isa(x, 'gpuArray')
                x_a = gpuArray(x);
            else
                x_a = x;
            end
            
            y = 0;
            if obj.adjoint
                x_a = reshape(x_a, prod(obj.kSize), []);
            else
                x_a = reshape(x_a, prod(obj.imSize), []);
            end   
            for nc = 1:obj.useSegmentedNdiv
                if obj.adjoint
                    subinds1 = 1 + (nc-1)*obj.SegNpntsAdj;
                    subinds2 = min(obj.SegNpntsAdj + (nc-1)*obj.SegNpntsAdj, prod(obj.kSize));
                    if obj.useSegmented > 1
                        phs = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                            obj.sampTimes,[],[],[subinds1,subinds2]);
                    else
                        if obj.useGPU
                            phs = gpuArray(obj.SegPhs{nc,2});
                        else
                            phs = obj.SegPhs{nc,2};
                        end
                    end
                    x_sub = x_a(subinds1:subinds2,:);
                    x_sub = exp(-1i*phs)*x_sub;
                    y = y + reshape(x_sub, [obj.imSize, szx(obj.NDimk+1:end)]);
                else
                    subinds1 = 1 + (nc-1)*obj.SegNpnts;
                    subinds2 = min(obj.SegNpnts + (nc-1)*obj.SegNpnts, prod(obj.imSize));
                    if obj.useSegmented > 1
                        phs = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                            obj.sampTimes,[],[subinds1,subinds2],[]);
                    else
                        if obj.useGPU
                            phs = gpuArray(obj.SegPhs{nc,1});
                        else
                            phs = obj.SegPhs{nc,1};
                        end
                    end
                    x_sub = x_a(subinds1:subinds2,:);
                    if isempty(obj.ksphaDiv)
                        x_sub = exp(1i*phs.')*x_sub;
                    else
                        tmp = 1i*exp(1i*phs).*prepForDirect(obj,obj.ksphaDiv,obj.kconcDiv,1,[],[subinds1,subinds2],[]);
                        tmp = tmp.';
                        x_sub = tmp * x_sub;
                        clear tmp
                    end
                    y = y + reshape(x_sub, [obj.kSize, szx(obj.NDim+1:end)]);
                end
            end
            % Normalization so that a cartesian Fourier transform would have
                % the adjoint equal to the inverse
            y = y/sqrt(prod(obj.imSize));
            if ~isa(x,'gpuArray')
                y = gather(y);
            end
        end
        
        function SegPhs = prepForSegmented(obj)
            try
                % Precompute segmented encoding matrix
                SegPhs = cell(obj.useSegmentedNdiv, 2);
                for nc = 1:obj.useSegmentedNdiv
                    for nm = 1:2
                        if nm==2
                            Npnts = obj.SegNpntsAdj;
                            NpntsTot = prod(obj.kSize);
                        else
                            Npnts = obj.SegNpnts;
                            NpntsTot = prod(obj.imSize);
                        end
                        subinds1 = 1 + (nc-1)*Npnts;
                        subinds2 = min(Npnts + (nc-1)*Npnts, NpntsTot);
                        if subinds1 <= NpntsTot
                            if nm==2
                                tmp = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                                    obj.sampTimes,[],[],[subinds1,subinds2]);
                            else
                                tmp = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                                    obj.sampTimes,[],[subinds1,subinds2],[]);
                            end
                            SegPhs{nc,nm} = gather(tmp);
                        end
                    end
                end
            catch ME
                if contains(ME.message, 'memory') || strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
                    msg = 'Set sampHighOrder option useSegmented = 2 if out of memory';
                    causeException = MException('MATLAB:sampHighOrder:memory',msg);
                    ME = addCause(ME,causeException);
                end
                rethrow(ME);
            end
        end

		function kbase = prepForDirect(obj,phs_spha_a,phs_conc_a,sampTimes_a,sphaInds,subIndsSpace,subIndsTime,subIndsSpaceDirect,subIndsTimeDirect)
            if nargin<5 || isempty(sphaInds)
                sphaInds = 1:size(obj.phs_spha,1);
            end
            if nargin<6 || isempty(subIndsSpace)
                subIndsSpace = [];
            end
            if nargin<7 || isempty(subIndsTime)
                subIndsTime = [];
            end
            if nargin<8 
                subIndsSpaceDirect = [];
            end
            if nargin<9 
                subIndsTimeDirect = [];
            end
            if ~isempty(subIndsTime)
                phs_spha_a = phs_spha_a(:,subIndsTime(1):subIndsTime(2));
                phs_conc_a = phs_conc_a(:,subIndsTime(1):subIndsTime(2));
                if numel(sampTimes_a)>1
                    sampTimes_a = sampTimes_a(subIndsTime(1):subIndsTime(2));
                end
            elseif ~isempty(subIndsTimeDirect)
                phs_spha_a = phs_spha_a(:,subIndsTimeDirect);
                phs_conc_a = phs_conc_a(:,subIndsTimeDirect);
                if numel(sampTimes_a)>1
                    sampTimes_a = sampTimes_a(subIndsTimeDirect);
                end
            end
            if ~isempty(subIndsSpace) 
                if ~isempty(obj.b0mask)
                    b0mask_a = obj.b0mask(subIndsSpace(1):subIndsSpace(2));
                    b0mask_a = b0mask_a(:);
                end
                b0_a = obj.b0(subIndsSpace(1):subIndsSpace(2));
                b0_a = b0_a(:);
                x = obj.phs_grid.x;
                y = obj.phs_grid.y;
                z = obj.phs_grid.z;
            elseif ~isempty(subIndsSpaceDirect)
                if ~isempty(obj.b0mask)
                    b0mask_a = obj.b0mask(subIndsSpaceDirect);
                end
                b0_a = obj.b0(subIndsSpaceDirect);
                x = obj.phs_grid.x(subIndsSpaceDirect);
                y = obj.phs_grid.y(subIndsSpaceDirect);
                z = obj.phs_grid.z(subIndsSpaceDirect);
            else
                if ~isempty(obj.b0mask)
                    b0mask_a = obj.b0mask(:);
                end
                b0_a = obj.b0(:);
                x = obj.phs_grid.x;
                y = obj.phs_grid.y;
                z = obj.phs_grid.z;
            end
			% Vectorize time dims to a row vector to enable implicit replication when multiplying spatial dims by time dims
            if numel(sampTimes_a) > 1
                sampTimes_a = sampTimes_a(:).';
            end
            phs_spha_a = phs_spha_a(:,:);
            phs_conc_a = phs_conc_a(:,:);
			% Precompute spatial variation at all times (high memory demand, but very fast)
            phs = 0;
            for n=sphaInds
				% Add all spatially varying spherical harmonic terms
                bfunc = basisFuncHarm(x,y,z,n,subIndsSpace);
                bfunc = bfunc(:);
                if n>4 && ~isempty(obj.b0mask)
                    bfunc = bfunc.*b0mask_a;
                end
                phs_a = bfunc .* phs_spha_a(n,:);
				phs = phs + phs_a;
            end
            for n=1:size(phs_conc_a,1)
				% Add all concomitant gradient terms
                bfunc = basisFuncConc(x,y,z,n,subIndsSpace);
                bfunc = bfunc(:);
                if ~isempty(obj.b0mask)
                    bfunc = bfunc.*b0mask_a;
                end
                phs_a = bfunc .* phs_conc_a(n,:);
				phs = phs + phs_a;
            end
            if numel(sampTimes_a)>1
                if ~isempty(obj.b0mask)
                    b0_a = b0_a.*b0mask_a;
                end
                phs_a = b0_a.*sampTimes_a;
                phs = phs + phs_a;
            end
			kbase = phs;
		end

        function obj = setPhiDiv(obj,ksphaDiv,kconcDiv)
            if (obj.useSegmented ~= 2) || ~obj.useInterp
                % Check memory again, since PhiDiv greatly increases memory
                % demand.
                G = gpuDevice;
                avMem = G.AvailableMemory;
                numbytes = numel(obj.b0)*numel(obj.sampTimes)*8;
                if numbytes < 0.15*avMem
                    obj.useGPU = 1;
                    obj.useSingle = 0;
                    obj.useInterp = 0;
                elseif numbytes < 0.3*avMem
                    obj.useGPU = 1;
                    obj.useSingle = 1;
                    obj.useInterp = 0;
                else
                    obj.useGPU = 1;
                    obj.useSingle = 1;
                    obj.useInterp = 0;
                    obj.useSegmented = 1;
                    obj.useSegmentedNdiv = ceil(numbytes / (0.1*avMem));
                    obj.kbase = [];
                    obj.SegNpnts = ceil(prod(obj.imSize)/obj.useSegmentedNdiv);
                    obj.SegNpntsAdj = ceil(prod(obj.kSize)/obj.useSegmentedNdiv);
                    obj.SegPhs = prepForSegmented(obj);
                end
            end
                
            % ksphaDiv and kconcDiv are temporal derivatives of obj.phs_spha and obj.phs_conc
            if obj.useSingle
                ksphaDiv = single(ksphaDiv);
				kconcDiv = single(kconcDiv);
            end
            % Move variables to GPU
			if obj.useGPU
				ksphaDiv = gpuArray(ksphaDiv);
				kconcDiv = gpuArray(kconcDiv);
			end
            obj.ksphaDiv = ksphaDiv;
            obj.kconcDiv = kconcDiv;
        end
		
		function [svdTime,svdSpace,traj,phsShft] = prepForInterp(obj)
			% Create nufft object
			% TODO: below assumes perfectly axial slices. To do this properly, need to:
			% 1. have normal vector to slice as an optional input
			% 2. find linear combination of terms 2:4 in kspha for in-plane to slice
			% 3. set that to kloc, and substract from kspha. Then can still have all kspha terms in sum below
            if (obj.NDim ~= 2)
                error('Interpolation currently only supported for 2D in object domain')
            end
            sz_spha = size(obj.phs_spha);
            kloc = zeros([2,sz_spha(2:end)], 'like', obj.phs_spha);
            if abs(obj.phs_grid.z(2,2)-obj.phs_grid.z(1,1)) ~= 0
                error('non-axial slices not supported for interpolated approach')
            end
            if abs(obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1)) > eps
                if abs(obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1)) ~= 0
                    error('non-axial slices not supported for interpolated approach')
                end
                res1 = abs(obj.phs_grid.y(2,1)-obj.phs_grid.y(1,1));
                res2 = abs(obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1));
                if obj.phs_grid.y(2,1)-obj.phs_grid.y(1,1) > 0
                    kloc(1,:) = obj.phs_spha(3,:)/2/pi*res1;
                else
                    kloc(1,:) = -obj.phs_spha(3,:)/2/pi*res1;
                end
                if obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1) > 0
                    kloc(2,:) = obj.phs_spha(2,:)/2/pi*res2;
                else
                    kloc(2,:) = -obj.phs_spha(2,:)/2/pi*res2;
                end
                % Create phase ramp to center object domain properly, since
                % nufft assumes isocenter is at matrix center
                phsShft = obj.phs_grid.y(end/2+1,1)*obj.phs_spha(3,:) + ...
                    obj.phs_grid.x(1,end/2+1)*obj.phs_spha(2,:);
            else
                res1 = abs(obj.phs_grid.x(2,1)-obj.phs_grid.x(1,1));
                res2 = abs(obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1));
                if obj.phs_grid.x(2,1)-obj.phs_grid.x(1,1) > 0
                    kloc(1,:) = obj.phs_spha(2,:)/2/pi*res1;
                else
                    kloc(1,:) = -obj.phs_spha(2,:)/2/pi*res1;
                end
                if obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1) > 0
                    kloc(2,:) = obj.phs_spha(3,:)/2/pi*res2;
                else
                    kloc(2,:) = -obj.phs_spha(3,:)/2/pi*res2;
                end
                % Create phase ramp to center object domain properly, since
                % nufft assumes isocenter is at matrix center
                phsShft = obj.phs_grid.x(end/2+1,1)*obj.phs_spha(2,:) + ...
                    obj.phs_grid.y(1,end/2+1)*obj.phs_spha(3,:);
            end
            phsShft = reshape(exp(1i*phsShft), size(obj.sampTimes));
			traj = nufftOp(size(obj.b0), kloc(:,:)',[],obj.useGPU);
			clear kloc
			% Determine full non-linear encoding matrix, subsampling along
			% time since phase is slowly varying in time. 
            if obj.subFact>1
                inds = 1:obj.subFact:numel(obj.sampTimes);
                if inds(end) ~= numel(obj.sampTimes)
                    % Keep the ends to avoid extrapolation
                    inds = [inds, numel(obj.sampTimes)]';
                end
            end
            % We also interp in space via obj.subFactSpc, but this should
            % be only 1 or 2.
            indsSpc = [];
            if obj.subFactSpc>1
                sz = [obj.imSize,1,1];
                indsSpc = false(sz);
                if sz(3) > 1
                    error('this part untested for 3D')
                end
                for n2 = 1:sz(2)
                    for n3 = 1:sz(3)
                        strt = mod(n3+mod(n2,obj.subFactSpc), obj.subFactSpc);
                        if strt==0
                            strt = obj.subFactSpc;
                        end
                        indsSpc(strt:obj.subFactSpc:end,n2,n3) = true;
                    end
                end
                % Keep the ends to avoid extrapolation
                indsSpc(1,:,:) = true;
                indsSpc(end,:,:) = true;
                if sz(2)>1
                    indsSpc(:,1,:) = true;
                    indsSpc(:,end,:) = true;
                end
                if sz(3)>1
                    indsSpc(:,:,1) = true;
                    indsSpc(:,:,end) = true;
                end
            end
			b = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes,[1,4:size(obj.phs_spha,1)],[],[],indsSpc,inds);
            b = exp(1i*b);
            % Find largest singular values and vectors
            S = 1;
            ntry = 0;
            if obj.svdThresh < 0.055
                delTry = 30;
            else
                delTry = 20;
            end
            subspcFact = 3;
            while (min(diag(S))/max(S(:)) > obj.svdThresh) 
                if ntry > 200
                    error('many large singular values. Try providing mask or adjusting svdThresh.')
                end
                ntry = ntry + delTry;
                if obj.useSingle
                    b = double(b);
                end
                [U,S,V,FLAG] = svds(b,ntry,'largest','SubspaceDimension',subspcFact*ntry);
                if obj.useSingle
                    U = single(U);
                    S = single(S);
                    V = single(V);
                end
                if FLAG
                    warning('svd failure to converge. Increasing subspace.')
                    ntry = ntry - delTry;
                    subspcFact = subspcFact+1;
                end
            end
            Ns = find(diag(S)/max(S(:))<obj.svdThresh,1,'first');
            svdTime = conj(V(:,1:Ns)*S(1:Ns,1:Ns));
            if obj.subFact > 1
                % Fill back in values if interpolation was used
                svdTime = interp1(inds,gather(svdTime),1:inds(end),'pchip');
            end
            svdSpace = U(:,1:Ns);
            if obj.subFactSpc > 1
                % Fill back in values if interpolation was used in space
                % This is actually pretty slow, because griddata does not
                % support GPU.
                if length(obj.imSize) == 1
                    warning('1D untested')
                    indsSpc = find(indsSpc);
                    svdSpace = interp1(indsSpc,gather(svdSpace),1:indsSpc(end),'pchip');
                elseif length(obj.imSize) == 2
                    [X,Y] = ndgrid(1:obj.imSize(1),1:obj.imSize(2));
                    svdSpaceOut = zeros(prod(obj.imSize),size(svdSpace,2));
                    for n=1:size(svdSpace,2)
                        tmp = griddata(X(indsSpc),Y(indsSpc),gather(svdSpace(:,n)),X(:),Y(:),'natural');
                        svdSpaceOut(:,n) = tmp(:);
                    end
                    svdSpace = svdSpaceOut;
                    clear X Y svdSpaceOut
                elseif length(obj.imSize) == 3
                    warning('3D untested')
                    [X,Y,Z] = ndgrid(1:obj.imSize(1),1:obj.imSize(2),1:obj.imSize(3));
                    svdSpaceOut = zeros(prod(obj.imSize),size(svdSpace,2));
                    for n=1:size(svdSpace,2)
                        tmp = griddata(X(indsSpc),Y(indsSpc),Z(indsSpc),gather(svdSpace(:,n)),X(:),Y(:),Z(:),'natural');
                        svdSpaceOut(:,n) = tmp(:);
                    end
                    svdSpace = svdSpaceOut;
                    clear X Y svdSpaceOut
                end
            end
            clear U V
			% Compute error wrt direct approach
			if (0)
				erVal = gather(svdSpace)*gather(svdTime.');
				b = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes,[1,4:size(obj.phs_spha,1)]);
				b = reshape(b, numel(obj.b0), numel(obj.sampTimes));
                b = exp(1i*b);
                erVal = erVal(:) - gather(b(:));
                erVal = norm(erVal)/norm(b(:))
			end
        end

        function  res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
        function out = subArray(obj, x, n, out, doAdd)
            if nargin < 4
                out = [];
            end
            if nargin < 5
                doAdd = 0;
            end
            if obj.adjoint
                NdimIn = obj.NDimk;
                NdimOut = obj.NDim;
            else
                NdimIn = obj.NDim;
                NdimOut = obj.NDimk;
            end
            if (nargin < 4) || isempty(out)
                switch NdimIn
                    case 1
                        out = x(:,n);
                    case 2
                        out = x(:,:,n);
                    case 3
                        out = x(:,:,:,n);
                    otherwise
                        error('NdimIn not allowed');
                end
                if obj.useGPU && ~isa(out, 'gpuArray')
                    out = gpuArray(out);
                end
            else
                if ~isa(out,'gpuArray')
                    out_a = gather(x);
                end
                switch NdimOut
                    case 1
                        if doAdd
                            out(:,n) = out(:,n) + out_a;
                        else
                            out(:,n) = out_a;
                        end
                    case 2
                        if doAdd
                            out(:,:,n) = out(:,:,n) + out_a;
                        else
                            out(:,:,n) = out_a;
                        end
                    case 3
                        if doAdd
                            out(:,:,:,n) = out(:,:,:,n) + out_a;
                        else
                            out(:,:,:,n) = out_a;
                        end
                    otherwise
                        error('NdimOut not allowed');
                end
            end
        end
        
	end
end


