classdef nufftOp
% Perform nufft using method of Beatty et al (doi.org/10.1109/TMI.2005.848376)
%
% Usage: Define nufft object using N = nufftOp(imN, kloc, {dcf, useGPU, os,
% kwidth}), where the inputs in {} are optional. Then, use N*im to apply
% type 2 nufft (i.e., Cartesian to non-Cartesian) and use N'*kvals to apply
% type 1 nufft (i.e., non-Cartesian to Cartesian). im can have higher
% number of dims than length(imN) or kvals have more than 1 dim, in which
% case the nufft is applied separately to each entry of the higher dims.
%
% Convolution matrices are completely precomputed, which gives large speed
% advantages for iterative approaches, but slightly higher memory
% requirements. Matrix operations are used when possible, for a highly
% efficient implementation in MATLAB.
%
% Inputs:
%   imN [1 x NumDims]: number of voxels in each dimension of cartesian
%       domain (typically image domain). length(imN) is the number of
%       dimensions nuffts are performed over.
%   kloc [N x NumDims]: non-Cartesian k-space locations, scaled to be
%       between -0.5 and 0.5. Values outside [-0.5,0.5) will be circularly
%       shifted into [-0.5,0.5), akin to the DFT.
%   dcf [N x 1], default = []: density compensation to be used for
%       non-iterative regridding. Speeds up convergence at the expense of SNR
%       for iterative solutions. See Baron et al, https://doi.org/10.1002/mrm.26928
%   useGPU, default = false: whether to use GPU
%   os, default = 1.5: oversampling ratio. See Beatty et al (doi.org/10.1109/TMI.2005.848376)
%   kwidth, default = 4: kernel width. See Beatty et al (doi.org/10.1109/TMI.2005.848376)
%
% Member Functions:
%   N = N.findDcf(numIterations,isIt): finds an optimal dcf. 
%       numIterations, default = 10: number of iterations. Usually
%           converges in less than 10 iterations anyway.
%       isIt, default = false: set to true to use a near optimal dcf for
%           iterative solutions as specified in Baron et al, https://doi.org/10.1002/mrm.26928
%           Note that as implemented, W from the paper is equal to diag(obj.dcf.^2).
%           e.g., to solve the problem argmin sqrt(W)(Ax-y) with lsqr, where y is acquired k-space data, use: 
%               N = nufftOp(...); 
%               N = N.findDcf([],true);
%               x0 = N'*y;
%               opt.imNFull = size(x0); % required for vectorized outputs needed for lsqr
%               Nfun = @(x,transp) mrSampFunc(x,transp,N,[],opt); % mrSampFunc creates the net MRI sampling operation
%               rootWy = N.dcf.*y(:);
%               x = lsqr(Nfun,rootWy,tol,maxIt,[],[],x0);
%               x = reshape(x, N.imNFull);
%
%   N = N.prepToep(os,kwidth): preps for Toeplitz usage, which causes N*x
%       to actually give N'*(N*x). Uses Toeplitz properties to avoid
%       convolutions. See Baron et al, https://doi.org/10.1002/mrm.26928
%       - os and kwidth are optional inputs. Otherwise based on values in N
%       object that is used to call prepToep
%       - interestingly, with a precomputed convolution matrix, as used
%       here, there seems to be no speed benefit to using Toeplitz with
%       kwidth<=6 and os<=1.5. However, memory usage will be slightly lower.
%         - with Toeplitz, you can increase kwidth and os for no speed
%         penalty during iterations (prepToep will be slower, though)
%
% (c) Corey Baron 2020

% TODO: should make a "low memory" option where the entire convolution
% matrix is not precomputed. Instead, it could loop through the number of
% kernel shifts (i.e., loop through numel(sampOffsets{nD})

	properties
		kwidth = 6;     % kernel width. Must be even
		os = 1.5;       % oversampling factor
		useGPU = 0;	    % whether to use gpu. 
        useSingle = 0;  % whether to use single precision to save memory. NB: matlab currently does not support single sparse arrays, so this is not yet possible...
        loopDim = 4;    % maximum number of dims to do matrix based mtimes rather than loops. Trade-off between speed and memory requirements
        dcf = [];		% density compensation
    end

	properties (SetAccess = protected)
		adjoint = 0;
		kloc = []; % [N dim] k-space trajectory, scaled to [-0.5, 0.5] for all dims 
		ksize = [];	      % size of k-matrix that was inputted
		dcfsize = [];
		dcfMask = [];
        osN = [];
        comp = [];
        convMat = [];       % sparse matrix for convolution operation. If toeplitz requested, holds transfer fcn
        imgN = [];   		% image matrix size
        isToep = false;
        nofftShift = 1; % do not explicitely use fftshifts to save time. 
        fftShiftArray = [];
    end

	methods
		function obj = nufftOp(varargin) %nufftOp(imgN_in, k_in, dcf_in, useGPU, os, kwidth)
			if nargin < 1
				% Perform tests if no inputs
				obj.tests;
				return;
			end
			% Inputs: imgN_in, res_in, k_in, dcf_in, useGPU, os, kwidth
			noptions = length(varargin);
			obj.imgN = varargin{1};
			obj.kloc = varargin{2};
            if (noptions > 2) && ~isempty(varargin{3})
                obj.dcf = varargin{3};
            end
			if (noptions > 3) && ~isempty(varargin{4})
				obj.useGPU = varargin{4};
			end
			if (noptions > 4) && ~isempty(varargin{5})
				obj.os = varargin{5};
			end
            if (noptions > 5) && ~isempty(varargin{6})
				obj.kwidth = varargin{6};
            end
            % Check if gpu is possible
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected or GPU not supported. Using CPU.')
                obj.useGPU = 0;
            end
            try gpuDevice;
            catch
                obj.useGPU = 0;
                warning('CUDA version likely not supported. Using CPU.')
            end
            % Convert data types
            if obj.useSingle
                obj.kloc = single(obj.kloc);
                if ~isempty(obj.dcf)
                    obj.dcf = single(obj.dcf);
                end
            elseif ~isa(obj.kloc,'double')
                obj.kloc = double(obj.kloc);
                if ~isempty(obj.dcf)
                    obj.dcf = double(obj.dcf);
                end
            end
            if obj.useGPU
                obj.kloc = gpuArray(obj.kloc);
                if ~isempty(obj.dcf)
                    obj.dcf = gpuArray(obj.dcf);
                end
            end
			% Do precalculations
            obj = obj.doPrep;
        end

        function obj = doPrep(obj)
            if obj.nofftShift  
                % Perform fftshift in k-space by simply adjusting k-space indices
                kloc_a = obj.kloc + 0.5;
            else
                kloc_a = obj.kloc;
            end
			% Set oversampled grid size 
			obj.osN = ceil(obj.imgN*obj.os/2)*2; % Make a factor of 2
			% Minimize number of prime factors
            for nD=1:length(obj.osN)
				while (max(factor(obj.osN(nD))) > 7)
					obj.osN(nD) = obj.osN(nD) + 2;
				end
            end
			% Precompute gridding kernel. 
            precision = 10^-6;
			bet = pi*sqrt(obj.kwidth^2/obj.os^2*(obj.os-0.5)^2-0.8);
			Nkern = ceil(0.5*(obj.kwidth+1)*sqrt(0.37/obj.os^2/precision))*2;
			kx = linspace(-1,1,Nkern)';
            if obj.useSingle
                kx = single(kx);
            end
            kern = cell(length(obj.imgN), 1);
            for nD=1:length(obj.imgN)
                % Note: extra zero added for edge case where kernInds = Nkern (kernFrac should be exactly 0 in this case)
                kern{nD} = [besseli(0,bet*sqrt(1-kx.^2)); 0]; 
                kern{nD} = kern{nD}/sum(kern{nD})*(length(kern{nD})-2)/obj.kwidth; % normalize
                if obj.useGPU
                    kern{nD} = gpuArray(kern{nD});
                end
            end
			% Precompute compensation function
            for nD=1:length(obj.imgN)
				x = (-obj.imgN(nD)/2:obj.imgN(nD)/2-1)*pi*obj.kwidth/obj.osN(nD);
                if obj.useSingle
                    x = single(x);
                end
				x = sqrt(x.^2-bet^2);
				obj.comp{nD} = x./sin(x);
                obj.comp{nD} = obj.comp{nD}/min(obj.comp{nD}); % Normalize
                % Account for different fft scalings for different os (so
                % that choice of os does not scale result)
                obj.comp{nD} = obj.comp{nD}*sqrt(obj.osN(nD)/obj.imgN(nD));
                obj.comp{nD} = reshape(obj.comp{nD}(:), [ones(1,nD-1), length(obj.comp{nD}), 1]);
                if obj.useGPU
                    obj.comp{nD} = gpuArray(obj.comp{nD});
                end
            end
            % Precompute kernel weighting terms using linear interpolation
            C = ones(size(kloc_a,1),obj.kwidth^length(obj.imgN), 'like', kloc_a);
            sampOffsets = cell(length(obj.imgN),1);
            switch length(obj.imgN)
                case 1
                    sampOffsets{1} = (-obj.kwidth/2+1:obj.kwidth/2)';
                case 2
                    [sampOffsets{1}, sampOffsets{2}] =...
                        ndgrid((-obj.kwidth/2+1:obj.kwidth/2)');
                case 3
                    [sampOffsets{1}, sampOffsets{2}, sampOffsets{3}] =...
                        ndgrid((-obj.kwidth/2+1:obj.kwidth/2)');
            end
            if obj.useSingle
                obj.osN = single(obj.osN);
                Nkern = single(Nkern);
                obj.kwidth = single(obj.kwidth);
                for nD = 1:length(obj.imgN)
                    sampOffsets{nD} = single(sampOffsets{nD});
                end
            end
            if obj.useGPU
                obj.osN = gpuArray(obj.osN);
                Nkern = gpuArray(Nkern);
                obj.kwidth = gpuArray(obj.kwidth);
                for nD = 1:length(obj.imgN)
                    sampOffsets{nD} = gpuArray(sampOffsets{nD});
                end
            end
            for nD = 1:length(obj.imgN)
                kernInds_a = kloc_a(:,nD)*obj.osN(nD);
                kernInds_a = floor(kernInds_a)-kernInds_a;
                for nS = 1:numel(sampOffsets{1})
                    kernInds = kernInds_a+sampOffsets{nD}(nS); % in units of Cartesian k-space grid samples
                    kernInds = kernInds*(Nkern-1)/obj.kwidth+Nkern/2+0.5; % Convert to units of kern samples
                    kernFrac = kernInds-floor(kernInds);
                    kernInds = floor(kernInds);
                    % Linearly interpolate
                    C(:,nS) = C(:,nS) .* ( (1-kernFrac).*kern{nD}(kernInds) + kernFrac.*kern{nD}(kernInds+1) );
                end
            end
            clear kernInds_a kernInds kernFrac 
            % Determine Cartesian indices that correspond to each kernel term
            cartInds = cell(length(obj.imgN), 1);
            for nD=1:length(obj.imgN)
                cartInds{nD} = floor(kloc_a(:,nD)*obj.osN(nD));
                cartInds{nD} = repmat(sampOffsets{nD}(:)', size(kloc_a,1),1) +...
                    cartInds{nD} + obj.osN(nD)/2 + 1;
                % Wrap indices around (i.e., account for circular shift property of fft)
                cartInds{nD}(cartInds{nD} > obj.osN(nD)) = cartInds{nD}(cartInds{nD} > obj.osN(nD)) - obj.osN(nD);
                cartInds{nD}(cartInds{nD} < 1) = cartInds{nD}(cartInds{nD} < 1) + obj.osN(nD);
            end
            switch length(obj.imgN)
                case 1
                    cartInds = cartInds{1};
                case 2
                    cartInds = sub2ind(obj.osN,cartInds{1},cartInds{2}); 
                case 3
                    cartInds = sub2ind(obj.osN,cartInds{1},cartInds{2},cartInds{3});
            end
            % Create sparse matrix for convolution operation (Cartesian to
            % non-Cartesian)
            nonCartInds = (1:size(kloc_a,1))';
            if obj.useSingle
                nonCartInds = single(nonCartInds);
            end
            if obj.useGPU
                nonCartInds = gpuArray(nonCartInds);
            end
            nonCartInds = repmat(nonCartInds, [1 obj.kwidth^length(obj.imgN)]);
            obj.convMat = sparse(nonCartInds(:),cartInds(:),C(:),size(kloc_a,1),prod(obj.osN)); 
            if obj.nofftShift  
                % Perform image domain fftshift using phase ramp in k-space
                % We also account for non-symmetric zero padding used by
                %   fft and ifft fcns
                fftShiftMat = -((0:obj.osN(1)-1)' - obj.osN(1)/2) * obj.imgN(1)/obj.osN(1);
                if length(obj.imgN)>1
                    fftShiftMat = repmat(fftShiftMat, [1 obj.osN(2)]) -...
                        ((0:obj.osN(2)-1) - obj.osN(2)/2) * obj.imgN(2)/obj.osN(2);
                end
                if length(obj.imgN)>2
                    fftShiftMat = repmat(fftShiftMat, [1 1 obj.osN(3)]) -...
                        reshape((0:obj.osN(3)-1) - obj.osN(3)/2, [1 1 obj.osN(3)]) * obj.imgN(3)/obj.osN(3);
                end
                fftShiftMat = exp(1i*pi*fftShiftMat);
                obj.fftShiftArray = fftShiftMat;
            end
		end

		function obj = ctranspose(obj)
			obj.adjoint = xor(obj.adjoint,1);
		end

		function res = transpose(a)
			res = a';
        end
        
        function y = mtimes(obj,x)
            szx = [size(x),1,1,1,1];

            % Data type conversion. 
            inClass = [];
            if ~obj.useSingle 
                if ~isa(x,'gpuArray') && ~isa(x,'double')
                    inClass = ones(1,'like',x);
                    x = double(x);
                elseif isa(x,'gpuArray') && ~isaUnderlying(x,'double')
                    inClass = ones(1,'like',x);
                    tmp = ones(1,'gpuArray');
                    x = cast(x,'like',tmp); 
                    clear tmp
                end
            end

            % Prepare for repetitions
            if obj.isToep
                Nextra = szx(length(obj.imgN)+1:end);
                y = zeros([obj.imgN, Nextra], 'like', x);
                loopDim_pos = max(obj.loopDim,length(obj.imgN));
                loopDim_pre = loopDim_pos;
                NextraInloop = szx(length(obj.imgN)+1:loopDim_pre);
            else
                if obj.adjoint 
                    % Non-cartesian samples should always be a column vector
                    Nextra = szx(2:end);
                    y = zeros([obj.imgN, Nextra], 'like', x);
                    % Determine looping dims
                    loopDim_pos = max(obj.loopDim,length(obj.imgN));
                    loopDim_pre = loopDim_pos-length(obj.imgN)+1;
                    NextraInloop = szx(2:loopDim_pre);
                else
                    Nextra = szx(length(obj.imgN)+1:end);
                    y = zeros([size(obj.kloc,1), Nextra], 'like', x);
                    % Determine looping dims
                    loopDim_pre = max(obj.loopDim,length(obj.imgN));
                    loopDim_pos = loopDim_pre-length(obj.imgN)+1;
                    NextraInloop = szx(length(obj.imgN)+1:loopDim_pre);
                end
            end
            Nrep = prod(szx(loopDim_pre+1:end));
            if isempty(NextraInloop)
                NextraInloop = 1;
            end

            for nR = 1:Nrep
                % Extract repetition
                y_a = getSub(x,loopDim_pre,nR);
                wasNotGpu = false;
                if obj.useGPU && ~isa(x, 'gpuArray')
                    try
                        y_a = gpuArray(y_a);
                    catch ME
                        if contains(ME.message, 'memory') || strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
                            msg = 'Set nufftOp option loopDim to a smaller value if out of memory (e.g., S = nufftOp(...); S.loopDim = 2;)';
                            causeException = MException('MATLAB:sampHighOrder:memory',msg);
                            ME = addCause(ME,causeException);
                        end
                        rethrow(ME);
                    end
                    wasNotGpu = true;
                end
                if obj.isToep
                    % Toep adjoint and forward are equivalent
                    % TODO: this could probably be done with avoided explicit fftshifts
                    y_a = padcrop(fftnc(obj.convMat.*ifftnc(padcrop(y_a,2*obj.imgN),...
                        length(obj.imgN)),length(obj.imgN)),obj.imgN);
                else
                    if obj.adjoint
                        % non-Cartesian to Cartesian
                        if ~isempty(obj.dcf)
                            y_a = y_a.*obj.dcf;
                        end
                        y_a = obj.convMat'*y_a(:,:); 
                        y_a = reshape(y_a, [obj.osN,NextraInloop]);
                        if ~isempty(obj.fftShiftArray)
                            y_a = y_a.*conj(obj.fftShiftArray);
                        end
                        y_a = fftnc(y_a,length(obj.imgN),~obj.nofftShift);
                        sz_a = [size(y_a),1,1,1,1];
                        y_a = y_a(1:obj.imgN(1),:,:,:);
                        if length(obj.imgN)>1
                            y_a = y_a(:,1:obj.imgN(2),:,:);
                        end
                        if length(obj.imgN)>2
                            y_a = y_a(:,:,1:obj.imgN(3),:);
                        end
                        y_a = reshape(y_a, [obj.imgN, sz_a(length(obj.imgN)+1:end)]);
                        for nD=1:length(obj.imgN)
                            % Apply compensation for kernel transfer fcn
                            y_a = y_a.*obj.comp{nD};
                        end
                    else
                        % Cartesian to non-Cartesian
                        for nD=1:length(obj.imgN)
                            % Precompensation is separable
                            y_a = y_a.*obj.comp{nD};
                        end
                        y_a = ifftnc(y_a,length(obj.imgN),~obj.nofftShift,1,obj.osN);
                        if ~isempty(obj.fftShiftArray)
                            y_a = y_a.*obj.fftShiftArray;
                        end
                        y_a = reshape(y_a,prod(obj.osN),[]);
                        y_a = obj.convMat*y_a(:,:);
                        y_a = reshape(y_a,[size(obj.kloc,1),NextraInloop]);
                        if ~isempty(obj.dcf)
                            y_a = y_a.*obj.dcf;
                        end
                    end
                end
                if isa(y_a, 'gpuArray') && wasNotGpu
                    y_a = gather(y_a);
                end
                % Insert result into full array
                y = setSub(y,y_a,loopDim_pos,nR); 
            end

            % Revert class
            if ~isempty(inClass)
                y = cast(y,'like',inClass);
            end
            
        end % end mtimes
        
        function obj = prepToep(obj,tos,kwidth)
            if ~obj.isToep
                % Finding toeplitz only requires one nufft, so accuracy
                % requirement lower than nufft for S'*(S*x). Thus, reduce
                % os to save memory
                if nargin<2 || isempty(tos)
                    tos = 0.5*(obj.os-1)+1; 
                end
                if nargin<3 || isempty(kwidth)
                    kwidth = gather(obj.kwidth);
                end
                S = nufftOp(2*obj.imgN, gather(obj.kloc),gather(obj.dcf),obj.useGPU,tos,kwidth);
                if ~isempty(obj.dcf)
                    obj.convMat = S'*obj.dcf;
                else
                    obj.convMat = S'*ones(size(obj.kloc,1),1,'like',obj.kloc);
                end
                obj.convMat = (2^length(obj.imgN)) * ifftnc(obj.convMat);
                obj.isToep = true;
                % Toeplitz requires less memory for precomputations vs nufft, so allow for bigger arrays to save time
                obj.loopDim = obj.loopDim + 1; 
            end
        end
        
        function obj = findDcf(obj,nit,isIt)
            % Based on the least-squares algorithm in Bydder et al. MRM, 25(5):695-702, 2007
            % Also inspired by	N. Dwork et al. A Least Squares Optimal Density Compensation Function for Gridding. ISMRM Annual Meeting, Honolulu, April 2017 
            if nargin<3
                isIt = false;
            end
            if nargin<2 || isempty(nit)
                nit = 10;
            end
            
            % Find mask that defines the region of support for obj.convMat
            valAfter = prod(obj.imgN./obj.osN); % Note that obj.convMat'*ones = prod(obj.imgN./obj.osN) for perfect Cartesian sampling at Nyquist
            thresh = 0.1 * valAfter; 
            mask = obj.convMat'*ones([size(obj.kloc,1),1],'like',obj.kloc) > thresh;
            mask = double(mask);
            if obj.useGPU
                mask = gather(mask);
            end
            
            % Solve: argmin_dcf( ||mask.*(obj.convMat'*dcf = valAfter)||_2^2 )
            %        s.t. dcf >= 0 
            C = gather(obj.convMat');
            mask_a = sparse(1:numel(mask),1:numel(mask),mask(:));
            C = mask_a*C;
            d = valAfter*mask(:).*ones([prod(obj.osN),1],'like',obj.kloc); 
            lb = zeros([size(obj.kloc,1),1],'like',obj.kloc);
            if obj.useGPU
                %warning('nufftOp.findDcf: As of MATLAB 2020a, lsqlin does not support gpu input. Switching to CPU.')
                %C = gather(C);
                d = gather(d);
                lb = gather(lb);
            end
            options = optimoptions('lsqlin','MaxIterations',nit,'display','none');
            obj.dcf = lsqlin(C,d,[],[],[],[],lb,[],[],options);
            if obj.useGPU
                obj.dcf = gpuArray(obj.dcf);
            end
            
            if isIt
                % See https://doi.org/10.1002/mrm.26928
                  % Note that as implemented, W from the paper is equal to diag(obj.dcf.^2)
                obj.dcf = (obj.dcf).^0.25; 
            end
        end
        
	end
end

function out = getSub(in,ND,n)
	switch ND
	case 1
		out = in(:,n);
	case 2
		out = in(:,:,n);
	case 3
        out = in(:,:,:,n);
    case 4
        out = in(:,:,:,:,n);
    case 5
        out = in(:,:,:,:,:,n);
    case 6
        out = in(:,:,:,:,:,:,n);
    otherwise
        error('loopDim too high')
	end
end

function out = setSub(out,in,ND,n)
	switch ND
	case 1
		out(:,n) = in;
	case 2
		out(:,:,n) = in;
	case 3
        out(:,:,:,n) = in;
    case 4
        out(:,:,:,:,n) = in;
    case 5
        out(:,:,:,:,:,n) = in;
    case 6
        out(:,:,:,:,:,:,n) = in;
    otherwise
        error('loopDim too high')
	end
end