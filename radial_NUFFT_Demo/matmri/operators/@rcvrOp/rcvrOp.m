classdef rcvrOp
% Calibrate receivers from low res data, and apply either forward or
% adjoint receiver operation
%
% Usage: Define rcvr object using R = rcvrOp(input,isCalDat,{apN, kalph}).
% Inputs in {} are optional.
% Use R*im to apply receiver sensitivities to an image and use R'*rcvrIm to
% apply adjoint operation (conjugate sum over receivers).
%
% Inputs:
%   input: either object-domain receiver sensitivity data (when isCalDat=0) or
%       image domain calibration data (when isCalDat = 1). When calibration data is
%       supplied, a direct method is used to find the profiles
%       (McKenzie et al, https://doi.org/10.1002/mrm.10087). Apodization is
%       performed using a kaiser window. Last non-zero dim is assumed to be
%       receivers dimension.
%   isCalDat: see description for "input"
%   apN: dimensions for apodization edge (i.e., size of calibration
%       region). Default = size(input). 
%   kalph: alpha for kaiser window used for apodization. Default = 4.
%
% Parameters that can be set:
%   R.doLoop (default = false): set to true to loop over receivers to
%       reduce memory requirements
%   R.useGPU (default = true): set to false to use CPU. Will default to CPU
%       if no GPU is detected.
%
% (c) Corey Baron 2020



	properties
        kalph = 4;      % alpha parameter for apodization
        apN = [];       % size of calibration region
        maps = [];      % receiver maps
		useGPU = 1;	    % whether to use gpu. 
        doLoop = 0;     % Use loops over receivers to save memory
    end

	properties (SetAccess = protected)
		adjoint = 0;
        imgN = [];      % Size of image (including receiver dimension)
        nDim = [];      % number of dimensions (excluding receiver dims)
    end

	methods
		function obj = rcvrOp(varargin) %rcvrOp(input,isCalDat,{apN, kalph})
            if nargin < 1
				% Perform tests if no inputs
				obj.tests;
				return;
            end
			input = varargin{1};
			isCalDat = varargin{2};
            if length(varargin)>3 && ~isempty(varargin{4})
                obj.kalph = varargin{4};
            end
            obj.imgN = size(input);
            % Check if gpu is possible
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected or GPU not supported. Using CPU.')
                obj.useGPU = 0;
            elseif ~isa(input,'gpuArray')
                input = gpuArray(input);
            end
            if isCalDat
                obj.nDim = length(obj.imgN)-1; 
                if length(varargin)>2 && ~isempty(varargin{3})
                    obj.apN = varargin{3};
                else
                    obj.apN = obj.imgN(1:obj.nDim);
                end
                if length(obj.apN) > obj.nDim
                    error('apN must have length <= to number of dimensions')
                else
                    fftDim = length(obj.apN);
                end
                % Apodize
                win = filtNd(obj.apN,obj.kalph,length(obj.apN),'kb');
                if isa(input,'gpuArray')
                    win = gpuArray(win);
                end
                win = padcrop(win,obj.imgN(1:fftDim));
                obj.maps = ifftnc(input,fftDim);
                obj.maps = obj.maps.*win;
                obj.maps = fftnc(obj.maps,fftDim);
                % Compute maps
                obj.maps = obj.maps./sqrt(sum(obj.maps.*conj(obj.maps),...
                    length(obj.imgN)));
            else
                obj.nDim = length(obj.imgN)-1; % NB: this will need to be adjusted for ESPIRiT-style maps with more than one set of maps
                obj.maps = input;
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
            
            % Check gpu status
            inputGpu = isa(x,'gpuArray');
            
            % Prepare for repetitions
            Nrep = 1;
            if obj.doLoop
                Nrep = obj.imgN(obj.nDim+1);
            end
            
            if obj.doLoop
                if obj.adjoint
                    y = 0;
                    for n=1:Nrep
                        x_a = getSub(x,obj.nDim,n);
                        if ~inputGpu && obj.useGPU
                            x_a = gpuArray(x_a);
                        end
                        map_a = getSub(obj.maps,obj.nDim,n);
                        y = y + x_a.*conj(map_a);
                    end
                    if ~inputGpu && obj.useGPU
                        y = gather(y);
                    end
                else
                    if ~inputGpu && obj.useGPU
                        x = gpuArray(x);
                    end
                    y = zeros(size(obj.maps),'like',x);
                    for n=1:Nrep
                        map_a = getSub(obj.maps,obj.nDim,n);
                        y_a = x.*map_a;
                        if ~inputGpu && obj.useGPU
                            y_a = gather(y_a);
                        end
                        y = setSub(y,y_a,obj.nDim,n); 
                    end
                end
            else
                if ~inputGpu && obj.useGPU
                    x = gpuArray(x);
                end
                if obj.adjoint
                    y = x.*conj(obj.maps);
                    y = sum(y,obj.nDim+1);
                else
                    y = x.*obj.maps;
                end
                 if ~inputGpu && obj.useGPU
                    y = gather(y);
                end
            end
            
        end % end mtimes
        
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
	end
end