% Class of object that represents discrete wavelet transform. This performs an
% object domain implementation for haar wavelet, and frequency domain
% implementation for larger filter lengths.
%
% Corey Baron, 2016

classdef dwt
properties
  adjoint = 0;
  inverse = 0; % Set to 1 to compute (W^H W)^-1 W^H. Note that (W^H W)^-1 = I for normal decimated or undecimated wavelet, but not 1 depending on "A" input.
  doFourier = false;
  J;
  imSize;
  imSizeIn;
  family;
  decimation; % true or false
  APrep = [];     % Any precomputed data required for AFcn and/or AHAinvFcn
  AHAinvFcn = []; % Function for computing AHAinv. Inputs: (x, obj.APrep)
  AFcn = [];      % Function for computing Ax, when using an extra fcn appended on wavelet (i.e., DWTnet = [DWT; A]). Inputs: (x, obj.APrep, doAdjoint)
  Lo_D = [];
  Hi_D = [];
  Lo_R = [];
  Hi_R = [];
  useGPU=1;    % Does NOT gather back from GPU after decomposition. This makes for time savings in iterative methods that have a recon shortly after a decomp
end

methods
  function  obj = dwt(J,imSize, decimation, family, useGPUFlag, AFcn, AHAinvFcn, APrep)    % constructor
      if nargin<1
          % Perform unit tests if no inputs
          obj.tests;
          return;
      end
      if nargin<6
        AFcn = [];
        AHAinvFcn = [];
        APrep = [];
      end
      if nargin<5
        useGPUFlag = 1;
      end
      if nargin<4
        family = 'db1';
      end
      if nargin<3
        decimation = false;
      end
      % Check if gpu is possible
      obj.useGPU = useGPUFlag;
      if ~(gpuDeviceCount>0) && obj.useGPU
        warning('No GPU detected or GPU not supported. Using CPU.')
        obj.useGPU = 0;
      end

      J(imSize==1) = 0; % Don't try to do convolution along size 1 dims
      obj.J = cat(1,J(:),zeros(20,1)); % Allow replicated inputs

      obj.decimation = decimation;
      obj = setFFTchoice(obj,family);

      obj.imSizeIn = imSize;
      J_a = obj.J(1:length(imSize))';
      if obj.decimation && any(mod(imSize./2.^(J_a-1),2))
        % For decimation it is not possible to maintain odd matrix sizes and have orthogonal transforms
        N = ceil(imSize./2.^(J_a-1));
        N(J==0) = 0;
        N = N + mod(N,2);
        imSize = N.*2.^(J_a-1);
      end
      obj.imSize = imSize;

      obj.family = family;
      if obj.useGPU
        obj = prepGPU(obj);
      end
      obj.AFcn = AFcn;
      obj.AHAinvFcn = AHAinvFcn;
      obj.APrep = APrep;
  end

  function obj = prepGPU(obj)
    if iscell(obj.Lo_D)
      for n1 = 1:size(obj.Lo_D,1)
        for n2 = 1:size(obj.Lo_D,2)
          obj.Lo_D{n1,n2} = gpuArray(single(obj.Lo_D{n1,n2}));
          obj.Hi_D{n1,n2} = gpuArray(single(obj.Hi_D{n1,n2}));
          obj.Lo_R{n1,n2} = gpuArray(single(obj.Lo_R{n1,n2}));
          obj.Hi_R{n1,n2} = gpuArray(single(obj.Hi_R{n1,n2}));
        end
      end
    else
      obj.Lo_D = gpuArray(single(obj.Lo_D));
      obj.Hi_D = gpuArray(single(obj.Hi_D));
      obj.Lo_R = gpuArray(single(obj.Lo_R));
      obj.Hi_R = gpuArray(single(obj.Hi_R));
    end
  end

  function obj = setFFTchoice(obj,family)
    tmp = wfilters(family);
    if length(tmp) > 2
      obj.doFourier = true;
    else
      obj.doFourier = false;
    end
  end

  function obj = set.family(obj,family)
    obj.family = family;
    [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(obj.family);
    if ~obj.decimation
      Lo_D = Lo_D/sqrt(2);
      Hi_D = Hi_D/sqrt(2);
      Lo_R = Lo_R/sqrt(2);
      Hi_R = Hi_R/sqrt(2);
    end

    if obj.doFourier
      obj.Lo_D = cell(max(obj.J), length(obj.imSize));
      obj.Hi_D = cell(max(obj.J), length(obj.imSize));
      obj.Lo_R = cell(max(obj.J), length(obj.imSize));
      obj.Hi_R = cell(max(obj.J), length(obj.imSize));
      for j = 1:max(obj.J)
        for d = 1:length(obj.imSize)
          if j <= obj.J(d)
            Lo_D_a = Lo_D;
            Hi_D_a = Hi_D;
            if obj.decimation
              N_a = ceil(obj.imSize(d)/2^(j-1));
            else
              N_a = obj.imSize(d);
              Lo_D_a = upsample(Lo_D_a, 2^(j-1));
              Hi_D_a = upsample(Hi_D_a, 2^(j-1));
            end
            Lo_D_a = circshift([Lo_D_a(:); zeros(N_a-length(Lo_D_a),1)], -length(Lo_D_a)/2);
            Hi_D_a = circshift([Hi_D_a(:); zeros(N_a-length(Hi_D_a),1)], -length(Hi_D_a)/2);
            obj.Lo_D{j, d} = fft(Lo_D_a);
            obj.Hi_D{j, d} = fft(Hi_D_a);
            obj.Lo_R{j, d} = conj(obj.Lo_D{j, d});
            obj.Hi_R{j, d} = conj(obj.Hi_D{j, d});
          end
        end
      end
    else
      obj.Lo_D = Lo_D;
      obj.Hi_D = Hi_D;
      obj.Lo_R = Lo_R;
      obj.Hi_R = Hi_R;
    end
  end % end set.family

  function x = dofft(obj,x,fftType)
    if nargin<3
      fftType = 'for';
    end

    gpuGather = 0;

    if ~iscell(x)
      if isa(obj.Lo_D{1}, 'gpuArray') && ~isa(x, 'gpuArray')
        x = gpuArray(single(x));
        gpuGather = 1;
      end
      if all(obj.J(1:ndims(x)))
        if strcmp(fftType,'inv')
          x = ifftn(x);
        else
          x = fftn(x);
        end
      else
        for n=1:ndims(x)
          if obj.J(n) > 0
            if strcmp(fftType,'inv')
              x = ifft(x,[],n);
            else
              x = fft(x,[],n);
            end
          end
        end
      end
      if 0 %gpuGather
        x = gather(x);
      end

    else
      if all(obj.J(1:ndims(x{1})))
        for nC = 1:length(x)
          if isa(obj.Lo_D{1}, 'gpuArray') && ~isa(x{nC}, 'gpuArray')
            x{nC} = gpuArray(single(x{nC}));
            gpuGather = 1;
          end
          if strcmp(fftType,'inv')
            x{nC} = ifftn(x{nC});
          else
            x{nC} = fftn(x{nC});
          end
          if 0 %gpuGather
            x{nC} = gather(x{nC});
          end
        end
      else
        for nC = 1:length(x)
          if isa(obj.Lo_D{1}, 'gpuArray') && ~isa(x{nC}, 'gpuArray')
            x{nC} = gpuArray(single(x{nC}));
            gpuGather = 1;
          end
          for n=1:ndims(x{1})
            if (obj.J(n) > 0)
              if strcmp(fftType,'inv')
                x{nC} = ifft(x{nC},[],n);
              else
                x{nC} = fft(x{nC},[],n);
              end
            end
          end
          if 0 %gpuGather
            x{nC} = gather(x{nC});
          end
        end
      end
    end
  end

end
end
