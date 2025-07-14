function y = mtimes(obj,x)

  if ~xor(obj.adjoint, obj.inverse)   %(0,0), (1,1)
    if obj.useGPU
      x = gpuArray(x);
    end
    reg2 = [];
    if ~isempty(obj.AFcn)
      if obj.inverse
        x = obj.AHAinvFcn(x,obj.APrep);
      end
      reg2 = obj.AFcn(x,obj.APrep,0);
    end
    y = decomp(obj,x);
    y.extra = reg2;
    y = waveletObj(y);
  else
    y = recon(obj,x);
    if ~isempty(obj.AFcn)
      y = y + obj.AFcn(x.extra,obj.APrep,1);
      if obj.inverse
        y = obj.AHAinvFcn(y,obj.APrep);
      end
    end
    if obj.useGPU
      y = gather(y);
    end
  end

end

function x_wav = decomp(obj, x_obj)
  if any(obj.imSize ~= obj.imSizeIn)
    x_obj = padImgW(x_obj, obj.imSize - obj.imSizeIn);
  end

  if obj.doFourier
    x_obj = dofft(obj,x_obj);
  end

  x_wav.low = x_obj;
  for j = 1:max(obj.J)
    [x_wav.low, x_wav.high{j}] = decompSep(obj, x_wav.low, j);
    if obj.doFourier
      x_wav.high{j} = dofft(obj,x_wav.high{j},'inv');
    end
  end
  if obj.doFourier
    x_wav.low = dofft(obj,x_wav.low,'inv');
  end
end

function x_obj = recon(obj, x_wav)
  if obj.doFourier
    x_wav.low = dofft(obj,x_wav.low);
  end
  for j = max(obj.J):-1:1
    if obj.doFourier
      x_wav.high{j} = dofft(obj,x_wav.high{j});
    end
    x_wav.low = reconSep(obj, x_wav.low, x_wav.high{j}, j);
  end
  x_obj = x_wav.low;

  if obj.doFourier
    x_obj = dofft(obj,x_obj,'inv');
  end

  if any(obj.imSize ~= obj.imSizeIn)
    x_obj = cropImgW(x_obj, obj.imSizeIn);
  end
end

function [x_wav_low, x_wav_high] = decompSep(obj, x_wav_low, j)
  % The dimensions to be decomposed
  DecomposeDim = obj.J>=j;

  % Allocate DFTs of the wavelet coefficients
  x_wav_high = cell(1, 2^sum(DecomposeDim)-1);

  % Filter-bank analysis
  ds = 0;
  for d = 1:numel(DecomposeDim)
  	if DecomposeDim(d)
  		if nargout > 1
  			x_wav_high{2^ds} = wavConv1D(obj, x_wav_low,obj.Hi_D,j,d,0);
  			for s = 1:2^ds-1
  				x_wav_high{2^ds+s} = wavConv1D(obj, x_wav_high{s}, obj.Hi_D, j, d, 0);
  				x_wav_high{s} = wavConv1D(obj, x_wav_high{s}, obj.Lo_D, j, d, 0);
  			end
  		end
  		x_wav_low = wavConv1D(obj, x_wav_low, obj.Lo_D, j, d, 0);
  		ds = ds + 1;
  	end
  end
end

function x_wav_low = reconSep(obj, x_wav_low, x_wav_high, j)
  % The dimensions do be reconstructed
  ReconstructDim = obj.J>=j;

  % Filter-bank synthesis
  ds = sum(ReconstructDim) - 1;
  for d = numel(ReconstructDim):-1:1
  	if ReconstructDim(d)
  		if isempty(x_wav_low)
  			x_wav_low = 0;
  		else
  			x_wav_low = wavConv1D(obj, x_wav_low, obj.Lo_R, j, d, 1);
  		end
  		if ~isempty(x_wav_high)
  			x_wav_low = x_wav_low + wavConv1D(obj, x_wav_high{2^ds}, obj.Hi_R, j, d, 1);
  			x_wav_high{2^ds} = [];
  			for s = 1:2^ds-1
  				x_wav_high{s} = wavConv1D(obj, x_wav_high{s}, obj.Lo_R, j, d, 1) + wavConv1D(obj, x_wav_high{2^ds+s}, obj.Hi_R, j, d, 1);
  				x_wav_high{2^ds+s} = [];
  			end
  		end
  		ds = ds - 1;
  	end
  end
end

function x = cropImgW(x,sz)
  for n=1:length(sz)
    switch n
    case 1
      x = x(1:sz(n),:,:,:,:,:);
    case 2
      x = x(:,1:sz(n),:,:,:,:);
    case 3
      x = x(:,:,1:sz(n),:,:,:);
    case 4
      x = x(:,:,:,1:sz(n),:,:);
    case 5
      x = x(:,:,:,:,1:sz(n),:);
    case 6
      x = x(:,:,:,:,:,1:sz(n));
    end
  end
end

function x = padImgW(x,pad)
  for n=1:length(pad)
    if pad(n)
      sz = size(x);
      sz(n) = pad(n);
      x = cat(n,x,zeros(sz,'like',x));
      %switch n
      %case 1
      %  x = cat(n,x,x(1:pad(n),:,:,:,:,:));
      %case 2
      %  x = cat(n,x,x(:,1:pad(n),:,:,:,:));
      %case 3
      %  x = cat(n,x,x(:,:,1:pad(n),:,:,:));
      %case 4
      %  x = cat(n,x,x(:,:,:,1:pad(n),:,:));
      %case 5
      %  x = cat(n,x,x(:,:,:,:,1:pad(n),:));
      %case 6
      %  x = cat(n,x,x(:,:,:,:,:,1:pad(n)));
      %end
    end % end if
  end % end for
end
