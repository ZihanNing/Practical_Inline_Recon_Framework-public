function y = wavConv1D(obj,x,filt,j,dim,dir)
  % dir = direction. 0 = decomposition, 1 = reconstruct.

  % Perform convolution
  if obj.doFourier
    if dim~=1
      filt_a = permute(filt{j,dim}(:), [2:dim, 1]);
    else
      filt_a = filt{j,dim}(:);
    end
    sz = size(x);
    sz(dim) = 1;
    filt_a = repmat(filt_a, sz);
    if obj.decimation && dir
      x = decRecF(x,dim,dir);
    end
    y = x.*filt_a;
    if obj.decimation && ~dir
      y = decRecF(y,dim,dir);
    end

  else

    if obj.decimation && dir
      x = decRecIm(x,dim,dir);
    end
    y = 0;
    for n = 1:length(filt)
      shift = n-1;
      if dir
        shift = shift - length(filt) + 1;
      end
      if ~obj.decimation
        shift = (2^(j-1))*shift;
      end
      if filt(n) ~= 0
        y = y + filt(n)*circshift(x,-shift,dim);
      end
    end
    if obj.decimation && ~dir
      y = decRecIm(y,dim,dir);
    end
  end
end

function x = decRecF(x,dim,dir)
  if ~dir
    phs = 0.5;
    switch dim
    case 1
      x = phs.*(x(1:end/2,:,:,:,:,:) + x(end/2+1:end,:,:,:,:,:));
    case 2
      x = phs.*(x(:,1:end/2,:,:,:,:) + x(:,end/2+1:end,:,:,:,:));
    case 3
      x = phs.*(x(:,:,1:end/2,:,:,:) + x(:,:,end/2+1:end,:,:,:));
    case 4
      x = phs.*(x(:,:,:,1:end/2,:,:) + x(:,:,:,end/2+1:end,:,:));
    case 5
      x = phs.*(x(:,:,:,:,1:end/2,:) + x(:,:,:,:,end/2+1:end,:));
    case 6
      x = phs.*(x(:,:,:,:,:,1:end/2) + x(:,:,:,:,:,end/2+1:end));
    end
  else
    x = cat(dim,x,x);
  end
end

function x = decRecIm(x,dim,dir)
  if ~dir
    switch dim
    case 1
      x = x(1:2:end,:,:,:,:,:);
    case 2
      x = x(:,1:2:end,:,:,:,:);
    case 3
      x = x(:,:,1:2:end,:,:,:);
    case 4
      x = x(:,:,:,1:2:end,:,:);
    case 5
      x = x(:,:,:,:,1:2:end,:);
    case 6
      x = x(:,:,:,:,:,1:2:end);
    end
  else
    sz_a = size(x);
    sz_a(dim) = 2*sz_a(dim);
    x_a = zeros(sz_a,'like',x);
    switch dim
    case 1
      x_a(1:2:end,:,:,:,:,:) = x;
    case 2
      x_a(:,1:2:end,:,:,:,:) = x;
    case 3
      x_a(:,:,1:2:end,:,:,:) = x;
    case 4
      x_a(:,:,:,1:2:end,:,:) = x;
    case 5
      x_a(:,:,:,:,1:2:end,:) = x;
    case 6
      x_a(:,:,:,:,:,1:2:end) = x;
    end
    x = x_a;
  end
end
