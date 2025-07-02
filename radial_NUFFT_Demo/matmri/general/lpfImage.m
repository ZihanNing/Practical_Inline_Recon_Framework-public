function x = lpfImage(x,type,alpha,dim)
  % win = filtNd(sz,alpha,dim, type)
  % alpha: alpha for type='kb', gaussian fwhm (in samples) for 'gauss'
  % No filtering along dims of k>dim
  %
  % (c) Corey Baron

  if nargin < 2
    type = 'kb';
  end
  if nargin < 3
    switch type
    case 'kb'
      alpha = 3;
    case 'gauss'
      alpha = 1/2;
    end
  end
  if nargin < 4
    dim = ndims(x);
  end

  win = filtNd(size(x),alpha,dim,type);

  x = ifftnc(win.*fftnc(x,dim),dim);
end
