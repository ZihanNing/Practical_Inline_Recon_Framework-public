function win = filtNd(sz,alpha,dim,type)
  % win = filtNd(sz,alpha,dim, type)
  % alpha: alpha for type='kb', gaussian fwhm (wrt full width of kspace) for 'gauss'
  % No filtering along dims of k>dim
  %
  % (c) Corey Baron

  switch type
  case 'kb'
    filtFcn = @(szIn,alphaIn) kaiser(szIn,alphaIn);
  case 'gauss'
    filtFcn = @(szIn,alphaIn) normpdf(-szIn/2:szIn/2-1,0,szIn*alphaIn/2.355);
  otherwise
    error('Unknown filter type')
  end

  if (dim > 1) && length(alpha)==1
    alpha = ones(dim,1)*alpha;
  end

  win = filtFcn(sz(1),alpha(1));
  win = win(:);
  for n=2:dim
    k2 = permute(reshape(filtFcn(sz(n),alpha(n)), [], 1), [2:n, 1]);
    win = repmat(win, [ones(1,n-1), sz(n)]) .* repmat(k2, [sz(1:n-1), 1]);
  end
  if dim<length(sz)
    win = repmat(win, [ones(1,dim), sz(dim+1:end)]);
  end

  win = win/max(win(:));

end
