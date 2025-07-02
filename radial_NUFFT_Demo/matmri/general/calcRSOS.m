function m_sos = calcRSOS(m, dim)
% Calculate root sum-of-squares over specified dimension 
%
% y = calcRSOS(x,dim)
% 
% (c) Corey Baron

if (nargin<2)
  dim = ndims(m);
end

m_sos = squeeze(sqrt(sum(m.*conj(m), dim)));
