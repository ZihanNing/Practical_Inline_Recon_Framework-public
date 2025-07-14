function a = ifftnc(a, N, doShift, doScale, sz)
% Calculates the multidimensional ifft of a matrix with DC at the
% center of the matrix. 
%
% Inputs:
%   a: input matrix
%   N: do fft's along first "N" dim
% 
% Outputs:
%   A: output matrix
% 
% (c) Corey Baron

if (nargin<5)
    sz = [];
end
if(nargin<4) || isempty(doScale)
  doScale = 1;
end
if(nargin<3) || isempty(doShift)
  doShift = 1;
end
if(nargin<2) || isempty(N)
  N = length(size(a));
end

if (N > length(size(a)))
  N = length(size(a));
end

scale = 1;
if doShift
  for n=1:N   
    a = ifftshift(a,n);
    if ~isempty(sz)
        a = ifft(a,sz(n),n);
    else
        a = ifft(a,[],n);
    end
    a = ifftshift(a,n);
    scale = scale*sqrt(size(a,n));
  end
else
  for n=1:N
    if ~isempty(sz)
        a = ifft(a,sz(n),n);
    else
        a = ifft(a,[],n);
    end
    scale = scale*sqrt(size(a,n));
  end
end

if doScale
  a = a*scale;
end	
