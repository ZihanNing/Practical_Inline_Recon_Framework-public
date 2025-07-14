function [R, theta, phi] = sphereGrid(sz)
% Like meshgrid, but for spherical coordinates
%
% (c) Corey Baron
  
theta = [];
phi = [];
  
X = (-sz(1)/2:sz(1)/2-1)/(sz(1)/2);
if length(sz)>1
  Y = (-sz(2)/2:sz(2)/2-1)/(sz(2)/2);
end
if length(sz)>2
  Z = (-sz(3)/2:sz(3)/2-1)/(sz(3)/2);  
end

switch length(sz)
  case 1
    R = abs(X);
  case 2
    [X,Y] = ndgrid(X,Y);
    Z = 0;
  case 3
    [X,Y,Z] = ndgrid(X,Y,Z);
end

if length(sz)>1
  R = sqrt(X.^2+Y.^2+Z.^2); 

  if nargout > 1
    theta = acos(Z./R);
  end

  if nargout > 2
    phi = atan2(Y,X);
  end
end
