function R = simRcvrSens(Nimg,NrcvrXY,NrcvrZ)
% Create receiver sensitivities with ~1/r^2 dependence and linear phase profile.
%   Some random fluctuations are introduced to make slightly more realistic.
%   Useful for reconstruction simulations.
%
% R = simRcvrSens(Nimg,NrcvrXY,NrcvrZ)
%
% Inputs:
%     Nimg: size of image. e.g., [32 32] for a 2D 32x32 matrix size.
%     NrcvrXY: number of receivers in the XY plane per Z-level. Distributed circularly in the plane.
%     NrcvrZ:  number of receivers along z-direction. If NrcvrZ = [] and length(Nimg)>2 (i.e., 3D), receivers are distributed in a sphere
%
% (c) Corey Baron 2017

% Check inputs
if length(Nimg) < 2 || length(Nimg) > 3
  error('Nimg must must be a length 2 (2D) or length 3 (3D) vector.')
end

% Distribute uniformly over a sphere, excluding bottom solid angle 'excl' (for neck hole).
% Otherwise, do a cylindrical placement of receivers.
if (nargin<3 || isempty(NrcvrZ)) && length(Nimg)>2
  NrcvrZ = 1;
  doSphere = 1;
  excl = pi/4;
else
  doSphere = 0;
end

% Define grid
if length(Nimg) == 2
  [X Y] = ndgrid(-Nimg(1)/2:Nimg(1)/2-1, -Nimg(2)/2:Nimg(2)/2-1);
  Z = 0;
  Nimg = [Nimg, 1];
  NrcvrZ = 1;
elseif length(Nimg) == 3
  [X Y Z] = ndgrid(-Nimg(1)/2:Nimg(1)/2-1, -Nimg(2)/2:Nimg(2)/2-1, -Nimg(3)/2:Nimg(3)/2-1);
else
  error('Need 2 or 3 img dim')
end

rng(1,'twister');
R = zeros([Nimg,NrcvrXY*NrcvrZ]);
if doSphere
  radFact = sqrt(3)*1.1;
  % Just spiral up the sphere (using phylotaxis) as a basic approximation
  n = 0:(NrcvrXY-1);
  rgold = 137.51/180*pi;
  phi = n*rgold;
  n0 = pi/2/(pi-excl)*(NrcvrXY-1);
  theta = pi/2*sqrt(n/n0);
  theta(n > n0) = pi - pi/2*sqrt((2*n0-n(n>n0))/n0);
  for n=1:NrcvrXY
    refpnt = mean(Nimg(1:2))*0.5*radFact*[sin(theta(n))*cos(phi(n)) sin(theta(n))*sin(phi(n)) cos(theta(n))];
    Rad2 = (X-refpnt(1)).^2+(Y-refpnt(2)).^2+(Z-refpnt(3)).^2;
    radmod = 0.85+rand*0.3;
    Phs = rand*X*5*pi/Nimg(1) + rand*Y*5*pi/Nimg(2) + rand*Z*5*pi/Nimg(3) + pi*rand;
    R(:,:,:,n) = (1./Rad2.^radmod).*exp(1i*Phs);
  end
else
  radFact = sqrt(2)*1.1;
  for nz = 1:NrcvrZ
    refpntz = nz*Nimg(3)/(NrcvrZ+1) - Nimg(3)/2;
    for n=1:NrcvrXY
      phi0 = pi/NrcvrXY;
      refpnt = mean(Nimg(1:2))*0.5*radFact*[sin(2*pi*n/NrcvrXY+phi0) cos(2*pi*n/NrcvrXY+phi0)];
      Rad2 = (X-refpnt(1)).^2+(Y-refpnt(2)).^2+(Z-refpntz).^2;
      radmod = 0.85+rand*0.3;
      Phs = rand*X*5*pi/Nimg(1) + rand*Y*5*pi/Nimg(2) + rand*Z*5*pi/Nimg(3) + pi*rand;
      R(:,:,:,n+(nz-1)*NrcvrXY) = (1./Rad2.^radmod).*exp(1i*Phs);
    end
  end
end

% Normalize
Rrsos = sqrt(abs(sum(R.*conj(R),4)));
for n=1:size(R,4)
  R(:,:,:,n) = R(:,:,:,n)./Rrsos;
end

R = squeeze(R);
