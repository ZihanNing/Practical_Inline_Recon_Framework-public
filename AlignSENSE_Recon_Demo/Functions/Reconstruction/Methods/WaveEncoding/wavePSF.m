
function [PSF, PSF_PE1, PSF_PE2] = wavePSF (N, RO, numHelix, gMax)

if length(gMax)==1; gMax = gMax * ones([1 2]);end
  
perm = 1:3;
perm = [ RO setdiff(perm, RO) ];
N  = N(perm);

[RO,PE1,PE2] = ndgrid(linspace(-1,1,N(1)) ,linspace(-1,1,N(2)),linspace(-1,1,N(3)));
[kRO,kPE1,kPE2] = ndgrid(linspace(-pi,pi,N(1)) ,linspace(-pi,pi,N(2)),linspace(-pi,pi,N(3)));

T = 2*pi/numHelix;
deltaK = kRO(2) - kRO(1);%linear

gPE1 = gMax(1)*cos(2*pi/T*kRO);
gPE2 = gMax(2)*sin(2*pi/T*kRO);%Phase lag
Cy = -2*pi*1i*cumsum(gPE1,1)*deltaK;
Cz = -2*pi*1i*cumsum(gPE2,1)*deltaK;

PSF_PE1 = exp(Cy.*PE1);
PSF_PE2 = exp(Cz.*PE2);

PSF = PSF_PE1 .* PSF_PE2;

%%% Permute back
PSF = ipermute(PSF, perm);

end