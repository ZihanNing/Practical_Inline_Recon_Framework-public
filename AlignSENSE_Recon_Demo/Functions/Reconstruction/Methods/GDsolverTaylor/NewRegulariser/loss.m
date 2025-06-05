

function [L] = loss(D,Dmeas, lambda, x)

N = size(D);

%%% create linear terms
Nx=N(1:3);
rangex = (1:Nx(1)) - ( floor(Nx(1)/2)); rangex = rangex / (Nx(1)/2);
rangey = (1:Nx(2)) - ( floor(Nx(2)/2)); rangey = rangey / (Nx(2)/2);
rangez = (1:Nx(3)) - ( floor(Nx(3)/2)); rangez = rangez / (Nx(3)/2);

rangex= pi * rangex;
rangey= pi * rangey;
rangez= pi * rangez;

[Kx, Ky, ~] = ndgrid(rangex, rangey, rangez);

%%% regulariser
temp = Ky .* fftn(fftshift(dynInd(D,1,6))) + Kx .* fftn(fftshift(dynInd(D,2,6)));
temp = sqrt(abs(x)) .* ifftshift(ifftn(temp));

Lreg = normm(temp);

%%% data consistency
Ldata = normm(sqrt(abs(x)).*(D - Dmeas));

%%% sum
L = lambda * Lreg + Ldata;
