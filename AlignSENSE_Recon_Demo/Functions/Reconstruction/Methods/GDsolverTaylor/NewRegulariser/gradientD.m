

function [grad] = gradientD(D, Dmeas, lambda,x)

N = size(D);

grad = zeros(N);

%%% create linear terms
Nx=N(1:3);
rangex = (1:Nx(1)) - ( floor(Nx(1)/2)); rangex = rangex / (Nx(1)/2);
rangey = (1:Nx(2)) - ( floor(Nx(2)/2)); rangey = rangey / (Nx(2)/2);
rangez = (1:Nx(3)) - ( floor(Nx(3)/2)); rangez = rangez / (Nx(3)/2);

rangex= pi * rangex;
rangey= pi * rangey;
rangez= pi * rangez;

[Kx, Ky, ~] = ndgrid(rangex, rangey, rangez);

K = cat(6, Kx, Ky);

%%% regualriser
part2 = dynInd(K,2,6) .* fftn(fftshift(dynInd(D,1,6))) + dynInd(K,1,6) .* fftn(fftshift(dynInd(D,2,6)));
part2 = fftn(fftshift(abs(x).*ifftshift(ifftn(part2))));

for i =1:N(6)
    temp = ifftshift( ifftn( conj(dynInd(K,3-i,6)).* part2)) ;
    temp = real(temp);
    grad = dynInd(grad, i, 6, temp);
end

%%% data consistency
for i =1:N(6)
    grad = dynInd(grad, i, 6, dynInd(lambda * grad +  abs(x).*(D-Dmeas), i,6) );
end


