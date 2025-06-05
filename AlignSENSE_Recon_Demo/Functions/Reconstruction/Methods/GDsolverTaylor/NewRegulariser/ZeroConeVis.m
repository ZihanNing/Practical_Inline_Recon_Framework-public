
clc
clear all
close all

N = 80*ones(1,3);
[Kx,Ky,Kz] = linearTerms(N, pi);

KMagn = sqrt((Kx.^2 + Ky.^2 + Kz.^2));

theta = 0;%degrees
kernel0 = (1/3  - (-Kx.*sind(theta) + Kz .* cosd(theta)).^2 ./(KMagn.^2)  );

theta = 30;%degrees
kernel1 = (1/3  - (-Kx.*sind(theta) + Kz .* cosd(theta)).^2 ./(KMagn.^2)  );

theta = 50;%degrees
kernel2 = (1/3  - (-Kx.*sind(theta) + Kz .* cosd(theta)).^2 ./(KMagn.^2)  );

threshold=0.03;
%plot_(kernel0, kernel1,kernel2,[],[],[],0,0)
%plot_(kernel0, abs(kernel1)<threshold,abs(kernel2)<threshold,[],[],[],0,0)

%% Kernel for susc. model estimation for linear model
% theta = 1;%degrees
% kernelM = (Kx*theta).*(Kz./(KMagn.^2));
%plot_(kernelM, kernelM,abs(kernelM)<threshold,[],[],[],0,0)

%% Checl assumption linearity
theta = 8;
L1 = (2*Kx .* Kz)     ./(KMagn.^2)    *sind(theta)*cosd(theta);   L1(KMagn==0) = 0;
L2 = (Kz.^2 - Kx.^2)  ./(KMagn.^2)    *sind(theta)^2;             L2(KMagn==0) = 0;

l1 = 10^6* real( fftshift(ifftn(ifftshift(L1))));
l2 = 10^6* real( fftshift(ifftn(ifftshift(L2))));
lcomb = l1+l2;

plot_(lcomb,[],[],l1,l2,[-40 40],[],1,[],[],224)

X = zeros(N);
X(KMagn<pi/6) = 1;
kernel = 1/3 - Kz.^2./(KMagn.^2); kernel(KMagn==0)=0;%X(KMagn==0);
kernel = kernel + L1 + L2;
B = real(ifftn( ifftshift(kernel) .* fftn(X)));
plot_(X,[],[],B,B,[-0.3 0.3],[],1,[],[],225)
