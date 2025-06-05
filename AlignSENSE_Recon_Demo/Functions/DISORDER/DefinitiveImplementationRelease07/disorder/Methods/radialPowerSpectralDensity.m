function [f,r]=radialPowerSpectralDensity(x,voxsiz,ou,fo)

% RADIALPOWERSPECTRALDENSITY computes the power spectral density of X at 
% different wavelenghts 
%   [F,R]=RADIALPOWERSPECTRALDENSITY(X,VOXSIZ,OU) 
%   * X is the data
%   * {VOXSIZ} is the voxel size, it defaults to all 1
%   * {OU} is the output type, 'li' for standard scale, 'dB' for decibels
%   * {FO} indicates the space, if 0, it is in image space, if 1 it is in
%   Fourier space, non shifted, but in power units, if 2 it is in Fourier
%   space, shifted and in power units
%   ** F is the power
%   ** R is the radious
%

if nargin<3 || isempty(ou);ou='li';end
if nargin<4 || isempty(fo);fo=0;end

gpu=isa(x,'gpuArray');
NF=size(x);NF(end+1:4)=1;N=NF(1:3);
kGrid=generateGrid(N,gpu,2./voxsiz,ceil((N+1)/2));
rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,kGrid{1}.^2,kGrid{2}.^2),kGrid{3}.^2));
if N(3)==1;N=N(1:2);end
ND=length(N);
spR=max(2./(voxsiz(1:ND).*N));%Spectral resolution, the factor 2 is arbitrary
rGrid=rGrid/spR;
rGrid=floor(rGrid)+1;
for n=1:3
    if fo==0;x=fftGPU(x,n);end
    if fo~=2;x=fftshift(x,n);end
end
if fo==0;x=abs(x).^2;end
x=x(:,:,:,:);
y=x(:,:,:,1);
y(:)=1;w=accumarray(rGrid(:),y(:));
NS=size(x,4);
f=zeros([length(w) NS],'like',y);
for n=1:size(x,4)
    y=x(:,:,:,n);
    f(:,n)=accumarray(rGrid(:),y(:));
end
f=bsxfun(@rdivide,f,w+eps);
r=(1:max(rGrid(:)))';
r=(r-1)*spR;
f=reshape(f,[],1,1,NF(4:end));
if strcmp(ou,'dB');f=10*log10(f);end
