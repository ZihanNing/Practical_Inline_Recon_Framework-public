function [x,k]=estimatePSD(x,di,sp,ou)

% ESTIMATEPSD estimates the power spectral density of an image along a
% given dimension
%   [X,K]=ESTIMATEPSD(X,DI,{SP},{OU})
%   * X is the array on which to estimate the PSD
%   * DI is the dimension over which to estimate the PSD
%   * {SP} is the spacing along that dimension (in mm). It defaults to 1 
%   (mm)
%   * {OU} is the output type, 'li' for standard scale, 'dB' for decibels
%   ** X is the power spectral density estimate
%   ** K are the spatial frequenties (in 1/mm)
%

if nargin<3 || isempty(sp);sp=1;end
if nargin<4 || isempty(ou);ou='li';end


gpu=isa(x,'gpuArray');
ND=numDims(x);
x=fftGPU(x,di);
x=fftshift(x,di);
x=abs(x.*conj(x));
perm=1:max(ND,di);perm(1)=di;perm(di)=1;
x=permute(x,perm);
N=size(x);
x=reshape(x,[N(1) prod(N(2:end))]);
x=sp*mean(x,2)/N(1);
if strcmp(ou,'dB');x=10*log10(x);end

kaux=generateGrid(N(1),gpu,1,ceil((N(1)+1)/2));
k=kaux{1}/sp;
