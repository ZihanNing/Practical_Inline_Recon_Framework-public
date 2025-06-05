function [x,n,A]=standardizeCoils(x,n,de)

%STANDARDIZECOILS   Applies standardization to the receiver channels 
%   X=STANDARDIZECOILS(X,N,{DE})
%   * X is an array
%   * N is an array of noise samples or a covariance matrix (considering
%   both quadrature receivers). Its meaning is determined from its
%   dimensionality
%   * {DE} is a flag to reverse the operation and destandardize the channels.
%   It defaults to 0 (standardization)
%   * X is the standardized set of measures (the method assumes equal dwell
%   times for noise and sampling and flat frequency response of noise samples)
%   * N as an output is the covariance matrix
%

if nargin<3 || isempty(de);de=0;end

gpu=isa(x,'gpuArray') || isa(n,'gpuArray');
if gpu;x=gpuArray(x);end

[N,M]=parUnaFun({n,x},@size);
M(end+1:5)=1;
if ~(length(N)==2 && N(1)==N(2) && N(1)==M(4))%Then n is assumed not to be a covariance matrix
    n=mapMat(n,4);
    n=cov(n);
end

[V,D]=eig(n);
if ~de;A=(V*diag(1./(diag(D/2).^(1/2)))*V');else A=(V*diag(diag(D/2).^(1/2))*V');end%Factor 2 comes from complex data (this way instead of unit power spectral density, we have unit variance in real and imaginary components)

N=size(x);N(end+1:4)=1;
x=resSub(x,5:14);
N5=size(x,5);
for l=1:N5
    xn=x(:,:,:,:,l);
    xn=reshape(xn,[prod(N(1:3)) N(4)]);
    xn=xn*A;
    xn=reshape(xn,N(1:4));
    x(:,:,:,:,l)=xn;
end
x=reshape(x,N);

