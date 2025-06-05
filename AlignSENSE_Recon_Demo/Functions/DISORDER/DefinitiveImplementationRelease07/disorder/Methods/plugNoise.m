function x=plugNoise(x,useR,useN)

% PLUGNOISE plugs complex circularly symmetric AWGN noise samples in an 
%array
%   X=PLUGNOISE(X,USER) 
%   * X is the input array
%   * USER enables the usage of real only noise
%   * USEN enables the usage of the normal distribution
%   ** X is the noise array with same dimensions as the input array
%

gpu=isa(x,'gpuArray');
if nargin<2 || isempty(useR);useR=0;end%Complex only, for back-compatibility!
if nargin<3 || isempty(useN);useN=1;end%Normal distribution

comp=~isreal(x) || ~useR;

N=size(x);
if useN
    if comp;x=single(randn(N))+1i*single(randn(N));else x=single(randn(N));end
else
    if comp;x=single(rand(N)-0.5)+1i*single(rand(N)-0.5);else x=single(rand(N)-0.5);end
end
if gpu;x=gpuArray(x);end
