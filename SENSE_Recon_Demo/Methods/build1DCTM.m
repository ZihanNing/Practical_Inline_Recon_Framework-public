function [F,FH]=build1DCTM(N,fn,gpu)

%BUILD1CFTM   Builds a standard DCT matrix
%   [F,FH]=BUILD1DFTM(N,{FN},{GPU})
%   * N is the dimension of the space
%   * {FN} indicates whether to generate fully unitary Fourier matrices. It
%   defaults to 0 (not used!!)
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   ** F is a discrete cosine transform matrix
%   ** FH is a an inverse discrete cosine transform matrix
%

if nargin<2 || isempty(fn);fn=0;end
if nargin<3 || isempty(gpu);gpu=useGPU;end


F=single(dctmtx(N));
if gpu;F=gpuArray(F);end
if nargout>1;FH=F';end
