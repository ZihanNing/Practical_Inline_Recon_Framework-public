function [F,FH]=buildStandardDCTM(N,fn,gpu)

%BUILDSTANDARDDCTM   Builds standard DCT matrices
%   [F,FH]=BUILDSTANDARDDCTM(N,{FN},{GPU})
%   * N are the dimensions of the space
%   * {FN} indicates whether to generate fully unitary Fourier matrices. It
%   defaults to 0 (not used!!)
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   ** F is a cell of discrete cosine transform matrices along the 
%   different dimensions
%   ** FH is a cell of inverse discrete cosine transform matrices along 
%   the different dimensions
%

if nargin<2 || isempty(fn);fn=0;end
if nargin<3 || isempty(gpu);gpu=useGPU;end

ND=length(N);
F=cell(1,ND);FH=cell(1,ND);
for m=1:ND
    if nargout>1;[F{m},FH{m}]=build1DCTM(N(m),fn,gpu);
    else F{m}=build1DCTM(N(m),fn,gpu);
    end
end
