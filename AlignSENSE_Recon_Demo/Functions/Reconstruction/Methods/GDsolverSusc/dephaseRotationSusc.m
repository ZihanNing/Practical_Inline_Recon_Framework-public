
function [x,p]=dephaseRotationSusc(T,fun,chi,kernelStruct, TE,returnField)

%DEPHASESUSC   Computes the dephasing
%   X=DEPHASESUSC(T,CHI,TE, RETURNFIELD)

if nargin<4 || isempty(kernelStruct); error('dephaseSusc:: kernelStruct muse be provided to calcaulate susceptibility induced B0.');end
if nargin<5 || isempty(TE); TE=1;end
if nargin<6 || isempty(returnField); returnField=0;end

ndT=numDims(T);
N = size(chi);

%%% FOURIER DOMAIN
for i=1:3; chi=fftGPU(chi, i);end

%%% CREATE DIFFERENTIAL DIPOLE KERNELS
rotPar = fun(T);%Rotation parameters used in the susc model

kernel = cat(ndT, kernelStruct.Kroll , kernelStruct.Kpitch); 
kernel = sum(bsxfun(@times,kernel,rotPar),ndT);
%kernel(isnan(kernel))=0;

%%% FILTER
kernel = bsxfun(@times, fftshift(kernelStruct.H), kernel);%Filter kernel to avoid gibs artefacts
for i=1:3; kernel = ifftshiftGPU(kernel,i);end %chi not shifted so shift back

%%% APPLY KERNEL
chi = kernel.*chi;

%%% IMAGE DOMAIN
for i=1:3; chi=ifftGPU(chi, i);end%YB: attention!! This can be a IFT on 100 volumes so incredebly slow!!!!!!

H0 = 7; %Tesla
gamma= 42.58; %MHz/T
x=2* H0 * gamma * chi; %chi is in ppm so the 10^6 of gamma has been included here
if returnField; return; end %In Hz
    
if nargout>1; p = 2*pi*TE * x;end %phase in radians
x=exp(+1i *2*pi*TE * x);% YB: 2piTE since x in Hz and now in radians




