function [y,yW]=sincKernel(origGr,destGr,gpu)

%SINCKERNEL   Builds a sinc kernel used for interpolation of unstructred grids
%   [Y,YW]=SINCKERNEL(ORIGGR,DESTGR,{GPU})
%   * ORIGGR is the grid where the data is sampled (No in points x No dims)
%   * DESTGR is the grid where to obtain the data (No out points x No dims)
%   * {GPU} indicates whether to generate a gpu array for the kernel
%   * Y is the sinc kernel of size No out points x No in points
%   * YW is a density estimation of size No out points x 1
%

if ~exist('gpu','var') || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
           
[M,N]=parUnaFun({destGr,origGr},@size,1);
NDims=size(origGr,2);
assert(size(destGr,2)==NDims,'The dimensions of the destination grid (%d) should be the same as those of the origin grid (%d)',size(destGr,2),NDims);
y=single(ones([N M]));
if gpu;y=gpuArray(y);origGr=gpuArray(origGr);destGr=gpuArray(destGr);end

destGr=permute(destGr,[2 1]);
for m=1:NDims    
    K=sinc(bsxfun(@minus,origGr(:,m),destGr(m,:)));
    y=y.*K; 
end
y=permute(y,[2 1]);
yW=sum(y,2);
yW=abs(yW);