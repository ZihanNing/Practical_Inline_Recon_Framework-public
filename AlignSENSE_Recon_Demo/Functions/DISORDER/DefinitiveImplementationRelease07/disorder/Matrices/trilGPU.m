function X=trilGPU(X)

%TRILGPU   Generates a set of lower triangular matrices from X using GPU computations
%   X=TRILGPU(X)
%   * X are the input matrices
%   * X are the output matrices
%

gpu=isa(X,'gpuArray');

N=size(X);N(end+1:3)=1;
ind1=(1:N(1))';
ind2=1:N(2);
if gpu;[ind1,ind2]=parUnaFun({ind1,ind2},@gpuArray);end
X=bsxfun(@times,X,bsxfun(@ge,ind1,ind2));
