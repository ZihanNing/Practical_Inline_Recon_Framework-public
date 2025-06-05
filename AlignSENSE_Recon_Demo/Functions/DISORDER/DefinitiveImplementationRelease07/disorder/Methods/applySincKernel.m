function x=applySincKernel(x,K,KW,factReg)

%APPLYSINCKERNEL   Applies a sinc kernel for interpolation / regridding
%   X=APPLYSINCKERNEL(X,K,KW,{FACTREG})
%   * X is the data in original locations
%   * K is the sinc kernel of size No out points x No in points
%   * KW is a density estimation of size No out points x 1
%   * {FACTREG} is a regularization factor for areas where the density is smaller
%   ** X is the data in destination locations
%

if ~exist('factReg','var');factReg=0.5;end

[N,M]=parUnaFun({K,x},@size);
M(end+1:2)=1;
C=prod(M)/N(2);
assert(mod(C)==0,'Data samples (%d) cannot be placed in original locations (%d)',prod(M),N(2));
x=reshape(x,[N(2) prod(C)]);
x=K*x;

KWM=mean(KW)*factReg;
KW=max(KW,KWM);
x=bsxfun(@times,x,1./KW);
x=reshape(x,[N(1) M(2:end)]);
