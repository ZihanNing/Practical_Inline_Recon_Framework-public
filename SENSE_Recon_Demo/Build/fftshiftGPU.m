function x=fftshiftGPU(x,m,n)

%FFTSHIFTGPU   Performs an fftshift using a for loop
%   X=FFTSHIFTGPU(X,M,{N})
%   * X is an array
%   * M is the dimension along which to perform the shift
%   * {N} is the minimum dimension along which to loop, it defaults to 4
%   ** X is the shifted result
%

if nargin<3 || isempty(n);n=4;end

NX=size(x);NX(end+1:n)=1;
ND=length(NX);
[x,N]=resSub(x,n:ND);N(end+1:n)=1;
for l=1:N(n);x=dynInd(x,l,n,fftshift(dynInd(x,l,n),m));end
x=reshape(x,NX);
