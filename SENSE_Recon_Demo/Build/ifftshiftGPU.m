function x=ifftshiftGPU(x,m,n)

%IFFTSHIFTGPU   Performs an ifftshift using a for loop
%   X=IFFTSHIFTGPU(X,M,{N})
%   * X is an array
%   * M is the dimension along which to perform the shift
%   * {N} is the minimum dimension along which to loop, it defaults to 4
%   ** X is the shifted result
%

if nargin<3 || isempty(n);n=4;end

NX=size(x);NX(end+1:n)=1;
ND=length(NX);
[x,N]=resSub(x,n:ND);N(end+1:n)=1;
for l=1:N(n);x=dynInd(x,l,n,ifftshift(dynInd(x,l,n),m));end
x=reshape(x,NX);
