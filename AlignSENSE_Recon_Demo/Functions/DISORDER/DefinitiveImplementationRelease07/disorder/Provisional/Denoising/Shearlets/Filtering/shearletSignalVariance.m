function [x2s,n2s,bhat]=shearletSignalVariance(x,n2,ND,sH,wi)

%SHEARLETSIGNALVARIANCE computes the shearlet signal variance
%   X2=SHEARLETNOISEVARIANCE(X,N2,{ND},{SH},{WI}) computes the empirical shearlet
%   variance on a per scale basis and uses noise variance estimates to 
%   obtain an estimate for the signal variance
%   X is the observed signal
%   N2 is the noise variance
%   ND are the number of dimensions, it defaults to 3
%   {SH} is a shearlet structure
%   {WI} is the width of the kernel for estimation, it defaults to 32
%   X2S is the estimated signal variance
%   N2S is the estimated noise variance
%   BHAT is the estimated beta parameter
%

if nargin<3 || isempty(ND);ND=3;end
if nargin<4;sH=[];end
if nargin<5 || isempty(wi);wi=[16 32];end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:3)=1;
if ~isempty(sH);x=shearletTransform(x,sH);end%Otherwise we assume it is already transformed
NDD=numDims(x);
x=cat(NDD+1,real(x),imag(x));
x2=x.^2;

wi(end+1:2)=wi(1);
H=cell(1,2);
for n=1:2;H{n}=buildFilter(N(1:ND),'tukeyIso',1/wi(n),gpu,1);end%Kernel for estimation
x2=sum(x2,NDD+1);
x2s=max(real(filtering(x2,H{1})),0);
n2s=max(real(filtering(n2,H{1})),0);
if nargout>2
    x2b=max(real(filtering(x2,H{2})),0);
    x1b=max(real(filtering(sqrt(x2),H{2})),0);    
    bhat=estbeta(x1b,x2b);
end
x2s=max(bsxfun(@minus,x2s,n2s),0);
x2s=sqrt(x2s);n2s=sqrt(n2s);
