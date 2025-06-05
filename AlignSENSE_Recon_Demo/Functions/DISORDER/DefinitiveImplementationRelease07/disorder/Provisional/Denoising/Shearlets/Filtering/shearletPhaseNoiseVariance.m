function [x2,n2]=shearletPhaseNoiseVariance(x,n2,ND,sH)

%SHEARLETPHASENOISEVARIANCE computes the shearlet phase noise variance
%   X2=SHEARLETNOISEVARIANCE(X,N2,{ND},{SH}) computes the empirical shearlet
%   variance on a per scale basis and uses noise variance estimates to 
%   obtain an estimate for the signal variance
%   X is the observed signal
%   N2 is the noise variance
%   ND are the number of dimensions, it defaults to 3
%   {SH} is a shearlet structure
%   X2 is the estimated signal variance
%

if nargin<3 || isempty(ND);ND=[];end
if nargin<4;sH=[];end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:3)=1;
if ~isempty(sH);x=shearletTransform(x,sH);end%Otherwise we assume it is already transformed
NDD=numDims(x);
x=cat(NDD+1,real(x),imag(x));
x2=x.^2;
spacing=1/8;
strength=1;
H=buildFilter(N(1:ND),'tukeyIso',spacing,gpu,strength);
%x2=multDimSum(x2,1:ND)/prod(N(1:ND));

x2=sum(x2,NDD+1);
x2=max(real(filtering(x2,H)),0);
n2=max(real(filtering(n2,H)),0);
x2=max(bsxfun(@minus,x2,n2),0);
x2=sqrt(x2);
n2=sqrt(n2);

