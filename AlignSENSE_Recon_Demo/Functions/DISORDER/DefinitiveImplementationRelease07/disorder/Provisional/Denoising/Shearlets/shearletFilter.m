function x=shearletFilter(x,sH,nF,bhat,wi,ph)

%SHEARLETFILTER performs a BayesShrink shearlet filtering
%   X=SHEARLETFILTER(X,SH,{NF},{BHAT},{WI},{PH}) performs a BayesShrink 
%   shearlet filtering
%   X is the observed signal
%   SH is a shearlet structure
%   {NF} are the noise factors
%   {BHAT} is the shape of the generalized Gaussian. If empty it is
%   estimated
%   {WI} is the width of the kernel used for parameter estimation
%   {PH} is a flag to estimate the phase noise model
%   X is the filtered signal
%

if nargin<3 || isempty(x);nF=ones(size(x),'like',x);end%noise factors
if nargin<4;bhat=[];end
if nargin<5 || isempty(wi);wi=[16 32];end%To estimate parameters for Bayesian denoising (sigma/beta)
if nargin<6 || isempty(ph);ph=0;end%noise factors

N=size(x,4);ND=numDims(x);
perm=1:ND+1;perm([4 ND+1])=[ND+1 4];
if N~=1;x=permute(x,perm);end

co=~isreal(x);
x(nF<1e-6)=0;
n2=(1+co)*nF.^2;%Noise variance
if ph
    x2=abs(x).^2;%Signal variance
    n2=phaseVariance(x2,n2,2,ph-1);%Phase variance
    if ph==1%Has not worked very well although is free from artifacts, there seems to be problems with low frequencies?
        x=x./(abs(x)+eps);    
        n2=2*n2{2};
    else%Has not worked well at all
        x=angle(x);
        n2=n2{1};
    end
end
n2=bsxfun(@times,n2,sH.rms.^2);%Shearlet noise variance

x=shearletTransform(x,sH);
if ~isempty(bhat);[xstd,nstd]=shearletSignalVariance(x,n2,numDims(sH.dfw),[],wi);%Shearlet signal variance
else [xstd,nstd,bhat]=shearletSignalVariance(x,n2,numDims(sH.dfw),[],wi);
end

T=(1./sqrt(bhat)).*bsxfun(@times,nstd,bsxfun(@rdivide,nstd,xstd).^(sqrt(bhat)));
%T=((2-beta)./(2*(1-beta))).*(bsxfun(@rdivide,bsxfun(@times,2*(1-beta),nstd.^2),xstd.^beta).^(1./(2-beta)));%From Wavelet thresholding for some classes of non-Gaussian noise, note xstd is not the std but the "beta-moment" 
%T=sqrt(2)*(beta.^1.8).*bsxfun(@rdivide,(nstd.^2),xstd.^beta);%Not using a prior
%T=sqrt(n2)*sqrt(2*log....);%Universal
%T=x2./(bsxfun(@plus,x2,n2)+eps);%Wiener
x=shearletTransform(x,sH,0,T);
if ph==2;x=exp(1i*x);end

if N~=1;x=permute(x,perm);end
 