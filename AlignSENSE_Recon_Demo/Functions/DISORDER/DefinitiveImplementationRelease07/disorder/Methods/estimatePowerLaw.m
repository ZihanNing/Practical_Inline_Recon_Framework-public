function [W,dV,xn,Wn,p]=estimatePowerLaw(x,useFinite,thUp,invW)

% ESTIMATEPOWERLAW estimates a logarithm fit to the spectral power assuming
% a log-normal distribution of the observations, which lie in a 3D space
%   [W,DV]=ESTIMATEPOWERLAW(X,{USEFINITE})
%   * X is a residuals or spectral power density array
%   * {USEFINITE} serves to use finite differences, defaults to 0
%   * {THUP} serves to use an upper threshold in the frequencies used to 
%   compute the law
%   * {INVW} serves to invert the weights
%   ** W is the estimated law for the power
%   ** DV are the exponents of the laws
%   ** XN are the residuals (for visualization)
%   ** WN are the spectral radii (for visualization)
%   ** P are the terms of the fit (for visualization)
%

if nargin<2 || isempty(useFinite);useFinite=0;end
if nargin<3;thUp=[];end
if nargin<4 || isempty(invW);invW=0;end

ND=numDims(x);
N=size(x);N(end+1:3)=1;N3=N(1:3);
gpu=isa(x,'gpuArray');

[x,NN]=resSub(x,4:ND);NN(end+1:4)=1;
if ~useFinite 
    W=WFiso;
    Wor=W;
elseif useFinite==1
    W=1./buildFilter(N3,'FractionalFiniteDiscreteIso',ones(1,3),gpu,-1);
    Wor=W;
elseif useFinite==2    
    Wor=WFiso;%We apply in a non-flat manner
    W=1./buildFilter(N3,'FractionalFiniteDiscreteIso',ones(1,3),gpu,-1);%We fit to flat spectrum for high frequencies
end

dV=[];
for n=1:NN(4)
    Wn=W;
    xn=dynInd(x,n,4);
    Wn(xn==0)=[];
    xn(xn==0)=[];
    if ~isempty(thUp)
        xn(Wn>=thUp)=[];
        Wn(Wn>=thUp)=[];
    end
    
    Wn=log2(Wn(:));
    xn=log2(xn(:));
    NP=length(xn);
    A=cat(2,ones([NP,1],'like',Wn),Wn);
    p=single(gather(double(A))\gather(double(xn)));
    dV=cat(4,dV,p(2));
    Wn=[];xn=[];A=[];
end;x=[];
if gpu;dV=gpuArray(dV);end
if invW;dV=-dV;end
W=bsxfun(@power,Wor,dV);
W=reshape(W,N);
