function [x,s]=augmentScale(x,factor)

%AUGMENTSCALE   Augments a given set of training images by applying a 
%random scaling
%   [X,S]=AUGMENTSCALE(X,FACTOR)
%   * X is the original data
%   * FACTOR is the scaling factor
%   * X is the augmented data
%   * S is the set of applied scalings
%

if nargin<2 || isempty(factor);factor=1.5;end

N=size(x);N(end+1:5)=1;
NSamp=N(5:end);NSampTotal=prod(NSamp);
NChan=N(4);
N=N(1:3);

s=factor.^(rand([prod(NSamp) 1])-0.5);
NFull=round(bsxfun(@times,N,s));
x=reshape(x,[N NChan NSampTotal]);
for l=1:NSampTotal;x(:,:,:,:,l)=resampling(resampling(x(:,:,:,:,l),NFull(l,:)),N,2);end
s=reshape(s,[ones(1,4) NSamp]);
x=reshape(x,[N NChan NSamp]);
