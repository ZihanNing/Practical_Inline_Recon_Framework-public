function [x,no]=augmentNoise(x,si)

%AUGMENTNOISE   Augments a given set of training images by applying random
%Gaussian noise
%   [X,NO]=AUGMENTNOISE(X,SI)
%   * X is the original data
%   * SI is the noise level, defaults to 1 
%   * X is the augmented data
%   * NO is the generated noise
%

if nargin<2 || isempty(si);si=1;end

comp=~isreal(x);
no=plugNoise(x,1);
if comp;no=no/sqrt(2);end
no=bsxfun(@times,no,si);
x=bsxfun(@plus,x,no);
