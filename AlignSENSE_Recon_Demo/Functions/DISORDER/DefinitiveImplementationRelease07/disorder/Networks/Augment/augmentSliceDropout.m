function [x,d]=augmentSliceDropout(x,p,ab)

%AUGMENTSLICEDROPOUT   Augments a given set of training images by applying 
%a random slice dropout
%   [X,D]=AUGMENTSLICEDROPOUT(X,P,AB)
%   * X is the original data
%   * P is the probability to apply a certain dropout to a slice
%   * AB are the parameters of the beta distribution for dropouts, if both
%   1, uniform distribution in 0-1, if a lower/greater than 1, then 
%   grows/drops to inf/0 at 0, if b lower/greater than 0 then grows/drops 
%   to inf/0 at 1. Also cases 2-1/1-2 are linear pdfs from 0 to 2 / 2 to 0,
%   default is 1.5/1
%   * X is the augmented data
%   * S is the set of applied scalings
%

if nargin<2 || isempty(p);p=0.5;end
if nargin<3 || isempty(ab);ab=[1.5 1];end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);

d=single(betarnd(ab(1),ab(2),[1 1 N(3) 1 NSamp]).*binornd(1,p,[1 1 N(3) 1 NSamp]));
if gpu;d=gpuArray(d);end
d(d==0)=1;
x=bsxfun(@times,x,d);


