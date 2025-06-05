function [x,f]=augmentFlip(x,p)

%AUGMENTFLIP   Augments a given set of training images by applying a 
%random axis flip
%   [X,F]=AUGMENTFLIP(X,P)
%   * X is the original data
%   * P is the probability that a given axis gets flipped
%   * X is the augmented data
%   * F is the set of applied flips
%

if nargin<2 || isempty(p);p=0.5;end
if length(p)==1;p=p*ones(1,3);end

N=size(x);N(end+1:5)=1;
NSamp=N(5:end);NSampTotal=prod(NSamp);
NChan=N(4);
N=N(1:3);

f=cell(1,3);
x=reshape(x,[N NChan NSampTotal]);
for n=1:3
    f{n}=binornd(1,p(n),[ones(1,4) NSamp]);
    x(:,:,:,:,f{n}(:)==1)=flip(x(:,:,:,:,f{n}(:)==1),n);
end
x=reshape(x,[N NChan NSamp]);
