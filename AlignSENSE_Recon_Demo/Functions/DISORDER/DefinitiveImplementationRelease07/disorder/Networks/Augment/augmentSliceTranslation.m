function [x,t]=augmentSliceTranslation(x,p,strength)

%AUGMENTSLICETRANSLATION   Augments a given set of training images by 
%applying a random translation independently for each slice
%   [X,T]=AUGMENTSLICETRANSLATION(X,P,STRENGTH)
%   * X is the original data
%   * P is the probability that a given slice gets translated
%   * STRENGTH is the strength of the translation
%   * X is the augmented data
%   * T is the set of applied translations
%

if nargin<2 || isempty(p);p=0.5;end
if nargin<3 || isempty(strength);strength=4;end
if length(strength)==1;strength=strength*ones(1,2);end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);

rndMot=single(binornd(1,p*ones([1 1 N(3) 1 NSamp])));
if gpu;rndMot=gpuArray(rndMot);end

t=cell(1,2);
for n=1:2
    t{n}=bsxfun(@times,rndMot,strength(n)*randn([1 1 N(3) 1 NSamp],'like',real(x)));
    if N(n)==1;t{n}(:)=0;end
end
x=shifting(x,t);
