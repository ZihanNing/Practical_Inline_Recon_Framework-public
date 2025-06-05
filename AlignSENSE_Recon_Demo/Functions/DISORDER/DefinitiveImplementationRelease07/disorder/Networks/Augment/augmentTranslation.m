function [x,t]=augmentTranslation(x,Nrange,zf)

%AUGMENTTRANSLATION   Augments a given set of training images by applying a 
%random integer translation
%   [X,T]=AUGMENTTRANSLATION(X,NRANGE)
%   * X is the original data
%   * NRANGE is the range of pixels for random uniform translation
%   * X is the augmented data
%   * T is the set of applied translations
%

if nargin<3 || isempty(zf);zf=0;end

N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);

if nargin<2 || isempty(Nrange);Nrange=N;end%Uniform location at random
if length(Nrange)==1;Nrange=Nrange*ones(1,3);end

t=cell(1,3);
for n=1:3
    t{n}=randi(round(Nrange(n)),[ones(1,4) NSamp],'like',real(x))-round(Nrange(n)/2);
    if N(n)==1;t{n}(:)=0;end
end
x=shifting(x,t,[],zf);
