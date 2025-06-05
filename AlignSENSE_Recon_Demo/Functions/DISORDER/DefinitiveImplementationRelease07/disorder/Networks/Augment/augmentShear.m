function [x,s]=augmentShear(x,only,slope)

%AUGMENTSHEAR   Augments a given set of training images by applying a 
%random shear 
%   [X,S]=AUGMENTSHEAR(X,ONLY,SLOPE)
%   * X is the original data
%   * ONLY applies the shear only along a specific axis
%   * SLOPE serves to select the shear slope. It defaults to 1 (meaning
%   that the maximum difference in translations is of 1 pixel per pixel
%   * X is the augmented data
%   * S is the set of applied shears
%

if nargin<3 || isempty(slope);slope=1;end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);

comp=~isreal(x);
per=[1 2 3 1 2 3;
     2 1 1 3 3 2];

if nargin>=2 && ~isempty(only)
    dper=sum(per==only,1);
    per=per(:,~dper);
end
NSh=size(per,2);
she=slope*(rand([ones(1,4) NSamp NSh],'like',x)-0.5);
for n=1:NSh;x=shearing(x,dynInd(she,n,6),per(:,n));end
if ~comp;x=real(x);end
