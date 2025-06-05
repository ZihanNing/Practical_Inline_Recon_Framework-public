function [x,r]=augmentSliceRotation(x,p,range,cent)

%AUGMENTSLICEROTATION   Augments a given set of training images by 
%applying a random rotation independently for each slice
%   [X,R]=AUGMENTSLICEROTATION(X,{P},{RANGE},{CENT})
%   * X is the original data
%   * {P} is the probability that a given slice gets rotated
%   * {RANGE} describes the range of random rotations. It defaults to pi/8
%   * {CENT} describes the center of rotation. It defaults to (N/2)+1
%   ** X is the augmented data
%   ** R is the set of applied translations
%

if nargin<2 || isempty(p);p=0.5;end
if nargin<3 || isempty(range);range=pi/8;end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);
comp=~isreal(x);

rndMot=single(binornd(1,p*ones([1 1 N(3) NSamp])));
if gpu;rndMot=gpuArray(rndMot);end

rot=bsxfun(@times,rndMot,range*(rand([1 1 N(3) NSamp],'like',x)-0.5));
rot=reshape(rot,[1 1 1 1 N(3)*NSamp]);
r=zeros([1 1 1 1 N(3)*NSamp 6],'like',rot);
r(1,1,1,1,:,4)=rot;

NT=N;
NT(3)=1;
if nargin<4 || isempty(cent);cent=(NT/2)+1;end
cent(end+1:3)=(NT(3)/2)+1;

[~,kGrid,rkGrid,~,cGrid]=generateTransformGrids(NT,gpu,NT,cent);
et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,r,[],[],[],[],cGrid);
x=reshape(x,[N NChan prod(NSamp)]);
x=permute(x,[1 2 4 3 5]);
x=reshape(x,[N(1:2) 1 NChan N(3)*prod(NSamp)]);
x=sincRigidTransform(x,et);
x=reshape(x,[N(1:2) NChan N(3) prod(NSamp)]);
x=permute(x,[1 2 4 3 5]);
x=reshape(x,[N NChan NSamp]);
if ~comp;x=real(x);end