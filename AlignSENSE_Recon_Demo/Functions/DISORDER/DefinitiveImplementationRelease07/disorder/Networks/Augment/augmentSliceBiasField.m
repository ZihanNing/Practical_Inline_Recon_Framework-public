function [x,b]=augmentSliceBiasField(x,p,order,strength)

%AUGMENTSLICEBIASFIELD   Augments a given set of training images by 
%applying a random (spin-history induced) bias field to slices
%   [X,B]=AUGMENTSLICEBIASFIELD(X,P,ORDER,STRENGTH)
%   * X is the original data
%   * P is the probability to apply a certain bias field to a slice
%   * ORDER is the order of the field in each dimension of space. It
%   defaults to 1
%   * STRENGTH is the strength of the field. It defaults to 0.2 (note this 
%   is applied to exponentials, so the effect grows very quickly)  
%   ** X is the augmented data
%   ** B are the generated fields
%

if nargin<2 || isempty(p);p=0.5;end
if nargin<3 || isempty(order);order=1;end%[1 4 1] could be interesting to model high variability in the PE direction
if nargin<4 || isempty(strength);strength=0.2;end%

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
N=N(1:3);
if length(order)==1;order=order*ones(1,2);end

rGrid=generateGrid(N(1:2),gpu,N(1:2),ones(1,2));
for n=1:2;rGrid{n}=(rGrid{n}/(order(n)+eps)).^2;end
idx=bsxfun(@plus,rGrid{1},rGrid{2})<=1;
NP=gather(sum(idx(:)));
rc=strength*randn([NP prod(NSamp)*N(3)]);

b=zeros([prod(N(1:2)) prod(NSamp)*N(3)],'like',x);
b(idx,:)=rc;
b=reshape(b,[N 1 NSamp]);
b=b*sqrt(prod(N(1:2)));
for n=1:2
    if N(n)~=1;b=ifctGPU(b,n);end
end
r=single(binornd(1,p,[1 1 N(3) 1 NSamp]));
if gpu;r=gpuArray(r);end
b=bsxfun(@times,r,b);
b=min(b,0);
b=exp(b);
x=bsxfun(@times,x,b);



