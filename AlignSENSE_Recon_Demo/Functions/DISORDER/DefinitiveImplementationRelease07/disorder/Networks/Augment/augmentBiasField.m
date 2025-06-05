function [x,b]=augmentBiasField(x,order,strength)

%AUGMENTBIASFIELD   Augments a given set of training images by applying a
%random bias field
%   [X,B]=AUGMENTBIASFIELD(X,ORDER,STRENGTH)
%   * X is the original data
%   * ORDER is the order of the field in each dimension of space. It
%   defaults to 2
%   * STRENGTH is the strength of the field. It defaults to 0.1 (note this 
%   is applied to exponentials, so the effect grows very quickly)  
%   ** X is the augmented data
%   ** B are the generated fields
%

if nargin<2 || isempty(order);order=2;end%[1 4 1] could be interesting to model high variability in the PE direction
if nargin<3 || isempty(strength);strength=0.1;end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);NSampTotal=prod(NSamp);
N=N(1:3);
if length(order)==1;order=order*ones(1,3);end

rGrid=generateGrid(N,gpu,N,ones(1,3));
for n=1:3;rGrid{n}=(rGrid{n}/(order(n)+eps)).^2;end
idx=bsxfun(@plus,bsxfun(@plus,rGrid{1},rGrid{2}),rGrid{3})<=1;
NP=gather(sum(idx(:)));
rc=strength*randn([NP NSampTotal]);

b=zeros([prod(N) NSampTotal],'like',x);
b(idx,:)=rc;
b=reshape(b,[N 1 NSamp]);
b=b*sqrt(prod(N));
for n=1:3
    if N(n)~=1;b=ifctGPU(b,n);end
end
b=exp(b);
x=bsxfun(@times,x,b);