function x=augmentSlicePermuting(x,p,NSl)

%AUGMENTSLICEPERMUTING   Augments a given set of training images by
%applying a random permutation in groups of NSl slices
%   [X,T]=AUGMENTSLICEPERMUTING(X,P,NSL)
%   * X is the original data
%   * P is the probability that a given group of slices get permuted
%   * NSL is the group of slices for permutation
%   * X is the augmented data
%

if nargin<2 || isempty(p);p=0.1;end
if nargin<3 || isempty(NSl);NSl=2;end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
NChan=N(4);
N=N(1:3);

indSl=(1:N(3))';
NBlocks=floor(N(3)/NSl);
NUse=NBlocks*NSl;
NTotalSamples=prod(NSamp);

indSl=resampling(indSl,NUse,2);
indSl=reshape(indSl,[NSl NBlocks]);
perm=cell2mat(arrayfun(@(x) randperm(NSl),(1:NBlocks*NTotalSamples)','un',0));
perm=permute(perm,[2 1]);
perm=reshape(perm,[NSl NBlocks NTotalSamples]);
indSl=repmat(indSl,[1 1 NTotalSamples]);
indSlApply=indDim(indSl,perm,1);
r=binornd(1,p*ones([1 NBlocks NTotalSamples]));
indSlApply=bsxfun(@times,(1-r),indSl)+bsxfun(@times,r,indSlApply);
indSl=reshape(indSl,[NSl*NBlocks NTotalSamples]);
indSlApply=reshape(indSlApply,[NSl*NBlocks NTotalSamples]);

x=reshape(x,[N NChan NTotalSamples]);
for n=1:NTotalSamples;x(:,:,indSl(:,n),:,n)=x(:,:,indSlApply(:,n),:,n);end
x=reshape(x,[N NChan NSamp]);