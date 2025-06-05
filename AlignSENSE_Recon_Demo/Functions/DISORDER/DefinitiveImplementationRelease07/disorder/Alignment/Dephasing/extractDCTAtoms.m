function x=extractDCTAtoms(N,M,FOV,gpu)

%EXTRACTDCTATOMS   Extracts the lowest atoms of the DCT
%   X=EXTRACTDCTATOMS(N,M,FOV)
%   * N are the dimensions of the space
%   * M is the number of atoms
%   * {FOV} are the FOV factors
%   * {GPU} determines whether to use gpu computations 
%   * X are the extracted atoms
%

ND=length(N);
if nargin<3 || isempty(FOV);FOV=ones(1,ND);end
if nargin<4 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end

%INDEXES
xx=generateGrid(N,gpu,N./FOV,ones(1,ND));
r=xx{1}(1);
for n=1:ND;r=bsxfun(@plus,r,xx{n}.^2);end
[~,ir]=sort(r(:));
irM=ir(1:M);
irs=ind2subV(N,irM);

for n=1:ND
    D{n}=single(dctmtx(N(n)));
    if gpu;D{n}=gpuArray(D{n});end
end

x=zeros([N M],'like',D{1});
for m=1:M
    r=xx{1}(1);r=r+1;
    for n=1:ND
        perm=1:ND;perm([2 n])=perm([n 2]);
        r=bsxfun(@times,r,permute(D{n}(irs(m,n),:),perm));
    end
    x=dynInd(x,m,ND+1,r);
end