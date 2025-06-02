function x=linearRigidTransform(x,T,rGrid,di,sumR)

%LINEARRIGIDTRANSFORM rigidly transforms volumes using linear-based interpolation (both forwards and backwards)
%   [X,XB]=LINEARRIGIDTRANSFORM(X,T,{DI},{SUMR})
%   * X is a volume
%   * T are the transform factors
%   * RGRID is the input spatial grid of size 3xNP
%   * {DI} is a flag to indicate whether to perform direct (default) or inverse transform
%   * {SUMR} indicates whether to sum along the motion states when applying the inverse transform (defaults to 1)
%   * X is the rigidly transformed volume
%

gpu=isa(x,'gpuArray');if gpu;gpuF=1;else gpuF=0;end

if ~exist('di','var') || isempty(di);di=1;end
if ~exist('sumR','var') || isempty(sumR);sumR=1;end

%CONTROL VARIABLES
tr=[1 3 2;
    2 1 3
    3 2 1];
NDT=numDims(T);
NT=size(T);
NTP=[NT(1:3) prod(NT(4:NDT-1)) NT(NDT)];
T=reshape(T,NTP);
NX=size(x);NX(end+1:NDT-1)=1;
NXP=[NX(1:3) prod(NX(4:end))];
x=repmat(x,[ones(1,3) NTP(4)/NXP(4)]);

%ARRANGE THE CELL STRUCTURE
for n=1:3
    NR=NX;NR(n)=1;
    rGrid{n}=repmat(rGrid{n},NR);
    rGrid{n}=rGrid{n}(:)';
end
rGrid{4}=rGrid{3};rGrid{4}(:)=1;
rGrid=cat(1,rGrid{:});

%BUILD HOMOGENEOUS TRANSFORMS---THIS MAY BE MOVED TO PRECOMPUTEFACTORS IF
%USED FOR RECONSTRUCTION
R=cell(1,3);
for m=1:3
    R{m}=single(zeros([3 3 1 NT(4)]));    
    if gpu;R{m}=gpuArray(R{m});end
    th=dynInd(T,m+3,5);cth=cos(th);sth=sin(th);
    R{m}=dynInd(R{m},[tr(1,m) tr(1,m)],1:2,cth);R{m}=dynInd(R{m},[tr(1,m) tr(2,m)],1:2,-sth);
    R{m}=dynInd(R{m},[tr(2,m) tr(1,m)],1:2,sth);R{m}=dynInd(R{m},[tr(2,m) tr(2,m)],1:2,cth);
    R{m}=dynInd(R{m},[tr(3,m) tr(3,m)],1:2,1);
end
A=single(zeros([4 4 1 NTP(4)]));
if gpu;A=gpuArray(A);end
A=dynInd(A,{1:3,1:3},1:2,matfun(@mtimes,R{3},matfun(@mtimes,R{2},R{1})));

perm=1:5;perm(5)=1;perm(1)=5;
A=dynInd(A,{1:3,4},1:2,permute(dynInd(T,1:3,5),perm));
A=dynInd(A,[4 4],1:2,1);
if di;A=matfun(@inv,A);end%This may be accelerated probably

rIGrid=cell(1,3);rIoGrid=cell(1,3);
for m=1:3;rIGrid{m}=reshape(rGrid(m,:),NX(1:3));end

%INTERPOLATE
for n=1:NTP(4)
    roGrid=dynInd(A,n,4)*rGrid;
    for m=1:3;rIoGrid{m}=reshape(roGrid(m,:),NX(1:3));end
    %if di;x=dynInd(x,n,4,interpn(rIoGrid{1},rIoGrid{2},rIoGrid{3},dynInd(x,n,4),rIGrid{1},rIGrid{2},rIGrid{3},'linear',0));
    %else end
    x=dynInd(x,n,4,interpn(rIGrid{1},rIGrid{2},rIGrid{3},dynInd(x,n,4),rIoGrid{1},rIoGrid{2},rIoGrid{3},'linear',0));
end
if sumR;x=sum(x,4);else x=reshape(x,[NX(1:3) NT(4:NDT-1)]);end
