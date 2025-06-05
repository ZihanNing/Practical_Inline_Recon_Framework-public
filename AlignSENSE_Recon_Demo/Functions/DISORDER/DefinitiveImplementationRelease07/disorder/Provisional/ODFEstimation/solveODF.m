function x=solveODF(x,b,order,M,no,lambda)

%SOLVEODF   Estimates the ODFs up to a given order by nnLS using the diffusion
%measures. It adapts the methods in the fanDTasia ToolBox, from
%https://uk.mathworks.com/matlabcentral/fileexchange/26997-fandtasia-toolbox?requestedDomain=www.mathworks.com,
%which uses [1] A Barmpoutis and BC Vemuri, "A Unified Framework for 
%Estimating Diffusion Tensors of any order with Symmetric 
%Positive-Definite Constraints," Proc IEEE ISBI, 1385-1388, 2010.
%   Y=SOLVEODF(X,B,{ORDER},{M},{NO},{LAMBDA})
%   * X is the DWI data
%   * B are the corresponding diffusion measures
%   * {ORDER} is the ODF estimation order, defaults to 2 (DTI)
%   * {M} is a mask to restrict the computation to a specific ROI, defaults
%   to all 1s
%   * {NO} indicates the norm to be minimized for estimation (2 for LS, 1 for LAD). If singleton, same norm is applied for all bvals. Otherwise,
%   first component is used for 0 b-vals and the other components for non 0 b-vals
%   * {LAMBDA} is a Tikhonov regularization termp, defaults to 0
%   * X are the estimated ODFs
%

if ~exist('order','var') || isempty(order);order=2;end
if ~exist('lambda','var') || isempty(lambda);lambda=0;end
if ~exist('no','var') || isempty(no);no=[2 2];end
useM=exist('M','var') && ~isempty(M);

gpu=isa(b,'gpuArray');
NX=size(x);NDX=numDims(x);
NB=size(b);NDB=numDims(b);
assert(NDX==4,'The data should be 4th-dimensional, with different acquired b-vals arranged along the 4th dimension, but it has %d dimensions',NDX);
assert(NDB==5,'The gradients should be 5th-dimensional, with different acquired b-vals arranged along the 4th dimension, and their values along the 5th, but they have %d dimensions',NDB);
b=repmat(b,NX./NB(1:NDX));%A table of b-vals per pixel
NB=size(b);

%We operate on a pixelwise basis
x=reshape(x,[prod(NX(1:3)) NX(4)]);
b=reshape(b,[prod(NB(1:3)) NB(4:5)]);

if useM;x=x(logical(M(:)),:);b=b(logical(M(:)),:,:);end%Only over the mask
bun=unique(gather(b(1,:,4)));
bun(bun==0)=[];
%Construct set of polynomial coefficients C
C=constructSetOf321Polynomials(order)'; %computes C from section 5.1 (ISBI'10)
D=constructMatrixHighOrderTensorToSphericalHarmonic(order);
UnitVectors;
if gpu;g=gpuArray(g);end
NB=length(bun);
NC=size(C,1);
[x,b]=parUnaFun({x,b},@permute,[2 3 1]);

%B0 calculations
xb=x(b(:,4,1)==0,:,:);
xb=double(gather(xb));
yM=xb(1,:,:);
A=ones(size(xb,1),1);
noi=no(1);
parfor n=1:size(xb,3);yM(:,:,n)=LADODF(A,xb(:,:,n),noi);end
if gpu;yM=gpuArray(yM);end
%yM=mean(x(:,b(1,:,4)==0),2);%B0
x=bsxfun(@times,x,1./yM);%S/S0
yM=repmat(yM,[NB*NC+1 1 1]);

dev=gpuDevice;
wait(dev);tic

for bval=1:length(bun)
    xb=x(b(:,4,1)==bun(bval),:,:);
    bb=b(b(:,4,1)==bun(bval),1:3,:);
    fprintf('%d voxels for bval %.3f\n',size(bb,3),bun(bval)/1000);   
    bS=1e4;
    O=size(xb,3);
    for o=1:bS:O;vO=o:min(o+bS-1,O);        
        if ~isempty(vO)  
            [bbo,xbo]=parUnaFun({bb,xb},@dynInd,vO,3);        
            B=constructMatrixOfIntegralsV(bbo,order,100,g);
            B=matfun(@mtimes,B,C);
            xbr=double(zeros([size(C,2) 1 size(B,3)]));
            [B,xbo]=parUnaFun({B,xbo},@double);
            [B,xbo]=parUnaFun({B,xbo},@gather);
            noi=no(2);
            parfor n=1:length(vO);xbr(:,:,n)=LADODF(B(:,:,n),xbo(:,:,n),noi);end
            xbr=single(xbr);
            if gpu;xbr=gpuArray(xbr);end
            yM(2+NC*(bval-1):1+NC*bval,:,vO)=matfun(@mtimes,D*C,xbr);
        end
    end
end

wait(dev);toc

yM=permute(yM,[3 1 2]);
NO=size(yM,2);
y=single(zeros([prod(NX(1:3)) NO]));
if gpu;y=gpuArray(y);end
if useM;y(logical(M),:)=yM;else y=yM;end
y=reshape(y,[NX(1:3) NO]);
x=[];
x=cell(1,length(bun)+1);
x{1}=dynInd(y,1,4);
for bval=1:length(bun);x{bval+1}=dynInd(y,2+NC*(bval-1):1+NC*bval,4);end
