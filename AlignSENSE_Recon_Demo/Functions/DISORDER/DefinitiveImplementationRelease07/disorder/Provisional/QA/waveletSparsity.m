function waSp=waveletSparsity(x,db)

%WAVELETSPARSITY  Computes the wavelet sparsity of a 3D array
%   WASP=WAVELETSPARSITY(X,{DB})
%   * X is the 3D array
%   * {DB} is the Daubechies wavelet type, it defaults to 4
%   ** GREN is the gradient entropy
%

if nargin<2 || isempty(db);db=4;end

po=0.25;%l1 norm

wty=sprintf('db%d',db);         
NL=3;

%N=size(x);
%gpu=isa(x,'gpuArray');
%W=buildWaveletMatrix(N,wty,NL,[],gpu);
%x=waveletTransform(x,W);
%waSp=gather(sum(abs(x(:)).^po));

xd=wavedec3(x,NL,wty);
waSp=0;
for w=1:length(xd.dec);waSp=waSp+sum(abs(xd.dec{w}(:)).^po);end
