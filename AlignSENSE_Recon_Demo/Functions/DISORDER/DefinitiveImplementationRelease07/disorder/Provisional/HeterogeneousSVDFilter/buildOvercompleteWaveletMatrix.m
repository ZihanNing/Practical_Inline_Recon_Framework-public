function [W,WH,I]=buildOvercompleteWaveletMatrix(N,wname,L,gpu)

%BUILDOVERCOMPLETEWAVELETMATRIX  builds a set of overcomplete separable 
%wavelet matrices for multidimensional wavelet decomposition. This code is 
%based on https://www.mathworks.com/matlabcentral/fileexchange/65476-n-level-wavelet-matrix
%   [W,WH]=BUILDOVERCOMPLETEWAVELETMATRIX(N,{WNAME},{L},{GPU})
%   * N are the dimensions of the space
%   * {WNAME} is the wavelet name, it defaults to 'db4'
%   * {L} are the number of levels for wavelet decomposition. It defaults
%   to 1
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   ** W is the cell of wavelet decomposition matrices
%   ** WH is the cell of wavelet synthesis matrices
%   ** I is an array with the decomposition levels, not implemented
%

ND=length(N);

if nargin<2 || isempty(wname);wname='db4';end
if nargin<3 || isempty(L);L=ones(1,ND);end
if nargin<4 || isempty(gpu);gpu=single(gpuDeviceCount) && ~blockGPU;end

if length(L)==1;L=repmat(L,[1 ND]);end
I=[];

W=cell(1,ND);
if nargout>1;WH=cell(1,ND);end

for n=1:ND
    if iscell(wname);[ha,hd]=wfilters(wname{n});
    else [ha,hd]=wfilters(wname);
    end
    [ha,hd]=parUnaFun({ha,hd},@single);    
    H=cat(1,ha,hd);
    Wn=zeros([N(n) N(n) 2 L(n)],'single');      
  
    for l=1:L(n);Wn(:,:,:,l)=waveletMatrix(N(n),H,l);end 
    %To match order of swt, first detail from level 1 to L, then
    %approximation from level 1 to L
    Wn=permute(Wn,[1 2 4 3]);Wn=flip(Wn,4);
    
    Wn=permute(Wn,[1 3 4 2]);
    Wn=reshape(Wn,[2*L(n)*N(n) N(n)]);

    if gpu;Wn=gpuArray(Wn);end
    W{n}=Wn;
    if nargout>1;WH{n}=mldivide(Wn'*Wn,Wn');end
end

function W=waveletMatrix(N,H,L)
ha=H(1,:);hd=H(2,:);haor=ha;
NN=size(ha,2)/2;
for l=1:L-1 
    hd=conv(dyadup(hd,'m',0),haor);
    ha=conv(dyadup(ha,'m',0),haor);
end
H=cat(1,ha,hd);
H=permute(H,[3 1 2]);
W=zeros([N 2 N],'single');
assert(size(W,3)>=size(H,3),'Matrix not orthogonal, try to increase N or reduce L');
W(1,:,1:size(H,3))=flip(H,3);
W(1,:,:)=circshift(W(1,:,:),[0 0 -((NN-1)*2^L-(NN-1))]);
for m=2:N;W(m,:,:)=circshift(W(m-1,:,:),[0 0 1]);end
W=permute(W,[1 3 2]);
