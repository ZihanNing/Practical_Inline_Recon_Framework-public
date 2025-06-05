function [W,WH,I]=buildWaveletMatrix(N,wname,L,packet,gpu)

%BUILDWAVELETMATRIX  builds a set of separable wavelet matrices for
%multidimensional wavelet decomposition. This code is based on 
%https://www.mathworks.com/matlabcentral/fileexchange/65476-n-level-wavelet-matrix
%   [W,WH]=BUILDWAVELETMATRIX(N,{WNAME},{L},{PACKET},{GPU})
%   * N are the dimensions of the space
%   * {WNAME} is the wavelet name, it defaults to 'db4'
%   * {L} are the number of levels for wavelet decomposition. It defaults
%   to 1
%   * {PACKET} is a flag to perform a wavelet packet decomposition instead
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   ** W is the cell of wavelet decomposition matrices
%   ** WH is the cell of wavelet synthesis matrices
%   ** I is an array with the decomposition levels
%

ND=length(N);

if nargin<2 || isempty(wname);wname='db4';end
if nargin<3 || isempty(L);L=ones(1,ND);end
if nargin<4 || isempty(packet);packet=0;end
if nargin<5 || isempty(gpu);gpu=single(gpuDeviceCount) && ~blockGPU;end

if length(L)==1;L=repmat(L,[1 ND]);end

if strcmp(wname,'dct');packet=1;end

W=cell(1,ND);
if nargout>1;WH=cell(1,ND);end
if nargout>2 && ~packet
    I=ones(1,'single');
    if gpu;I=gpuArray(I);end
end

idCont=1;
for n=1:ND
    if nargout>2 && ~packet
        NN=ones(1,ND+1);NN(n)=N(n);
        I=repmat(I,NN);
    end
    if ~strcmp(wname,'dct')
        if iscell(wname);[ha,hd]=wfilters(wname{n});
        else [ha,hd]=wfilters(wname);
        end
        [ha,hd]=parUnaFun({ha,hd},@single);
        H=cat(2,permute(ha,[1 3 2]),permute(hd,[1 3 2]));

        Wn=eye(N(n),'single');
        for l=1:L(n)
            if ~packet
                temp=eye(N(n),'single');       
                temp(1:N(n)/2^(l-1),1:N(n)/2^(l-1))=waveletMatrix(N(n)/2^(l-1),H);
            else
                temp=[];
                aux=waveletMatrix(N(n)/2^(l-1),H);
                for s=1:2^(l-1);temp=blkdiag(temp,aux);end
            end       
            Wn=temp*Wn;
            if nargout>2 && ~packet
                tempI=zeros(NN,'single');
                tempI(1:NN(n)/2^l)=1;
                if gpu;tempI=gpuArray(tempI);end          
                I=bsxfun(@plus,I,idCont*tempI);
            end
        end 
        idCont=idCont*(L(n)+1);
    else
        Wn=single(dctmtx(N(n)));
    end        
    if gpu;Wn=gpuArray(Wn);end
    W{n}=Wn;
    if nargout>1;WH{n}=Wn';end
end
if nargout>2 && packet && ~strcmp(wname,'dct')
    I=single(1:prod(2.^L));
    I=reshape(I,[2.^L 1]);
    perm=1:2*ND;perm(1:2:2*ND)=(ND+1):2*ND;perm(2:2:2*ND)=1:ND;
    I=permute(I,perm);
    NN=ones(1,2*ND);NN(1:2:2*ND)=N./(2.^L);
    I=repmat(I,NN);
    I=reshape(I,[N 1]);
    if gpu;I=gpuArray(I);end
end
if strcmp(wname,'dct')
    IG=generateGrid(N,gpu,pi*ones(1,3),ones(1,3));
    for n=1:ND
        NN=[N 1];NN(n)=1;
        IG{n}=repmat(IG{n},NN);
    end
    IG=cat(ND+1,IG{:});
    IG=sum(abs(IG).^2,ND+1);
    [~,indS]=sort(IG(:));
    indS=reshape(indS,[],prod(2.^L));
    I=zeros([prod(N) 1],'like',real(W{1}));
    for s=1:prod(2.^L);I(indS(:,s))=s;end
    I=reshape(I,[N 1]);
end

function W=waveletMatrix(N,H)

W=zeros([N/2 2 N],'single');
assert(size(W,3)>=size(H,3),'Matrix not orthogonal, try to increase N or reduce L');
W(1,:,1:size(H,3))=flip(H,3);
for m=2:N/2;W(m,:,:)=circshift(W(m-1,:,:),[0 0 2]);end
W=reshape(W,[N N]);
