function [X,amse,U,S,V]=heterogeneousSVShrinkage(Y,L,I,H,MNeff,orth)

%HETEROGENEOUSSVSHRINKAGE  applies a shrinkage of the singular values of Y 
%according to [1] W Leeb, "Matrix denoising for weighted loss functions 
%and heterogeneous signals", arXiv, 1902.0947v1, 2019
%   [AMSE,X]=HETEROGENEOUSSVSHRINKAGE(Y,LR,LN,IM,IN)
%   * Y is the matrix to be denoised
%   * {L} are noise generation matrices for the rows (size Mx1 or MxW) and
%   the columns (size Nx1 or NxW)
%   * {I} are indicator functions for localizing along specific indexes
%   on the rows, (size MxI with I the number of subspaces). The columns
%   should add up to 1 (orthogonality would be another condition, but it is
%   not checked) / columns, (size JxN with J the number of subspaces). The 
%   rows should add to 1 (orthogonality would be another condition, but it 
%   is not checked)
%   * {H} are (orthogonal) decompositions for the rows / columns 
%   * {MNEFF} are the effective dimensions of the row /column space
%   * {ORTH} indicates whether to make the orthogonal decomposition
%   approximation when computing the traces
%   ** X is the estimated matrix
%   ** AMSE is the estimated asymptotic mean square error (per subspace)
%   ** U are the estimated left singular vectors
%   ** S are the estimated singular values
%   ** V are the estimated right singular vectors
%

[M,N]=size(Y);
MN=[M,N];
if nargin<2;L=cell(1,2);end
if nargin<3 || isempty(I);I={ones(M,1),ones(N,1)};end
for co=1:2
    if isempty(I{co});I{co}=ones(MN(co),1);end
end
if nargin<4;H=cell(1,2);end
if nargin<5;MNeff={M,N};end
if nargin<6 || isempty(orth);orth=0;end
for co=1:2
    if isempty(MNeff{co});MNeff{co}=MN(co);end
end

%WE FLIP THE ENTRIES
flEn=0;
if M>N
    flEn=1;
    Y=Y';
    L(2:-1:1)=L(1:2);I(2:-1:1)=I(1:2);H(2:-1:1)=H(1:2);MNeff(2:-1:1)=MNeff(1:2);
    [M,N]=size(Y);
    MN=[M,N];
end

gpu=isa(Y,'gpuArray');sing=isa(Y,'single');
Y=double(Y);
for co=1:2
    I{co}=double(I{co});
    for n=1:length(L{co});L{co}{n}=double(L{co}{n});end
    for n=1:length(H{co});H{co}{n}=double(H{co}{n});end
end
%NOISE WHITENING
for co=1:2;Y=whiten(Y,L{co},1,co,MNeff{co});end
%SVD
comp=~isreal(Y);
Ynorm=sqrt(N*(1+comp));
[S,U{1},U{2}]=svdm(Y/Ynorm);%S column vector

%WEIGHTED VECTORS
for co=1:2;Uw{co}=whiten(U{co},L{co},0,1,MNeff{co});end

%WE HAVE IMPLEMENTED AN APPROXIMATION FOR BIG MATRICES
for co=1:2
    if MN(co)<=2^12 && ~orth
        WL{co}=whiten(eye(MN(co)),L{co},0,1,MNeff{co});
    else        
        if ~isempty(L{co});WL{co}=L{co}{1};else WL{co}=[];end
    end
end
dev=gpuDevice;

%SHRINKAGE OVER THE DIFFERENT SUBSPACES
IJ=zeros(1,2);for co=1:2;IJ(co)=size(I{co},2);end
X=zeros([M N],'like',Y);
amse=zeros(1,'like',Y);
Ii=cell(1,2);Uiw=cell(1,2);WLi=cell(1,2);mu2=zeros(1,2);
for i=1:IJ(1);co=1;
    Ii{co}=I{co}(:,i)==1;
    if any(Ii{co})
        [Uiw{co},WLi{co},mu2(co)]=buildUmu(H{co},Uw{co},WL{co},Ii{co},MN(co),MNeff{co});
        for j=1:IJ(2);co=2;
             Ii{co}=I{co}(:,j)==1;
             if any(Ii{co})
                [Uiw{co},WLi{co},mu2(co)]=buildUmu(H{co},Uw{co},WL{co},Ii{co},MN(co),MNeff{co});
                
                %0-> unweighted, 1-> asymptotically orthogonal, 2->general
                ca=0;for co=1:2;ca=ca+single(~(isvector(WLi{co}) && all(WLi{co}==1)));end

                %SINGULAR VALUE ESTIMATION
                [Sij,amseij,Rij]=weightedShrinkage(gather(S),gather(Uiw{1}),gather(Uiw{2}),M/N,mu2(1),mu2(2),ca);
                if gpu;Sij=gpuArray(Sij);amseij=gpuArray(amseij);end

                %SIGNAL SYNTHESIS
                enRat=(sum(Ii{1})*sum(Ii{2}))/Ynorm^2;
                X=X+bsxfun(@times,Uiw{1},Ynorm*Sij)*Uiw{2}';
                amse=amse+amseij/enRat;%NOTE THIS MAY BE SPLIT SPATIALLY              
            end
        end
    end
end

if nargout>=3;[S,U,V]=svdm(X/Ynorm);end%S column vector

if sing
    [X,amse]=parUnaFun({X,amse},@single);
    if nargout>=3;[S,U,V]=parUnaFun({S,U,V},@single);end
end

if gpu
    [X,amse]=parUnaFun({X,amse},@gpuArray);
    if nargout>=3;[S,U,V]=parUnaFun({S,U,V},@gpuArray);end
end

if flEn
    X=X';amse=amse';
    if nargout>=3;aux=U;U=V';V=aux';end
end

function [Uiw,WLi,mu2]=buildUmu(H,Uw,WL,Ii,MN,MNeff)
    Uiw=whiten(Uw,H,[0 1],1,MNeff,Ii,1);    
    if MN<=2^12 && ~orth        
        WLi=whiten(WL,H,[0 1],1,MNeff,Ii,1);
        mu2=gather(sum(abs(WLi(:)).^2)/MN);%Normalized squared trace
    else%WE ASSUME ORTHOGONAL DECOMPOSITION
        if isempty(WL)
            mu2=gather(sum(Ii)/MN);
            WLi=1;
        else
            mu2=gather(sum(abs(WL).^2)*sum(Ii)/MN^2);
            WLi=WL;
        end
        %mu2=sum(abs(WL(Ii)).^2)/MN;%This is if one can operate spatially      
    end
end

end