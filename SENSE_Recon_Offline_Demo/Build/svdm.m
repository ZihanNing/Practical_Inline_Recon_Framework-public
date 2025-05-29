function [S,U,V]=svdm(X,NC,useParallel)

%SVDM   Performs accelerated multiple svds for wide and square matrices 
%based on a Schur decomposition of the normal matrix X*X'
%   [S,U,V]=SVDM(X,{NC},{USEPARALLEL})
%   * X is the input matrix
%   * {NC} indicates to return a given number of components
%   * {USEPARALLEL} indicates whether to use parallel computing, it
%   defaults to 1
%   ** S are the singular values of X
%   ** U is the column space (range) of X
%   ** V is the row space (range) space of X
%

if nargin<2;NC=[];end
if nargin<3 || isempty(useParallel);useParallel=1;end

gpu=isa(X,'gpuArray');
NX=size(X);NX(end+1:3)=1;
X=resSub(X,3:length(NX));
[M,N,O]=size(X);
if nargout>1
    Xp=matfun(@ctranspose,X);
    U=matfun(@mtimes,X,Xp);
    U=(U+matfun(@ctranspose,U))/2;%Force the matrix to be Hermitian  
    if isempty(NC) || nargout>2
        U=gather(U);S=U;  
        if O>=8 && useParallel
            parfor (o=1:O,32);[U(:,:,o),S(:,:,o)]=schur(U(:,:,o));end
        else
            for o=1:O;[U(:,:,o),S(:,:,o)]=schur(U(:,:,o));end
        end
        if gpu;[U,S]=parUnaFun({U,S},@gpuArray);end
        S=diagm(sqrt(abs(S)));
        S=flip(S,2);%We use decreasingly sorted singular values
        U=flip(U,2);
        [S,iS]=sort(S,2,'descend');  
        U=indDim(U,iS,2);
        D=S;D(D<1e-9)=1;
        if nargout>2;V=bsxfun(@times,matfun(@mtimes,Xp,U),1./D);end%PROBLEM WITH SINGLETONS IS THAT THIS MAY DEPART A LOT FROM ORTHONORMALITY!!
    else
        U=double(gather(U));
        UU=zeros([M NC O],'like',U);
        S=zeros([NC NC O],'like',U);
        if O>=8 && useParallel
            parfor (o=1:O,32);[UU(:,:,o),S(:,:,o)]=eigs(U(:,:,o),NC);end
        else
            for o=1:O;[UU(:,:,o),S(:,:,o)]=eigs(U(:,:,o),NC);end
        end             
        
        S=diagm(S);           
        [U,S]=parUnaFun({UU,S},@single);
        if gpu;[U,S]=parUnaFun({U,S},@gpuArray);end
    end
else
    X=matfun(@mtimes,X,matfun(@ctranspose,X));
%    S=sqrt(abs(eigm(X)),useParallel);
    S=sqrt(abs(eigm(X)));
    S=sort(S,1,'descend');
    if ~isempty(NC);S=S(1:NC,:,:);end
    S=permute(S,[2 1 3]);
end
S=resSub(S,3,NX(3:end));
U=resSub(U,3,NX(3:end));
if nargout>2;V=resSub(V,3,NX(3:end));end
