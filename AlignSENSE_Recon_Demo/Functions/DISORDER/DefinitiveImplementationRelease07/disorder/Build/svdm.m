function [S,U,V]=svdm(X,NC)

%SVDM   Performs accelerated multiple svds for wide and square matrices 
%based on a Schur decomposition of the normal matrix X*X'
%   [S,U,V]=SVDM(X,{NC})
%   * X is the input matrix
%   * {NC} indicates to return a given number of components
%   ** S are the singular values of X
%   ** U is the column space (range) of X
%   ** V is the row space (range) space of X
%

if nargin<2;NC=[];end

gpu=isa(X,'gpuArray');
[M,N,O]=size(X);
%assert(M<=N,'Only wide and square matrices are allowed but input size is %dx%d',M,N);

if nargout>1
    Xp=matfun(@ctranspose,X);
    U=matfun(@mtimes,X,Xp);
    U=(U+matfun(@ctranspose,U))/2;%Force the matrix to be Hermitian  
    if isempty(NC) || nargout>2
        U=gather(U);S=U;  %YB: Need to gather since shur.m does not support GPU functionality
        if O>=8
            parfor o=1:O;[U(:,:,o),S(:,:,o)]=schur(U(:,:,o));end
%             for o=1:O;[U(:,:,o),S(:,:,o)]=schur(U(:,:,o));end
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
        %U=bsxfun(@times,U,1./sqrt(dot(U,U,1)));
        if nargout>2;V=bsxfun(@times,matfun(@mtimes,Xp,U),1./D);end%PROBLEM WITH SINGLETONS IS THAT THIS MAY DEPART A LOT FROM ORTHONORMALITY!!
    else
        U=double(gather(U));
        UU=zeros([M NC O],'like',U);
        S=zeros([NC NC O],'like',U);
        if O>=8
            parfor o=1:O;[UU(:,:,o),S(:,:,o)]=eigs(U(:,:,o),NC);end
        else
            for o=1:O;[UU(:,:,o),S(:,:,o)]=eigs(U(:,:,o),NC);end
        end                
        S=diagm(S);           
        [U,S]=parUnaFun({UU,S},@single);
        if gpu;[U,S]=parUnaFun({U,S},@gpuArray);end

        %V=zeros([NC N O],'like',X);
        %if O>=8
        %    parfor o=1:O;[U(:,:,o),S(:,:,o),V(:,:,o)]=svds(X(:,:,o),NC);end
        %else
        %    for o=1:O;[U(:,:,o),S(:,:,o),V(:,:,o)]=svds(X(:,:,o),NC);end
        %end                
        %S=diagm(S);
        %[U,S,V]=parUnaFun({U,S,V},@single);
        %if gpu;[U,S,V]=parUnaFun({U,S,V},@gpuArray);end   
    end
else
    X=matfun(@mtimes,X,matfun(@ctranspose,X));
    S=sqrt(abs(eigm(X)));
    S=sort(S,1,'descend');
    if ~isempty(NC);S=S(1:NC,:,:);end
    S=permute(S,[2 1 3]);
end


