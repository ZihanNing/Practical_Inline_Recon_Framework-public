function x=whiten(x,A,diV,co,R,I,ish)

%WHITEN  applies noise whitening operations for approach inspired on [1] 
%W Leeb, "Matrix denoising for weighted loss functions and heterogeneous 
%signals", arXiv, 1902.0947v1, 2019
%   X=WHITEN(X,A,DI,CO,{R},{I})
%   * X is the array to be whitened
%   * A is a set of matrices to whiten. Cell with length given by the 
%   dimensions of the space, we assume they have been damped beforehand to
%   avoid numerical issues
%   * DIV is the direction: 1 for whiten, 0 for dewhiten, vector to apply
%   one and then another
%   * CO is the coordinate over which to apply whiten: 1 for rows, 2 for
%   columns
%   * {R} are the underlying dimensions of the space
%   * {I} are some indexes to be extracted
%   * {ISH} indicates whether we are dealing with Hermitian matrices, it
%   defaults to 0
%   * X is the whitened/dewhitened array
%

if nargin<6;I=[];end
if nargin<7 || isempty(ish);ish=0;end

dev=gpuDevice;
if ~isempty(A) || ~isempty(I)
    if co==2;x=x';end
    for di=diV
        %if di==1 && ~isempty(I);x=bsxfun(@times,x,I);end%This may be slightly quicker but perhaps riskier
        if di==1 && ~isempty(I);x(I,:)=x;x(~I,:)=0;end
        if ~isempty(A)
            [M,N]=size(x);
            if nargin<5 || isempty(R);R=M;end

            ND=length(R);%Number of dimensions
            NA=length(A);%Number of operations
        
            x=reshape(x,[R N]);
            if di==1;nV=1:NA;else nV=NA:-1:1;end
            for n=nV
                %if ~isempty(I);wait(dev);tic;end
                NDA=numDims(A{n});
                if numel(A{n})==M && NDA==ND+1;x=bsxfun(@times,reshape(A{n},[R 1]).^((-1)^di),x);
                elseif numel(A{n})==R(NDA);x=bsxfun(@times,A{n}.^((-1)^di),x);
                elseif numel(A{n})==R(NDA)^2
                    NNA=size(A{n});
                    A{n}=reshape(A{n},[R(NDA) R(NDA)]);
                    perm=1:ND+1;perm([1 NDA])=[NDA 1];
                    x=permute(x,perm);
                    if di==1
                        if ~ish;x=matfun(@mldivide,A{n},x);else x=matfun(@mtimes,A{n}',x);end
                    else                        
                        x=matfun(@mtimes,A{n},x);
                    end
                    x=permute(x,perm);
                    A{n}=reshape(A{n},NNA);
                else
                    error('Element %d of whitening could not be interpreted: %d',n);
                end
                %if ~isempty(I);wait(dev);toc;end%THIS COULD POTENTIALLY BE A BIT QUICKER IF USING THE FFT CODE
            end
            x=reshape(x,[M N]);  
        end
        %if di==0 && ~isempty(I);x=bsxfun(@times,x,I);end%This may be slightly quicker but perhaps riskier
        if di==0 && ~isempty(I);x=x(I,:);end      
    end
    if co==2;x=x';end
end