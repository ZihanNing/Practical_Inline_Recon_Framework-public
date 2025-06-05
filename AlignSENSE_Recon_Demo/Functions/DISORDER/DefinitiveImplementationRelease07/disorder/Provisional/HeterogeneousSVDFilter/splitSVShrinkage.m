function [X,amse]=splitSVShrinkage(Y,gamma,overlapFactor,Lm,Ln,Im,In,recomputeSVD)

%SPLITVSHRINKAGE  denoises a volume by splitting (or overlapping) over
%intensitiess in the column space
%   [AMSE,X]=SPLITSVSHRINKAGE(Y,LR,LN,IM,IN)
%   * Y is the matrix to be denoised
%   * {GAMMA} is the desired aspect ratio of the colummn space. It defaults
%   to 0 (global denoising)
%   * {OVERLAPFACTOR} is the overlap ratio. It defaults to 1
%   * {LM} are the variances of the rows (size Mx1)
%   * {LN} are the variances of the columns (size 1xN)
%   * {IM} is an indicator function for localizing along specific indexes
%   on the rows, (size MxI with I the number of subspaces). The columns
%   should add up to 1 (orthogonality would be another condition, but it is
%   not checked)
%   * {IN} is an indicator function for localizing along specific indexes
%   on the columns, (size JxN with J the number of subspaces). The rows
%   should add to 1 (orthogonality would be another condition, but it is
%   nnot checked)
%   * {RECOMPUTESVD} is a flag to recompute the SVD in for thresholding the
%   singular vectors, it defaults to 0
%   ** X is the estimated matrix
%   ** AMSE is the estimated asymptotic mean square error (per subspace)
%

[M,N]=size(Y);MN=[M N];

if nargin<2 || isempty(gamma);gamma=0;end
if nargin<3 || isempty(overlapFactor);overlapFactor=1;end
if nargin<4 || isempty(Lm);Lm=ones([M 1],'like',Y);end
if nargin<5 || isempty(Ln);Ln=ones([1 N],'like',Y);end
if nargin<6 || isempty(Im);Im=ones([M 1],'like',Y);end
if nargin<7 || isempty(In);In=ones([1 N],'like',Y);end
if nargin<8 || isempty(recomputeSVD);recomputeSVD=0;end
verbosity=0;

%RUN A GLOBAL SVSHRINKAGE (WHICH GIVES DECENT RESULTS VERY
%QUICKLY)
Xg=heterogeneousSVShrinkage(Y,Lm,Ln,Im,In,recomputeSVD);

%SORT THE COLUMNS ACCORDING TO THE ENERGY
%Xgn=bsxfun(@times,Xg,1./sqrt(sum(abs(Xg).^2,2)));%Normalize the columns
%F=sum(abs(Xgn).^2,1);
F=sum(abs(Xg).^2,1);
%F=sum(bsxfun(@times,V,S).^2,2)';
%F=sum(abs(V).^2,2)';

%figure
%histogram(10*log10(F),ceil(N/40));
%pause

[~,iS]=sort(F,2,'descend');
Y=indDim(Y,iS,2);
Ln=Ln(iS);
if any(In(:)~=0);In=indDim(In,iS,2);end

%INITIALIZE THE OUTPUT
W=zeros([1 N],'like',Y);%CUMMULATIVE WEIGHTS MATRIX
X=zeros([M N],'like',Y);%RESULTING MATRIX
amse=zeros([M N],'like',Y);%ERROR MATRIX

%DEFINE A CONVOLUTION MATRIX FOR SPLITTING
Ms=min(M/gamma,N);%SUBMATRIX SIZE
overlapFactor=max(overlapFactor,1);
if Ms==N;overlapFactor=1;end
vS=(1:Ms);%SUBMATRIX INDEXES
NB=floor(overlapFactor*N/Ms);%NUMBER OF BLOCKS
sh=Ms/overlapFactor;%SHIFT AMOUNT
indW=sh*(0:NB-1);%STARTING INDEXES
for b=1:NB    
    vSn=floor(vS+indW(b));vSn(vSn>N)=[];%For prevention although don't think it will happen
    if b==NB;vSn=[vSn vSn(end)+1:N];end%To process the last bits
    Yb=Y(:,vSn);
    Lnb=Ln(:,vSn);
    if any(In(:)~=0);Inb=In(:,vSn);else Inb=0;end
    [Xb,amseb]=heterogeneousSVShrinkage(Yb,Lm,Lnb,Im,Inb,recomputeSVD);
    X(:,vSn)=X(:,vSn)+Xb;
    amse(:,vSn)=amse(:,vSn)+amseb;
    W(vSn)=W(vSn)+1;
end
X=bsxfun(@times,X,1./W);
amse=bsxfun(@times,amse,1./W);

%REVERSE SORTING
X(:,iS)=X;
amse(:,iS)=amse;

