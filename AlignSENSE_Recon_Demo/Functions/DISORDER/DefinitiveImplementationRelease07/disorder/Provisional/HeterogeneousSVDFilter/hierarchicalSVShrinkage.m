function [X,amse]=hierarchicalSVShrinkage(Y,minGammaSplit,minGammaLocalize,Lm,Ln,Im,In,Hm,Hn,xx)

%HIERARCHICALSVSHRINKAGE  denoises a volume by splitting (or overlapping) over
%intensitiess in the column space
%   [X,AMSE]=HIERARCHICALSVSHRINKAGE(Y,LR,LN,IM,IN)
%   * Y is the matrix to be denoised
%   * {MINGAMMASPLIT) defines the minimum aspect ratio for matrix
%   splitting, it defaults to 1
%   * {MINGAMMALOCALIZE) defines the minimum aspect ratio for matrix
%   localization, it defaults to 0.2
%   * {OVERLAPFACTOR} is the overlap ratio. It defaults to 1
%   * {LM} are the variances of the rows (size Mx1)
%   * {LN} are the variances of the columns (size 1xN)
%   * {IM} is an indicator function for localizing along specific indexes
%   on the rows, (size MxI with I the number of subspaces). The columns
%   should add up to 1
%   * {IN} is an indicator function for localizing along specific indexes
%   on the columns, (size JxN with J the number of subspaces). The rows
%   should add up to 1
%   * {HM} is an orthogonal decomposition for the rows
%   * {HN} is an orthogonal decomposition for the columns
%   * {RECOMPUTESVD} is a flag to recompute the SVD in for thresholding the
%   singular vectors, it defaults to 0
%   ** X is the estimated matrix
%   ** AMSE is the estimated asymptotic mean square error (per subspace)
%

[M,N]=size(Y);MN=[M N];

if nargin<2 || isempty(minGammaSplit);minGammaSplit=1;end
if nargin<3 || isempty(minGammaLocalize);minGammaLocalize=0.2;end
if nargin<4 || isempty(Lm);Lm=ones([M 1],'like',Y);end
if nargin<5 || isempty(Ln);Ln=ones([1 N],'like',Y);end
if nargin<6 || isempty(Im);Im=ones([M 1],'like',Y);end
if nargin<7 || isempty(In);In=ones([1 N],'like',Y);end
if nargin<8;Hm=[];end
if nargin<9;Hn=[];end
if nargin<10 || isempty(xx);xx=zeros(3,N);end
verbosity=0;

typFeat='NormEnerg';
%typFeat='Alpha';
%typFeat='AlphaW';

%RUN A GLOBAL SVSHRINKAGE (WHICH GIVES DECENT RESULTS VERY
%QUICKLY)
Inaux=ones([1 N],'like',Y);
[X,amse,~,S{1},V]=heterogeneousSVShrinkage(Y,Lm,Ln,Im,Inaux,[],[],minGammaLocalize);

C{1}=1:N;%Initial cluster
Cconv=0;
nh=0;
NMmin=min(MN)*minGammaSplit;
if ~isempty(xx)
    xx=bsxfun(@minus,xx,mean(xx,2));
    xx=bsxfun(@times,xx,1./(std(xx,0,2)+eps));
end
dimSub=2;%Dimension of subdivision
while 1
    if strcmp(typFeat,'NormEnerg')
        [CParent,CconvParent,C,Cconv,Parent,Childs,S]=subdivideClusters(C,Cconv,X,xx,dimSub,NMmin,typFeat,S);    
    elseif strcmp(typFeat,'Alpha') || strcmp(typFeat,'AlphaW')
        [CParent,CconvParent,C,Cconv,Parent,Childs,S]=subdivideClusters(C,Cconv,V',xx,dimSub,NMmin,typFeat,S);
    end
    nh=nh+1;
    if all(Cconv)
        fprintf('Number of levels: %d\n',nh);
        break;
    end
    NC=length(C);
    Xcand=cell(1,NC);amsecand=cell(1,NC);Vcand=cell(1,NC);Scand=cell(1,NC);
    for n=1:length(C)
        if ~Cconv(n);[Xcand{n},amsecand{n},~,Scand{n},Vcand{n}]=heterogeneousSVShrinkage(Y(:,C{n}),Lm,Ln(:,C{n}),Im,Inaux(:,C{n}),[],[],minGammaLocalize);end
    end
    NC=length(CParent);
    for n=1:NC
        if ~CconvParent(n)
            if multDimSum([amsecand{Childs{n}(1)} amsecand{Childs{n}(2)}],1:2)<multDimSum(amse(:,CParent{n}),1:2)
                for m=1:2
                    amse(:,C{Childs{n}(m)})=amsecand{Childs{n}(m)};
                    X(:,C{Childs{n}(m)})=Xcand{Childs{n}(m)};
                    V(C{Childs{n}(m)},:)=Vcand{Childs{n}(m)};
                    S{Childs{n}(m)}=Scand{Childs{n}(m)};
                end
            else
                Cconv(Childs{n})=1;
            end
        end
    end
end

if all(In==0)
    NC=length(C);
    for n=1:NC;[X(:,C{n}),amse(:,C{n})]=heterogeneousSVShrinkage(Y(:,C{n}),Lm,Ln(:,C{n}),Im,In,Hm,Hn(C{n},C{n}),minGammaLocalize);end
end
