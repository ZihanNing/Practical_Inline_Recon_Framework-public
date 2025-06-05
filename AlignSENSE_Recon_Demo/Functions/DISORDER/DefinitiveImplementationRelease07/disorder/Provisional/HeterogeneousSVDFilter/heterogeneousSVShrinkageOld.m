function [X,amse,U,S,V]=heterogeneousSVShrinkageOld(Y,Lm,Ln,Im,In,Hm,Hn,M,N)

%HETEROGENEOUSSVSHRINKAGE  applies a shrinkage of the singular values of Y 
%according to [1] W Leeb, "Matrix denoising for weighted loss functions 
%and heterogeneous signals", arXiv, 1902.0947v1, 2019
%   [AMSE,X]=HETEROGENEOUSSVSHRINKAGE(Y,LR,LN,IM,IN)
%   * Y is the matrix to be denoised
%   * {LM} are noise generation matrices for the rows (size Mx1 or MxW)
%   * {LN} are noise generation matrices for the columns (size 1xN or WxN)
%   * {IM} is an indicator function for localizing along specific indexes
%   on the rows, (size MxI with I the number of subspaces). The columns
%   should add up to 1 (orthogonality would be another condition, but it is
%   not checked)
%   * {IN} is an indicator function for localizing along specific indexes
%   on the columns, (size JxN with J the number of subspaces). The rows
%   should add to 1 (orthogonality would be another condition, but it is
%   nnot checked)
%   * {HM} is a (orthogonal) decomposition for the rows
%   * {HN} is a (orthogonal) decomposition for the columns
%   * {MEFF} is the effective dimension of the row space, in case there is
%   masking that constraints this space
%   * {NEFF} is the effective dimension of the column space, in case there
%   is masking that constraints this space
%   ** X is the estimated matrix
%   ** AMSE is the estimated asymptotic mean square error (per subspace)
%   ** U are the estimated left singular vectors
%   ** S are the estimated singular values
%   ** V are the estimated right singular vectors
%

[M,N]=size(Y);
if nargin<2 || isempty(Lm);Lm=ones([M 1],'like',Y);end
if nargin<3 || isempty(Ln);Ln=ones([1 N],'like',Y);end
if nargin<4 || isempty(Im);Im=ones([M 1],'like',Y);end
if nargin<5 || isempty(In);In=ones([1 N],'like',Y);end
if nargin<6;Hm=[];end
if nargin<7;Hn=[];end
verbosity=1;

%WE FLIP THE ENTRIES
flEn=0;
if M>N
    flEn=1;
    Y=Y';
    aux=Ln;Ln=Lm';Lm=aux';
    aux=In;In=Im';Im=aux';
    aux=Hn;Hn=Hm';Hm=aux';
    [M,N]=size(Y);
end
MN=[M,N];

gpu=isa(Y,'gpuArray');
[Y,Im,In]=parUnaFun({Y,Im,In},@gather);
sing=isa(Y,'single');
[Y,Im,In]=parUnaFun({Y,Im,In},@double);
if ~isa(Hm,'function_handle');Hm=double(gather(Hm));end
if ~isa(Hn,'function_handle');Hn=double(gather(Hn));end
if ~isa(Lm,'function_handle');Lm=double(gather(Lm));end
if ~isa(Ln,'function_handle');Ln=double(gather(Ln));end

%NOISE WHITENING
if ~isa(Lm,'function_handle');if ~isvector(Lm);Y=mldivide(Lm,Y);else Y=bsxfun(@times,1./Lm,Y);end
else Y=Lm(Y,0);
end
if ~isa(Ln,'function_handle');if ~isvector(Ln);Y=mrdivide(Y,Ln);else Y=bsxfun(@times,1./Ln,Y);end
else Y=Ln(Y,0);
end

%SVD
comp=~isreal(Y);
Ynorm=sqrt(N*(1+comp));
[S,U,V]=svdm(Y/Ynorm);%S column vector

%WEIGHTED VECTORS
if isvector(Lm);Uw=bsxfun(@times,Lm,U);else Uw=Lm*U;end
if isvector(Lm);WLm=diag(Lm);else WLm=Lm;end

if isvector(Ln);Vw=bsxfun(@times,Ln',V);else Vw=Ln'*V;end
if isvector(Ln);WLn=diag(Ln);else WLn=Ln';end


%SHRINKAGE OVER THE DIFFERENT SUBSPACES
I=size(Im,2);J=size(In,1);

X=zeros([M N],'like',Y);
amse=zeros([M N],'like',Y);
for i=1:I%%%STILL PROBLEMATIC PROBLEM WHEN ACTIVATING ROTATION!!!
    Imi=Im(:,i)==1;
    if any(Imi)
        if ~isempty(Hm);Uiw=Hm(Imi,:)'*(Hm(Imi,:)*Uw);WLmi=Hm(Imi,:)'*(Hm(Imi,:)*WLm);
        else Uiw=Uw(Imi,:);WLmi=WLm(Imi,:);
        end
        for j=1:J
            Inj=In(j,:)==1;    
            if any(Inj)
                if ~isempty(Hn);Vjw=Hn(Inj,:)'*(Hn(Inj,:)*Vw);WLnj=Hn(Inj,:)'*(Hn(Inj,:)*WLn);
                else Vjw=Vw(Inj,:);WLnj=WLn(Inj,:);
                end

                %NORMALIZED SQUARED TRACES
                mum2=sum(diag(WLmi'*WLmi))/M;
                mun2=sum(diag(WLnj'*WLnj))/N;
                %mum2=sum(abs(WLmi(:)).^2)/M;
                %mun2=sum(abs(WLnj(:)).^2)/N;

                if (isvector(WLmi) && all(WLmi==1)) && (isvector(WLnj) && all(WLnj==1));ca=1;%Unweighted case
                elseif (isvector(WLmi) && all(WLmi==1)) || (isvector(WLnj) && all(WLnj==1));ca=2;%Asymptotically orthogonal with respect to the weights, either one sided or random signal vectors
                else ca=3;
                end

                %SINGULAR VALUE ESTIMATION
                [Sij,amseij,Rij]=weightedShrinkage(S,Uiw,Vjw,M/N,mum2,mun2,ca);               

                %SIGNAL SYNTHESIS
                enRat=(sum(Imi)*sum(Inj))/Ynorm^2;
                if isempty(Hm) && isempty(Hn)
                    X(Imi,Inj)=bsxfun(@times,Uiw,Ynorm*Sij)*Vjw';                                            
                    amse(Imi,Inj)=amseij/enRat;                
                elseif ~isempty(Hm) && isempty(Hn)
                    X(:,Inj)=X(:,Inj)+bsxfun(@times,Uiw,Ynorm*Sij)*Vjw';                
                    amse(:,Inj)=amse(:,Inj)+repmat(sum(abs(Hm(:,Imi)).^2,2)*amseij,[1 sum(Inj)])/enRat;
                elseif ~isempty(Hn) && isempty(Hm)
                    X(Imi,:)=X(Imi,:)+bsxfun(@times,Uiw,Ynorm*Sij)*Vjw';
                    amse(Imi,:)=amse(Imi,:)+repmat(sum(amseij*(abs(Hn(:,Inj)).^2)',1),[sum(Imi) 1])/enRat;
                else
                    X=X+bsxfun(@times,Uiw,Ynorm*Sij)*Vjw';
                    amse=amse+amseij/enRat;
                end
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
