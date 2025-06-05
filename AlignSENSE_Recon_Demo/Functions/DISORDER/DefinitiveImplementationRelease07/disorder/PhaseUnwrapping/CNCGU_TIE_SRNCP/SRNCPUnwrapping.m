function [x,B]=SRNCPUnwrapping(x,voxSiz,weightFunc)

%SRNCPUNWRAPPING calls a 3D version of [1] M Arevallilo Herraez et al. 
%"Fast two-dimensional phase-unwrapping algorithm based on sorting by 
%reliability following a noncontinuous path," Apl Opt, 41(35):7437-44
%   [X,B]=SRNCPUNWRAPPING(X,{VOXSIZ},{WEIGHTFUNC})
%   * X is the complex data used to estimate the field
%   * {VOXSIZ} are the voxel sizes of the data
%   * {WEIGHTFUNC} weighting function
%   * X is the estimation of the unwrapped field
%   * B is the distance between the unwrapped and the wrapped field
%

if nargin<2 || isempty(voxSiz)
    ND=numDims(x);
    voxSiz=ones(1,ND);
elseif numDims(voxSiz)==0
    ND=numDims(x);
else
    ND=length(voxSiz);
end
if nargin<3 || isempty(weightFunc);weightFunc='Magnitude';end

gpu=isa(x,'gpuArray');

%UNWRAPPING PARAMETERS
w=voxSiz;w(end+1:ND)=1;
%fprintf('Used weights:%s\n',sprintf(' %.2f',w));

assert(ND==3,'Method only implemented for 3 dimensions');
N=size(x);N=N(1:ND);

NO=size(x);
x=resSub(x,ND+1:numDims(x));
NL=size(x,ND+1);
B=angle(x);x=abs(x);
for l=1:NL
    l
    xl=dynInd(x,l,ND+1);
    Bl=dynInd(B,l,ND+1);

    %COMPUTE EDGES AND SORT BY RELIABILITY
    EV=[];ER=[];
    for n=1:ND
        vn=zeros(1,ND);vn(n)=1;
        NN=N-vn;NNP=prod(NN);
        vV=ones(1,'like',x);
        if ~isempty(strfind(weightFunc,'Magnitude'));vV=vV.*(dynInd(xl,1:N(n)-1,n).*dynInd(xl,2:N(n),n));end
        %if ~isempty(strfind(weightFunc,'Gradient'));vV=vV./abs(wrapToPi(dynInd(B,1:N(n)-1,n)-dynInd(B,2:N(n),n)));end
        if ~isempty(strfind(weightFunc,'Gradient'));vV=vV.*abs(cos(wrapToPi(dynInd(Bl,1:N(n)-1,n)-dynInd(Bl,2:N(n),n))/2));end

        ER=cat(1,ER,vV(:));%Reliability of edge
        vV=(1:NNP)';
        eva=ind2subV(NN,vV);evb=bsxfun(@plus,eva,vn);
        eva=sub2indV(N,eva);evb=sub2indV(N,evb);
        EV=cat(1,EV,cat(2,eva,evb));%Voxel of edge 
    end;eva=[];evb=[];
    if gpu;EV=gpuArray(EV);end
    [~,is]=sort(ER,'descend');ER=[];%Sorting index
    EV=EV(is,:);%Sort by reliability
    NP=prod(N);
    G=(1:NP)';%Groups
    NG=G;NG(:)=1;


    %UNWRAP
    Bl=Bl(:);n=1;nn=1;
    while nn<NP
        indE=EV(n,:);
        g=G(indE);
        if g(1)~=g(2)
            Ges=NG(g);
            b=diff(Bl(indE));
            b=b-wrapToPi(b);
            if Ges(1)>Ges(2);b=-b;else g=flip(g);end
            iG=find(G==g(2));%%%%%%PROBABLY SOME COMPLEXITY HERE AS WELL
            if abs(b)>pi;Bl(iG)=Bl(iG)+b;end
            G(iG)=g(1);%%%%%THIS IS WHERE THE COMPUTATIONAL COMPLEXITY LIES
            NG(g(1))=sum(Ges);
            nn=nn+1;     
        end
        n=n+1;
    end
    Bl=reshape(Bl,N);
    x=dynInd(x,l,ND+1,Bl);
end

x=reshape(x,NO);B=reshape(B,NO);
if nargout>=2;B=angle(exp(1i*(x-B)));end
