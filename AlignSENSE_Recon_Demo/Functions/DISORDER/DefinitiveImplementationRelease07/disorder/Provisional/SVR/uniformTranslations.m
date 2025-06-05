function T=uniformTranslations(Ntr,trLim,T,bid,randtr)

%UNIFORMTRANSLATIONS generates uniform translations evenly distributed 
%within a given interval
%   T=UNIFORMROTATIONS(NTR,TRLIM,{T},{BID})
%   * NRO is the number of translations to be generated (per dimension)
%   * TRLIM is the translation limit
%   * {T} is a previous transform to which to add translations
%   * {BID} indicates that the translations are for a 2D domain, it 
%   defaults to 0
%   * {RANDTR} serves to generate a random translation origin, it defaults
%   to 0
%   ** T is the new transform arranged as 3xN
%

if nargin<4 || isempty(bid);bid=0;end
if nargin<5 || isempty(randtr);randtr=0;end

if Ntr==0
    x=[0 0 0];%Zero translation only
else
    if mod(Ntr,2)==0;Ntr=Ntr+1;end%To force using the zero translation    
    if bid;ND=2;else ND=3;end
    y=generateGrid(Ntr*ones(1,ND),0,2*trLim*ones(1,ND),ceil((Ntr*ones(1,ND)+1)/2));   
    if randtr
       for n=1:ND;y{n}=y{n}+2*trLim*(rand(1)-0.5)/Ntr;end
    end
    if bid;[x{1},x{2}]=ndgrid(y{1}(:)',y{2}(:)');else [x{1},x{2},x{3}]=ndgrid(y{1}(:)',y{2}(:)',y{3}(:)');end
    for n=1:ND;x{n}=x{n}(:);end    
    y=cat(2,x{:});
    y(:,ND+1:3)=0;
end

if nargin>=3 && ~isempty(T)
    nDT=numDims(T);    
    perm=1:nDT+1;perm(nDT:nDT+1)=[nDT+1 nDT];%We extend by one dimension
    T=permute(T,perm);
    Tr=dynInd(T,1:3,nDT+1);    
    NY=size(y);
    y=reshape(y,[ones(1,nDT-1) NY]);   
    Tr=bsxfun(@plus,Tr,y);
    rep=ones(1,nDT+1);rep(nDT)=NY(1);
    T=repmat(T,rep);
    T=dynInd(T,1:3,nDT+1,Tr);
else
    T=y;
end
