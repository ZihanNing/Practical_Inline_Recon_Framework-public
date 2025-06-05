function T=uniformRotations(Nor,Nro,roLim,T,bid,randor)

%UNIFORMROTATIONS generates uniform rotations evenly distributed in a subset
%of SO(3)
%   T=UNIFORMROTATIONS(NOR,NRO,ROLIM,{T},{BID})
%   * NOR are a number of orientations uniformly covering the sphere using
%   the Vogel's method
%   * NRO is the number of rotations around previous orientations to be
%   generated
%   * ROLIM is the rotation limit (in degrees)
%   * {T} is a previous transform to which to add rotations
%   * {BID} indicates that the rotations are for a 2D domain, it defaults
%   to 0
%   * {RANDOR} serves to generate a random orientation origin, it defaults
%   to 0
%   ** T is the new transform or the Euler rotation angles arranged as 3xN
%

if nargin<5 || isempty(bid);bid=0;end
if nargin<6 || isempty(randor);randor=0;end

if Nor==0 && Nro==0
    x=[0 0 1 0];%Zero rotation only
else
    if bid==1
        Nro=round((Nor*Nro).^(1/3));
        Nor=1;
    end        
    x=vogelSphereGrid(Nor,1);
    if mod(Nro,2)==0;Nro=Nro+1;end%To force using the zero rotation (and centered coordinates), note it was 0 before    
    y=generateGrid(Nro,0,2*roLim,ceil((Nro+1)/2));
    y=y{1}(:)';    
    if Nro>1;y(y==0)=[];end
    y(y<0)=y(y<0)+360;
    y=repmat(y,[Nor 1]);
    if randor;y=bsxfun(@plus,y,2*roLim*(rand([Nor 1])-0.5)/Nro);end
    if Nro>1;x=repmat(x,[Nro-1 1]);else x=repmat(x,[Nro 1]);end
    x=cat(2,x,y(:));
    if Nro>1;x(end+1,:)=0;x(end,3)=1;end%To add the zero rotation    
end

if nargin>=4 && ~isempty(T)
    nDT=numDims(T);    
    perm=1:nDT+1;perm(nDT:nDT+1)=[nDT+1 nDT];%We extend by one dimension
    T=permute(T,perm);
    NT=size(T);
    Tr=dynInd(T,4:6,nDT+1);    
    Tr=reshape(Tr,[prod(NT(1:nDT)) 3]);   
    Tr=-convertRotation(Tr,'rad','deg');
    HT=permute(SpinCalc('EA321toDCM',Tr,1e-3,0),[2 1 3]);
    perm=1:4;perm(3:4)=[4 3];    
    M=SpinCalc('EVtoDCM',x,1e-3,0);
    M=permute(M,perm);
    HT=emtimes(M,HT);HT(end+1:4)=1;%First the old rotations, then the new ones
    NHT=size(HT);NHT(end+1:4)=1;
    HT=reshape(HT,[NHT(1:2) prod(NHT(3:4))]);
    Tr=SpinCalc('DCMtoEA321',permute(HT,[2 1 3]),1e-3,0);
    Tr=-convertRotation(Tr,'deg','rad');
    Tr=reshape(Tr,[NT(1:nDT-1) NHT(4) 3]);
    rep=ones(1,nDT+1);rep(nDT)=NHT(4);
    T=repmat(T,rep);
    T=dynInd(T,4:6,nDT+1,Tr);
else  
    T=SpinCalc('EVtoEA321',x,1e-3,0);
    T=convertRotation(T,'deg','rad');
end
