function T=randomRigidTransforms(Ntr,trStd,rotStd,T,bid)

%RANDOMRIGIDTRANSFORMS generates random rigid transforms uniformly
%distributed in the spheres of rotations and translations and Gaussian
%distributed translations and rotations
%   T=RANDOMRIGIDTRANSFORMS(NTR,TRSTD,ROTSTD,{T},{BID})
%   * NTR is the number of transforms to be generated
%   * TRSTD is the translation std
%   * ROTSTD is the rotation std (in degrees)
%   * {T} is a previous transform to which to add transformations
%   * {BID} indicates that the transformations are for a 2D domain, it 
%   defaults to 0
%   ** T is the new transform arranged as 6xN
%

if nargin<5 || isempty(bid);bid=0;end

%RANDOM ORIENTATION
if ~bid
    x=randn(Ntr,3);
    x=bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));
else
    x=zeros(Ntr,3);
    x(:,3)=1;
end
%RANDOM ROTATION
x(:,4)=rotStd*randn(Ntr,1);%In degrees
x(x(:,4)<0,4)=x(x(:,4)<0,4)+360;

%RANDOM TRANSLATION
x(:,5:7)=trStd*randn(Ntr,3);
if bid;x(:,7)=0;end

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
    M=SpinCalc('EVtoDCM',x(:,1:4),1e-3,0);
    M=permute(M,perm);
    HT=emtimes(M,HT);HT(end+1:4)=1;%First the old rotations, then the new ones
    NHT=size(HT);NHT(end+1:4)=1;
    HT=reshape(HT,[NHT(1:2) prod(NHT(3:4))]);
    Tr=SpinCalc('DCMtoEA321',permute(HT,[2 1 3]),1e-3,0);
    Tr=-convertRotation(Tr,'deg','rad');
    Tr=reshape(Tr,[NT(1:nDT-1) NHT(4) 3]);
    
    Tt=reshape(x(:,5:7),[ones(1,nDT-1) Ntr 3]);
    Tt=bsxfun(@plus,Tt,dynInd(T,1:3,nDT+1));
    
    rep=ones(1,nDT+1);rep(nDT)=NHT(4);  
    T=repmat(T,rep);
    T=dynInd(T,4:6,nDT+1,Tr);
    T=dynInd(T,1:3,nDT+1,Tt);
else  
    T=SpinCalc('EVtoEA321',x(:,1:4),1e-3,0);
    T=convertRotation(T,'deg','rad');
    T=cat(2,x(:,5:7),T);
end
