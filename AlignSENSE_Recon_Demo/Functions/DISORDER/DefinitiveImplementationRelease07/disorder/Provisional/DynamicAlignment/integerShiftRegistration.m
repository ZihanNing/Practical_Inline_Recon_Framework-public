function [x,T,y]=integerShiftRegistration(x,W,T,kmax,dk,lev,fracOrd)

%INTEGERSHIFTREGISTRATION   Per-volume groupwise rigid registration of dynamic studies
%   [X,T]=INTEGERSHIFTREGISTRATION(X,W,{T},{KMAX},{DK},{LEV},{FRACORD},{PARCOMP})
%   * X is the data to be registered
%   * W is a mask for the ROI to compute the metric
%   * {T} are the initial motion parameters. They default to all 0
%   * {KMAX} are the translation search radii. They default to [3 3]
%   * {DK} is the quantization step. It defaults to [2 1]
%   * {LEV} are the resolution levels for correction. They default to [2 2]
%   * {FRACORD} is the fractional finite difference metric order
%   * X is the registered data
%   * T are the estimated motion parameters
%   * Y is the average in time
%

NX=size(x);NX(end+1:8)=1;ND=max(numDims(x),3);

if nargin<3 || isempty(T);T=single(zeros([ones(1,ND-1) NX(ND) 6]));end
NXV=NX(1:3);
if nargin<4 || isempty(kmax);kmax=[3 3];end
if nargin<5 || isempty(dk);dk=[2 1];end
if nargin<6 || isempty(lev);lev=[2 2];end
if nargin<7 || isempty(fracOrd);fracOrd=0;end

gpu=isa(x,'gpuArray');
NT=size(T);ndT=ndims(T);
%ONLY INTEGER SHIFTS
T=dynInd(T,1:3,ndT,round(dynInd(T,1:3,ndT)));

L=length(lev);
assert(length(kmax)==L,'Number of levels of translation search radious (%d) do not match number of levels of multirresolution (%d)',length(kmax),L);
assert(length(dk)==L,'Number of levels of quantization (%d) do not match number of levels of multirresolution (%d)',length(dk),L);


BlSz=10;
for l=1:L    
    fprintf('Resolution level: %d\n',lev(l));            
    NXres=round(NXV/lev(l));
    
    %GRID GENERATION AND RESAMPLING
    [~,kGrid,rkGrid]=generateTransformGrids(NXV,gpu,NXres,[],1);
    [FT,FTH]=buildStandardDFTM(NXres,0,gpu);
    xRes=resampling(x,NXres);WRes=abs(resampling(W,NXres));    
    H=buildFilter(2*NXres,'tukeyIso',ones(1,3),gpu,0.1,1);
    xRes=filtering(xRes,H,1); 
    
    %CANDIDATE TRANSLATIONS
    NGD=ceil((2*kmax(l)+1)*ones(1,3));
    xx=generateGrid(NGD,gpu,NGD,ceil((NGD+1)/2));    
    r=xx{1}(1);r(1)=double(0);
    for n=1:3
        r=bsxfun(@plus,r,(xx{n}*dk(l)).^2);
    end   
    [rir,ir]=sort(r(:));
    rmax2=(dk(l)*kmax(l))^2;
    M=gather(find(rir>rmax2,1)-1);
    fprintf('Number of candidate translations: %d\n',M);
    irM=ir(1:M);
    irs=ind2subV(NGD,irM);
    for n=1:3;xx{n}=xx{n}(irs(:,n));xx{n}=xx{n}(:);end    
    xx=cat(2,xx{:});
    xx=repmat(xx,[1 1 2]);xx(:,:,2)=0;xx=reshape(xx,[M 6]);
    perm=1:ndT;perm([2 ndT])=[ndT 2];
    xx=permute(xx,perm);
    Tcand=bsxfun(@plus,T,xx);
    
    %ENERGY
    E=single(zeros([M NT(ndT-1)]));
    if gpu;E=gpuArray(E);end

    %FRACTIONAL DERIVATIVE WEIGHTS    
    if fracOrd~=0;GRes=buildFilter(NXres(1:2),'FractionalFiniteDiscreteIso',NX(1:2)./NXres(1:2),gpu,fracOrd);end%For fractional finite difference motion estimation
    
    %TEMPLATE
    y=single(zeros([NXres NX(4:ndT-2)]));
    if gpu;y=gpuArray(y);end   
    for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));%This is the ''reconstruction'' step
        etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1);
        %if s==1;y=sincRigidTransform(dynInd(xRes,vS,ndT-1),etInv,0,FT,FTH,0);else y=cat(ndT-1,y,sincRigidTransform(dynInd(xRes,vS,ndT-1),etInv,0,FT,FTH,0));end            
        y=y+sum(sincRigidTransform(dynInd(xRes,vS,ndT-1),etInv,0,FT,FTH,0),ndT-1)/NT(ndT-1);
    end
    
        
    permN=1:ndT-1;permN([2 ndT-1])=[ndT-1 2];
    %ENERGY COMPUTATION
    for m=1:M
        %fprintf('Processing %d of %d\n',m,M);
        for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));                        
            %et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Tcand,{m,vS},[1 ndT-1]),1,0,[],1);
            %xT=sincRigidTransform(y,et,1,FT,FTH,[],0);
            %WT=abs(sincRigidTransform(WRes,et,1,FT,FTH,[],0));
            %xT=xT-dynInd(xRes,vS,ndT-1);
            %if fracOrd~=0;xT=filtering(xT,GRes);end
            %xT=bsxfun(@times,xT,WT);%IT SHOULD BE SQRT WT, THIS IS A DIRTY TRICK TO MAKE IT WORK BETTER... ONCE TRACKING IMPLEMENTED PROBABLY NOT NECESSARY ANYMORE
            %E(m,vS)=permute(multDimSum(real(xT.*conj(xT)),1:ndT-2),permN);          
            
            et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Tcand,{m,vS},[1 ndT-1]),0,0,[],1);
            xT=sincRigidTransform(dynInd(xRes,vS,ndT-1),et,0,FT,FTH,0,0);
            xT=bsxfun(@minus,xT,y);
            if fracOrd~=0;xT=filtering(xT,GRes);end
            xT=bsxfun(@times,xT,WRes);            
            E(m,vS)=permute(multDimSum(real(xT.*conj(xT)),1:ndT-2),permN);                      
        end
    end    
    [~,iMin]=min(E,[],1);
    iMin=permute(iMin,permN);
    T=indDim(Tcand,iMin,1);
end

%SHIFTED DATA
[~,kGrid,rkGrid]=generateTransformGrids(NXV,gpu,[],[],1);
[FT,FTH]=buildStandardDFTM(NXV,0,gpu);
for s=1:BlSz(1):NT(ndT-1);vS=s:min(s+BlSz(1)-1,NT(ndT-1));
    etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1);
    x=dynInd(x,vS,ndT-1,sincRigidTransform(dynInd(x,vS,ndT-1),etInv,0,FT,FTH,0));
end
y=mean(x,ndT-1);
