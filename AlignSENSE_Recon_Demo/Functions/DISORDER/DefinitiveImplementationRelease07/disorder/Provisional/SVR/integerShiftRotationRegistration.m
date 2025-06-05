function [x,T,y,Emin]=integerShiftRotationRegistration(x,W,T,kmax,dk,phmax,Nrot,Nori,lev,fracOrd,voxsiz,BlSz,useBest,metric,center)

%INTEGERSHIFTROTATIONREGISTRATION   Per-volume groupwise rigid registration of dynamic studies
%   [X,T]=INTEGERSHIFTREGISTRATION(X,W,{T},{KMAX},{DK},{PHMAX},{NROT},{NORI},{LEV},{FRACORD},{PARCOMP})
%   * X is the data to be registered
%   * W is a mask for the ROI to compute the metric
%   * {T} are the initial motion parameters. They default to all 0
%   * {KMAX} are the translation search radii. They default to [3 3]
%   * {DK} is the quantization step. It defaults to [2 1]
%   * {PHMAX} is the range of rotation angles. It defaults to [40 20]
%   * {NROT} are the number of rotation angles. It defaults to [10 10]
%   * {NORI} are the number of orientations to search. They default to 
%   [10 10]
%   * {LEV} are the resolution levels for correction. They default to [2 2]
%   * {FRACORD} is the fractional finite difference metric order
%   * {VOXSIZ} is the voxel size
%   * {BLSZ} is the block size
%   * {USEBEST} is a flag to use the volume that best approximates the 
%   average over volumes (1) or the average volume (0, default)
%   * {METRIC} indicates the metric. 'MS': mean squares (default). 'NC':
%   normalized correlation
%   * {CENTER} indicates the center of rotation
%   * X is the registered data
%   * T are the estimated motion parameters
%   * Y is the average in time
%   * Emin is the estimated energy
%

NX=size(x);NX(end+1:8)=1;ND=max(numDims(x),4);

if nargin<3 || isempty(T);T=single(zeros([ones(1,ND-1) NX(ND) 6]));end
NXV=NX(1:3);
if nargin<4 || isempty(kmax);kmax=[3 3];end
if nargin<5 || isempty(dk);dk=[2 1];end
if nargin<6 || isempty(phmax);phmax=[40 40];end
if nargin<7 || isempty(Nrot);Nrot=[10 10];end
if nargin<8 || isempty(Nori);Nori=[10 10];end
if nargin<9 || isempty(lev);lev=[2 2];end
if nargin<10 || isempty(fracOrd);fracOrd=0;end
if nargin<11 || isempty(voxsiz);voxsiz=ones(1,3);end
if nargin<12 || isempty(BlSz);BlSz=[10 10];end
if nargin<13 || isempty(useBest);useBest=0;end
if nargin<14 || isempty(metric);metric='MS';end
if nargin<15 || isempty(center);center=(NXV+1)/2;end
kmax=max(kmax,0.001);

gpu=isa(x,'gpuArray');
NT=size(T);ndT=ndims(T);

%%ONLY INTEGER SHIFTS
%T=dynInd(T,1:3,ndT,round(dynInd(T,1:3,ndT)));

L=length(lev);
assert(length(kmax)==L,'Number of levels of translation search radious (%d) do not match number of levels of multirresolution (%d)',length(kmax),L);
assert(length(dk)==L,'Number of levels of quantization (%d) do not match number of levels of multirresolution (%d)',length(dk),L);
assert(length(phmax)==L,'Number of levels of rotation search angle (%d) do not match number of levels of multirresolution (%d)',length(phmax),L);
assert(length(Nrot)==L,'Number of rotations (%d) do not match number of levels of multirresolution (%d)',length(Nrot),L);
assert(length(Nori)==L,'Number of orientations (%d) do not match number of levels of multirresolution (%d)',length(Nori),L);

centShift=center-(NXV/2+1);
for l=1:L    
    centShiftRes=centShift/lev(l);
    fprintf('Resolution level: %d\n',lev(l));            
    NXres=round(NXV/lev(l));
    
    %GRID GENERATION AND RESAMPLING
    [~,kGrid,rkGrid,~,cGrid]=generateTransformGrids(NXV.*voxsiz,gpu,NXres,NXres/2+1+centShiftRes,1);
    [FT,FTH]=buildStandardDFTM(NXres,0,gpu);
    xRes=resampling(x,NXres,0,2*ones(1,3));WRes=abs(resampling(W,NXres,0,2*ones(1,3)));%Using mirror boundary conditions has had a strong impact here
    H=buildFilter(2*NXres,'tukeyIso',ones(1,3),gpu,0.1,1);
    xRes=filtering(xRes,H,1);
    
    %CANDIDATE TRANSLATIONS
    NGD=ceil((2*kmax(l)+1)*ones(1,3));
    xx=generateGrid(NGD,gpu,NGD,ceil((NGD+1)/2));
    r=xx{1}(1);r(1)=double(0);
    for n=1:3
        xx{n}=xx{n}*dk(l);
        r=bsxfun(@plus,r,xx{n}.^2);
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
    Tcand=gather(Tcand);
    
    %CANDIDATE ROTATIONS
    Tcand=uniformRotations(Nori(l),Nrot(l),phmax(l),Tcand);
    perm=1:ND+2;perm([2 ND+1 ND+2])=[ND+1 ND+2 2];
    Tcand=permute(Tcand,perm);
    R=size(Tcand,2);
    fprintf('Number of candidate rotations: %d\n',R);
    NTcand=size(Tcand);
    Tcand=reshape(Tcand,[prod(NTcand(1:2)) 1 NTcand(3:end)]);
    M=size(Tcand,1);
    fprintf('Number of candidate motion states: %d\n',M);

    %FRACTIONAL DERIVATIVE WEIGHTS    
    if fracOrd~=0;GRes=buildFilter(NXres(1:2),'FractionalFiniteDiscreteIso',NX(1:2)./NXres(1:2),gpu,fracOrd);end%For fractional finite difference motion estimation
    
    %TEMPLATE
    yF=zeros([NXres NX(4:ndT-1)],'like',xRes);   
    for s=1:BlSz(2):NT(ndT-1);vS=s:min(s+BlSz(2)-1,NT(ndT-1));%This is the ''reconstruction'' step
        etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1,cGrid);        
        yF=real(sincRigidTransform(dynInd(xRes,vS,ndT-1),etInv,0,FT,FTH,0));
    end
    y=mean(yF,ndT-1);
    if strcmp(metric,'MS')
        yD=bsxfun(@minus,yF,y);
        if fracOrd~=0;yD=filtering(yD,GRes);end
        yD=bsxfun(@times,yD,WRes);
        yD=multDimSum(real(yD.*conj(yD)),1:ndT-2);
    elseif strcmp(metric,'NC')
        %yaux=bsxfun(@times,y,WRes);   
        %yaux=yaux-mean(yaux(:));
        %yauxstd=sqrt(mean(real(yaux(:).*conj(yaux(:)))));
        %NYF=size(yF);NYF(end+1:ndT-1)=1;
        %dM=1:ndT-2;
        %yFaux=bsxfun(@times,yF,WRes);
        %yFaux=bsxfun(@minus,yFaux,mean(resPop(yFaux,dM,prod(NYF(dM)),1),1));
        %num=-mean(resPop(real(bsxfun(@times,yFaux,conj(yaux))),dM,prod(NYF(dM)),1),1);
        %yFaux=resPop(yFaux,dM,prod(NYF(dM)),1);
        %yFauxstd=sqrt(mean(real(yFaux.*conj(yFaux)),1));
        %den=bsxfun(@times,yauxstd,yFauxstd);
        %yD=num./den;
        
        yaux=y(WRes>0.5);
        yauxstd=std(yaux,0,1);   
        yaux=bsxfun(@minus,yaux,mean(yaux,1));
        NYF=size(yF);NYF(end+1:5)=1;
        yFaux=reshape(yF,[prod(NYF(1:3)) NYF(4:5)]);
        yFaux=yFaux(WRes>0.5,:,:);
        yFauxstd=std(yFaux,0,1);
        yFaux=bsxfun(@minus,yFaux,mean(yFaux,1));
        num=-mean(bsxfun(@times,yFaux,yaux),1);
        den=bsxfun(@times,yauxstd,yFauxstd);
        yD=num./den;
    end        
    %if R>1
    %    figure
    %    plot(yD(:))
    %    pause
    %end
    if useBest==2 || useBest==4;[~,iYD]=max(yD);elseif useBest==1;[~,iYD]=min(yD);elseif useBest==3;y=median(yF,ndT-1);end
    if ~ismember(useBest,[0 3])
        fprintf('Template stack: %d\n',iYD);
        y=dynInd(yF,iYD,ndT-1);
        if useBest==4;y=median(cat(ndT-1,yF,y),ndT-1);end                
    end
    yF=[];
    
    if strcmp(metric,'NC')        
        %yaux=bsxfun(@times,y,WRes);
        %yaux=yaux-mean(yaux(:));
        %yauxstd=sqrt(mean(real(yaux(:).*conj(yaux(:)))));        
        yaux=y(WRes>0.5);
        yauxstd=std(yaux,0,1);
        yaux=bsxfun(@minus,yaux,mean(yaux,1));
    end
    
    perm=1:ndT+1;perm([1 ndT-1:ndT+1])=[ndT+1 ndT-1 1 ndT];
    Tcandprev=Tcand;
    Tcand=permute(Tcand,perm);
    
    %ENERGY
    E=single(zeros([M NT(ndT-1)]));
    if gpu;E=gpuArray(E);end
        
    permN=1:ndT;permN([1:2 ndT-1:ndT])=[ndT ndT-1 1:2];
    %ENERGY COMPUTATION
    for m=1:BlSz(1):M;vM=m:min(m+BlSz(1)-1,M);
        %fprintf('Processing %d of %d\n',m,M);
        for s=1:BlSz(2):NT(ndT-1);vS=s:min(s+BlSz(2)-1,NT(ndT-1));
            et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Tcand,{vS,vM},ndT-1:ndT),0,0,[],1,cGrid);
            if Nrot(l)==0 && Nori(l)==0;xT=real(sincRigidTransform(dynInd(xRes,vS,ndT-1),et,0,FT,FTH,0,0));
            else xT=real(sincRigidTransform(dynInd(xRes,vS,ndT-1),et,0,FT,FTH,0));
            end
            if strcmp(metric,'MS')
                xT=bsxfun(@minus,xT,y);
                if fracOrd~=0;xT=filtering(xT,GRes);end
                xT=bsxfun(@times,xT,WRes);
                E(vM,vS)=permute(multDimSum(real(xT.*conj(xT)),1:ndT-2),permN);
            elseif strcmp(metric,'NC')
%                 NYD=size(xT);NYD(end+1:ndT-1)=1;
%                 dM=1:ndT-2;
%                 xTaux=bsxfun(@times,xT,WRes);
%                 xTaux=bsxfun(@minus,xTaux,mean(resPop(xTaux,dM,prod(NYD(dM)),1),1));
%                 num=-mean(resPop(real(bsxfun(@times,xTaux,conj(yaux))),dM,prod(NYD(dM)),1),1);
%                 xTaux=resPop(xTaux,dM,prod(NYD(dM)),1);
%                 xTauxstd=sqrt(mean(real(xTaux.*conj(xTaux)),1));
%                 den=bsxfun(@times,yauxstd,xTauxstd);
%                 E(vM,vS)=permute(num./den,permN);
                
                NYD=size(xT);NYD(end+1:5)=1;
                xTaux=reshape(xT,[prod(NYD(1:3)) NYD(4:5)]);
                xTaux=xTaux(WRes>0.5,:,:);
                xTauxstd=std(xTaux,0,1);
                xTaux=bsxfun(@minus,xTaux,mean(xTaux,1));
                num=-mean(bsxfun(@times,xTaux,yaux),1);
                den=bsxfun(@times,xTauxstd,yauxstd);                
                E(vM,vS)=permute(num./den,[3 2 1]);
            end
        end
    end
    [Emin,iMin]=min(E,[],1);
    perm=1:ndT-1;perm([2 ndT-1])=[ndT-1 2];
    iMin=permute(iMin,perm);
    T=indDim(Tcandprev,iMin,1); 
    %permute(T,[4 5 1 2 3])
end

%SHIFTED DATA
[~,kGrid,rkGrid,~,cGrid]=generateTransformGrids(NXV.*voxsiz,gpu,NXV,NXV/2+1+centShift,1);
[FT,FTH]=buildStandardDFTM(NXV,0,gpu);
for s=1:BlSz(2):NT(ndT-1);vS=s:min(s+BlSz(2)-1,NT(ndT-1));
    etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,vS,ndT-1),0,0,[],1,cGrid);
    x=dynInd(x,vS,ndT-1,real(sincRigidTransform(dynInd(x,vS,ndT-1),etInv,0,FT,FTH,0)));
end
y=mean(x,ndT-1);
