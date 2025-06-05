function [x,par,F,xF]=brainSegmentation(x,c,rRange,K,S,r,useLike,lambdaU,diUse,kappa,extK,grTh,eroDilaFactor)

%BRAINSEGMENTATION performs a volumetric fetal brain detection based on an 
%extension of [1] L Cordero-Grande et al, "Unsupervised 4D myocardium 
%segmentation with a Markov Random Field based deformable model," MedIA, 
%15:283-301, 2011 to spherical coordinates
%   X=BRAINSEGMENTATION(X,ESTPAR,C,RRANGE,K,S,{R},{USELIKE},{LAMBDAU},{DIUSE},{KAPPA})
%   * X is the volume to segment
%   * C is an estimation of the center of the brain
%   * RRANGE is the range of valid radius in pixels
%   * K is the number of radii for discretizing the problem
%   * S is the number of points in the angular domain
%   * {R} is an estimate of the average radious for the starting conditions. 
%   It defaults to the average of rRange
%   * {USELIKE} controls the treatment of the likelihood terms. -1 implies that
%   the term is not used. 0 implies that it is used but not dinamically updated. 
%   1 implies that the parameters is dinamically updated using global estimates. 
%   2 that it is dinamically updated using local estimates. It defaults to all 0
%   * {LAMBDAU} is the regularization parameter
%   * {DIUSE} is the directionality gradient feature to be used
%   * {KAPPA} controls the conduction of anisotropic diffusion respectively for 
%   the intensity and the gradient
%   * {EXTK} is the number of radial locations used to extend the feature 
%   extraction sets for boundary features
%   * {ERODILAFACTOR} erosion/dilation factor in voxels
%   * X is a set of images with (1) the original data, (2) mask image with the 
%   MRF segmented brain, (3) MRF segmented brain after morphological operations 
%   for smoother mask boundaries, (4) an ellipse fitting the result.
%   * PAR are the 9 parameters corresponding to the fitted ellipsoid. They are
%   specified as c1,c2,c3,r1,r2,r3,th1,th2,th3 (center,radious,rotation, 
%   respectively in pixels and degrees). They correspond to the parameters used 
%   to generate the ellipsoid in the BUILDGEOMETRY method
%   * F are the spatial likelihood distribution arranged as a cell
%   * xF is the filtered data arranged as a cell
%

%DEFAULT VALUES
if nargin<6 || isempty(r);r=mean(rRange);end
if nargin<7 || isempty(useLike);useLike=[0 0 0 0];end


%PRIORS:
%1->SQUARED RADII DIFFERENCES
%2->SQUARED RADII DIFFERENCES NORMALIZED BY THE RADII
%3->SQUARED DISTANCES
%4->SQUARED DISTANCES NORMALIZED BY THE RADII
%5->ANGLE BETWEEN ADJACENT FACES
%6->CONVEXITY
if nargin<8 || isempty(lambdaU);lambdaU=[0 0 0 0 0 0];end
if nargin<9 || isempty(diUse);diUse=1;end
if nargin<10 || isempty(kappa);kappa=[5 5];end
if nargin<11 || isempty(extK);extK=[2 2];end
if nargin<12 || isempty(grTh);grTh=[0.2 0.8;%Threshold (relative to the  maximum gradient) for gradient feature extraction
                                    0.05 0.95];end%Percentile for gradient intensity feature extraction
if nargin<13 || isempty(eroDilaFactor);eroDilaFactor=[6 8];end


%CHECKS AND GENERAL PURPOSE DEFINES
assert(ndims(x)==3,'The image has to be 3D and it is %d',ndims(x));
N=size(x);
gpu=isa(x,'gpuArray');
on3=ones(1,3);

%HARDCODED PARAMETERS
useRobust=[0 0 0 0];%Indicates whether to use robust estimates of the parameters
nItAn=[32 32];%Number of iterations of anisotropic diffusion based feature extration respectively for the intensity and the gradient
nIt=20;%Maximum number of iterations to estimate the field
%typOrd='Ran';%typOrd='Seq';typOrd='Ran';typOrd='HCF';typOrd='LCF';Orders used for traversing the field (SEQ->Sequential, RAN->Random, HCF->Highest confidence first, LCF->Lowest confidence first
typOrd='Seq';%typOrd='Seq';typOrd='Ran';typOrd='HCF';typOrd='LCF';Orders used for traversing the field (SEQ->Sequential, RAN->Random, HCF->Highest confidence first, LCF->Lowest confidence first%RECENT CHANGE
erode=-eroDilaFactor(1)*on3;%Erosion in voxel units
dilat=eroDilaFactor(2)*on3;%Dilation in voxel units
inert=0;%1 to compute the inertia ellipsoid instead of doing an ellipsoid fitting

%DEFINITION OF CANDIDATE RADII
rCand=linspace(rRange(1)^1,rRange(2)^1,K);
rCand=rCand.^(1/1);
rCand=permute(rCand,[1 3 2]);%Along the third dimension

%DISTANCE BETWEEN CANDIDATE RADII
d=rCand(2)-rCand(1);

%DEFINITION OF GRID ANGLES. WE DISTRIBUTE THEM UNIFORMLY ON THE SPHERE 
%USING THE VOGEL'S METHOD
[pCand,fCand,fpCand,ppCand,nCand,cCand,ff2Cand,ff1Cand]=vogelSphereGrid(S);
pCand=pCand';

%AVERAGE DISTANCES AMONG NORMALS
nCandn=sqrt(sum(nCand.^2,2));
nCandn=bsxfun(@rdivide,nCand,nCandn);
NC=0;
mdr=0;
xdr=[];
for s=1:S
    sF=fpCand{s};%Faces the point belongs to
    vn=nCandn(sF,:);
    MM=length(sF);    
    aa=nchoosek(1:MM,2);    
    dr=dot(vn(aa(:,1),:),vn(aa(:,2),:),2);   
    dr(dr>1)=1;dr(dr<-1)=-1;
    dr=acos(dr)/pi;
    dr=dr.^2;
    xdr=cat(1,xdr,dr);
    mdr=mdr+sum(dr,1);    
    NC=NC+size(aa,1);
end
mdr=mdr/NC;

%CANDIDATE POINTS
xCand=bsxfun(@times,pCand,rCand);

%RADIAL COORDINATES OF THE IMAGE 
xGrid=generateGrid(N,gpu,N,[0 0 0],-c);
rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));

%ANGULAR COORDINATES OF THE IMAGE
aGrid=cell(1,3);
for n=1:3;aGrid{n}=bsxfun(@rdivide,xGrid{n},rGrid);end
aGrid=cat(4,aGrid{:});

%CLOSEST ANGULAR POINT FOR EACH VOXEL IN THE IMAGE
aInd=sum(bsxfun(@times,aGrid,permute(pCand,[3 4 5 1 2])),4);
aInd(aInd>1)=1;aInd(aInd<-1)=-1;
[~,aInd]=min(acos(aInd),[],5);

%RADIAL BOUNDARY POINT FOR EACH VOXEL IN THE IMAGE
rbInd=abs(bsxfun(@minus,rGrid,permute(rCand,[1 2 4 3])));
[rbVal,rbInd]=min(rbInd,[],4);
rbInd(rbVal>d)=0;%Points outside the search range

%RADIAL EXTERIOR POINTS FOR EACH VOXEL IN THE IMAGE
reInd=bsxfun(@le,rGrid,permute(rCand,[1 2 4 3]));
[veVal,reInd]=max(reInd,[],4);
reInd(veVal==0)=0;

%ROI COMPUTATION
fullMask=rbInd~=0 | reInd~=0;
ROI=computeROI(fullMask,1);
fprintf('ROI segmentation:\n%s',sprintf(' %d %d %d %d %d %d\n',ROI'));

[x,reInd,rbInd,rGrid,aInd,aGrid]=parUnaFun({x,reInd,rbInd,rGrid,aInd,aGrid},@extractROI,ROI,1,1:3);

%ANISOTROPIC DIFFUSION
xa=real(anidiffND(x,nItAn(1),kappa(1)));%Intensity features
xb=real(anidiffND(x,nItAn(2),kappa(2)));%Gradient features

xF{1}=xa;
xF{2}=xb;
%if lambdaU(6)==5
%
%writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'DiffANew','DiffBNew'},{xa,xb});
%pause
%
%end

%COMPUTATION OF THE GRADIENT
[xb,xbd]=geometricFeatures(xb,grTh(1,1));
xbd=dot(xbd,aGrid,4);
if diUse==1    
    xbd=abs(xbd)-xbd/2;
elseif diUse==2
    xbd=-xbd;
    %xbd=abs(xbd);
end

%si(si==0)=1;
xb=xb.*xbd;%(1-si);
%xbaux=sort(xb(:));xbth=xbaux(floor(end/2));
%xbth
xbth=grTh(1,:)*max(xb(:));

%VOXELS USED FOR INITIAL ESTIMATES
mask{1}=rGrid<r;%OBJECT VOXELS
mask{2}=rGrid<r;%OBJECT VOXELS
gmaux=xb;
gmaux(xb<xbth(1) | xb>xbth(2) | rbInd==0)=0;

%FEATURE FOR BOUNDARY INTENSITY EXTRACTION
mask{3}=rGrid<r;
xs=sort(x(mask{3}));
%maxs=xs(ceil(0.95*numel(xs)));

%mask{3}=x>maxs*grTh(2,1) & x<maxs*grTh(2,2);%(gmaux~=0;%BOUNDARY VOXELS
mask{3}=x>xs(max(ceil(numel(xs)*grTh(2,1)),1)) & x<xs(ceil(numel(xs)*grTh(2,2)));

mask{4}=gmaux~=0;%BOUNDARY VOXELS

%FEATURES
D{1}=xa;%INTENSITY
D{2}=xb;%GRADIENT
D{3}=x;%INTENSITY
D{4}=xb;%GRADIENT

%INITIAL ESTIMATES
mest=zeros([4 2],'like',xa);
for n=1:4;mest(n,:)=momentsEstimate(D{n}(mask{n}),useRobust(n));end

%LOGLIKELIHOOD COMPUTATIONS
F=cell(1,4);
Fnm=zeros(1,4,'like',xa);
for n=1:4
    F{n}=loggauss(D{n},mest(n,:));
    Fnm(n)=min(F{n}(:));
end

%SUBSETS
gamma=cell(4,S,K);%OBJECT/OBJECT/BOUNDARY/BOUNDARY
reInd0=reInd~=0;
rbInd0=rbInd~=0;
reIndK=cell(1,K);
rbIndK=cell(2,K);
for k=1:K
    reIndK{k}=reInd<=k & reInd0;
    for n=1:length(extK)
        rbIndK{n}{k}=abs(rbInd-k)<extK(n) & rbInd0;
    end
end
for s=1:S
    aIndsa=aInd==s;
    aIndsb=aIndsa;% | ismember(aInd,ppCand{s});
    for k=1:K
        gamma{1}{s}{k}=find(aIndsa & reIndK{k});
        gamma{4}{s}{k}=find(aIndsb & rbIndK{1}{k});
        gamma{3}{s}{k}=find(aIndsb & rbIndK{2}{k});
    end
end
reIndK=[];rbIndK=[];
gamma{2}=gamma{1};

%POTENTIALS
vY=cell(1,4);
for n=1:4
    if useLike(n)~=-1
        vY{n}=Fnm(n)*ones([S K],'like',xa);
        for s=1:S
            for k=1:K
                FH=F{n}(gamma{n}{s}{k});
                if ~isempty(FH);vY{n}(s,k)=mean(FH);end
            end
        end
    end
end

%INITIALIZATION TO ESTIMATE THE FIELD
[~,iF]=min(abs(bsxfun(@minus,rCand,r)),[],3);
xR=iF*ones([1 S],'like',xa);%xR=randi(K,[1 S],'like',xa);
noconf=xR;noconf(:)=0;%Inverse of confidence of results
xS=x;
lambda=lambdaU;
NP=length(lambda);%Number of priors
vC=cell(1,NP);for n=1:NP;vC{n}=zeros(1,K,'like',xa);end%Containers for the priors
vH=zeros(1,K,'like',xa);
Beta=inf;%Inverse of the temperature
kk=1:K;%Candidate radii

for n=1:nIt   
    %mest([1 4],:)
%     figure    
%     xCurr=indDim(xCand,xR,3)';
%     trimesh(fCand,xCurr(:,1),xCurr(:,2),xCurr(:,3))
%     pause
%     %xR
    
    %SOLVING
    if n==1 && any(useLike~=-1);lambda(:)=0;else lambda=lambdaU;end
    xRPrev=xR;
    if strcmp(typOrd,'Seq');sRand=1:S;
    elseif strcmp(typOrd,'Ran');sRand=randperm(S);
    elseif strcmp(typOrd,'HCF');[~,sRand]=sort(noconf);
    else [~,sRand]=sort(noconf,2,'descend');
    end
    kRand=rand(1,S);
    for t=1:S
        s=sRand(t);    
        sP=ppCand{s};%Neighboring points
        sF=fpCand{s};%Faces the point belongs to
        
        if any(lambda(1:4)~=0)
            rCurr=permute(rCand(xR),[1 3 2]);%Set of radii - S second dimension, K third dimension
            rCurrkk=rCand(kk);
            rCurr=rCurr(:,sP,:);            
        end
        if any(lambda(3:4)~=0)            
            xCurr=indDim(xCand,xR,3);%Set of points - S second dimension, K third dimension
            xCurrkk=xCand(:,s,kk);
            xCurr=xCurr(:,sP,:);            
        end
        if any(lambda(5:6)~=0)
            xxCurr=indDim(xCand,xR,3);
            xxCurrkk=xCand(:,s,kk);
            xxCurr=repmat(xxCurr,[1 1 K]);
            xxCurr(:,s,:)=xxCurrkk;
            psF=fCand(sF,:)';%Points corresponding to relevant faces             
            psF=bsxfun(@plus,psF,permute((kk-1)*S,[1 3 2]));
            NN=size(xxCurr);
            xxCurr=reshape(xxCurr,[NN(1) prod(NN(2:3))]);           
            MM=size(psF);
            psF=reshape(psF,[MM(1) prod(MM(2:3))]);
            [v,c]=triangleNormals(xxCurr',psF');v=v';c=c';            
            v=reshape(v,MM);c=reshape(c,MM);
            vn=sqrt(sum(v.^2,1));
            vn=bsxfun(@rdivide,v,vn);
            aa=nchoosek(1:MM(2),2);
        end            
        vH(:)=0;
        for p=1:NP%Priors
            if lambda(p)~=0
                if ismember(p,1:2)
                    dr=bsxfun(@minus,rCurrkk,rCurr);
                end
                if ismember(p,3:4)
                    dr=bsxfun(@minus,xCurrkk,xCurr);
                    dr=sqrt(sum(dr.^2,1));
                end
                if p==2 || p==4;dr=bsxfun(@rdivide,dr,rCurr);end 
                
                if ismember(p,5)
                    dr=dot(vn(:,aa(:,1),:),vn(:,aa(:,2),:),1);
                    dr(dr>1)=1;dr(dr<-1)=-1;
                    dr=acos(dr)/pi;                    
                end
                if ismember(p,6)
                    aac=cat(1,aa(:,1),aa(:,2));
                    aacf=fftshift(aac,1);
                    dr=c(:,aac,:)-c(:,aacf,:);
                    dr=bsxfun(@rdivide,dr,sqrt(sum(dr.^2,1)));
                    dr=dot(vn(:,aacf,:),dr);
                    dr(dr>1)=1;dr(dr<-1)=-1;
                    dr(dr<0)=0;
                end                    
                dr=dr.^2;        
                if p~=5
                    vC{p}=lambda(p)*permute(mean(dr,2),[1 3 2]);
                else
                    dr=mean(dr,2);   
                    dr=max(dr-mdr,0);
                    vC{p}=lambda(p)*permute(dr,[1 3 2]);
                end                    
                vH=vH+vC{p};
            end
        end 
        for v=1:4
            if useLike(v)~=-1;vH=vH-vY{v}(s,:);end
        end
        if ~isinf(Beta)%Gibbs sampling
            vH=vH.^Beta;
            vH=cumsum(vH);
            vH=bsxfun(@rdivide,vH,vH(end));
            xR(s)=find(kRand(s)<=vH,1);
        else 
            vH=bsxfun(@rdivide,vH,sum(vH));
            [noconf(s),xR(s)]=min(vH,[],2);            
        end
    end
    if all(xR==xRPrev);break;end    

    %PARAMETER ESTIMATION
    for v=1:4
        if useLike(v)==1
            xS(:)=0;
            for s=1:S;xS(gamma{v}{s}{xR(s)})=1;end
            mest(v,:)=momentsEstimate(D{v}(xS==1),useRobust(n));
            F{v}=loggauss(D{v},mest(v,:));
            vY{v}(:)=min(F{v}(:));
            for s=1:S
                for k=1:K
                    FH=F{v}(gamma{v}{s}{k});
                    if ~isempty(FH);vY{v}(s,k)=mean(FH);end
                end
            end 
        elseif useLike(v)==2
            for s=1:S
                xS(:)=0;
                xS(gamma{v}{s}{xR(s)})=1;
                mest(v,:)=momentsEstimate(D{v}(xS==1),useRobust(n));
                F{v}=loggauss(D{v},mest(v,:));
                vY{v}(:)=min(F{v}(:));
                for k=1:K
                    FH=F{v}(gamma{v}{s}{k});
                    if ~isempty(FH);vY{v}(s,k)=mean(FH);end            
                end
            end 
                
        end
    end 
end
fprintf('Number of iterations of segmentation: %d\n',n);

%CONSTRUCT THE SOLUTION
xS(:)=0;
for s=1:S;xS(gamma{1}{s}{xR(s)})=1;end
[x,xS]=parUnaFun({x,xS},@extractROI,ROI,0,1:3);

%MORPHOLOGICAL OPERATIONS
xa=morphFourier(xS,erode,on3);
xa=morphFourier(xa,dilat,on3);

%ELLIPSOID FIT
[par,xb]=ellipsoidFromImage(xa,[],inert);

%WE CONCATENATE THE OUTPUTS: Image / MRF segm / Morphological segm / Ellipsoid segm
x=cat(4,x,xS,xa,xb);

if nargout>=3
    %TO RETURN LIKELIHOOD ESTIMATES. NOTE THAT THEY DON'T MAKE MUCH SENSE FOR USELIKE=2
    for v=1:4
        if ~useLike(v)~=-1;F{v}=extractROI(F{v},ROI,0,1:3);end
    end
end
if nargout>=4
    for v=1:length(xF)
        xF{v}=extractROI(xF{v},ROI,0,1:3);
    end
end

function x=loggauss(x,mest)
    x=x-mest(1);
    x=x.^2;
    x=x/(2*mest(2)^2);
    x=-x-log(sqrt(2*pi)*mest(2)^2);
end

function mest=momentsEstimate(x,rob)
    if rob;mest(1)=median(x);mest(2)=1.4826*mad(x);else mest(1)=mean(x);mest(2)=std(x);end
end

end


