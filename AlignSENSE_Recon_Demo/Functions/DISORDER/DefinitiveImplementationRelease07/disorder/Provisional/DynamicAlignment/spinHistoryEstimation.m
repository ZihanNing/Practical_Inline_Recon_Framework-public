function x=spinHistoryEstimation(y,T,T1,TR,MB,FA,ups)

%SPINHISTORYESTIMATION performs intensity corrections due to spin history 
%in fMRI data
%   X=SPINHISTORYESTIMATION(Y,T,T1,TR,{MB},{FA})
%   * Y is an input sequence on which to perform the corrections (spatial
%   coordinates)
%   * T are the motion parameters
%   * T1 is an assumed longitudinal relaxation time for the data
%   * TR is the repeat time of the sequence
%   * {MB} is the multiband factor
%   * {FA} is the flip angle
%   * {UPS} is an upsampling factor along the third dimension
%   ** X is an output sequence with potential corrections
%   

if nargin<5 || isempty(MB);MB=1;end
if nargin<6 || isempty(FA);FA=pi/2;end
if nargin<7 || isempty(ups);ups=1;end

gpu=isa(y,'gpuArray');
NY=size(y);NY(end+1:8)=1;ND=min(numDims(y),3);
NYor=NY;
NSl=NY(3);
if ups>1
    NY(3)=NY(3)*ups;
    y=resampling(y,NY);
end
NYV=NY(1:3);
NT=size(T);ndT=ndims(T);

Mk=eye(NSl/MB,'like',real(y));
perm=1:ndT;perm([1:3 ndT-1])=[3 ndT-1 1:2];
Mk=repmat(Mk,[MB 1]);
Mk=permute(Mk,perm);
x=zeros(NYor,'like',y);

[~,kGrid,rkGrid]=generateTransformGrids(NYV,gpu,[],[],1);
[FT,FTH]=buildStandardDFTM(NYV,0,gpu);

%PULSE SHAPE GENERATION
MkUps=RFGeneration(NSl,MB,[],[],ups,gpu);
MkUps=permute(MkUps,perm);

td=TR/NT(ndT-1);%Time in between excitations
Ml=ones(NY(1:3),'like',real(y));%Longitudinal magnetization

%SPIN HISTORY SIMULATIONS
dummies=round(10*T1/TR);
cn=Mk;cn=repmat(cn,NY(1:2));
for s=1-dummies:NT(ndT-2)
    for t=1:NT(ndT-1)
        if s==1 && t==1
            cn=sum(cn,5);
            cn=mean(cn(:));
        end
        Ml=1-(1-Ml)*exp(-td/T1);
        if s>=1
            etDir=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,[s t],ndT-2:ndT-1),1,0,gpu,1);
            etInv=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,[s t],ndT-2:ndT-1),0,0,gpu,1);
            Ml=sincRigidTransform(Ml,etDir,1,FT,FTH);%Transform from material to spatial coordinates    
        end
        SlPr=dynInd(MkUps,t,5);
        Mt=bsxfun(@times,Ml,sin(FA*SlPr));%Transverse magnetization
        Ml=bsxfun(@times,Ml,cos(FA*SlPr));
        
        Mtres=resampling(Mt,NYor(1:3));
        Mtres=abs(Mtres);
        
        if s==0;cn=dynInd(cn,t,5,dynInd(cn,t,5).*(Mtres+1e-6));end
        if s>=1
            Mtres=(Mtres+1e-6)/cn;
            x=dynInd(x,s,ndT-2,dynInd(x,s,ndT-2)+sum(bsxfun(@times,Mk,Mtres),ndT-1));
            Ml=sincRigidTransform(Ml,etInv,0,FT,FTH);
            Ml=abs(Ml);
        end
    end
end
