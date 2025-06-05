function x=solveXSliceRegistration(x,y,T,Mk,Hsm,lambda,nX,toler,NZ,usePrecond)

%SOLVEXSLICEREGISTRATION   Computes the reconstruction x for a given T
%   X=SOLVEXSLICEREGISTRATION(X,Y,T,MK,HSM,LAMBDA,KGRID,RKGRID,NX,TOLER,NZ,USEPRECOND)
%   * X is the current reconstruction
%   * Y is the measured data
%   * T are the transform parameters
%   * MK is a mask for slice sampling
%   * HSM is a second order smoothness regularizer
%   * LAMBDA is a regularization parameter
%   * NX is the number of iterations of the CG procedure
%   * TOLER is the tolerance of the CG procedure
%   * {NZ} indicates whether to use slab extraction to accelerate
%   * {USEPRECOND} indicates whether to use a preconditioner
%   ** X is the result
%

if nargin<9;NZ=[];end
if nargin<10 || isempty(usePrecond);usePrecond=0;end

%fprintf('Slice registration\n')
%fprintf('Size of reconstruction:%s\n',sprintf(' %d',size(x)));
%fprintf('Size of measurement:%s\n',sprintf(' %d',size(y)));
%fprintf('Size of motion:%s\n',sprintf(' %d',size(T)));
%fprintf('Size of mask:%s\n',sprintf(' %d',size(Mk)));
%fprintf('Size of regularizer:%s\n',sprintf(' %d',size(Hsm)));

gpu=isa(x,'gpuArray');
NT=size(T);ndT=ndims(T);
NY=size(y);NY(end+1:8)=1;ND=min(numDims(y),3);
NYV=NY(1:3);

if ~isempty(NZ)
    [~,kGrid,rkGrid]=generateTransformGrids(NYV,gpu,[],[],1,[],NZ,[],5);
    [FT,FTH]=buildStandardDFTM([NYV(1:2) NZ],0,gpu);
else
    [~,kGrid,rkGrid]=generateTransformGrids(NYV,gpu,[],[],1);
    [FT,FTH]=buildStandardDFTM(NYV,0,gpu);
end

%BlSz=NT(5:6);
etDir=cell(1);etInv=cell(1);
BlSz=[1 NT(ndT-1)];
for s=1:BlSz(1):NT(ndT-2);vS=s:min(s+BlSz(1)-1,NT(ndT-2));
    for t=1:BlSz(2):NT(ndT-1);vT=t:min(t+BlSz(2)-1,NT(ndT-1));
        etDir{s}{t}=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,{vS,vT},ndT-2:ndT-1),1,0,gpu,1);
        etInv{s}{t}=precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(T,{vS,vT},ndT-2:ndT-1),0,0,gpu,1);
    end
end
       
yEnd=single(zeros(NY(1:ND)));
if gpu;yEnd=gpuArray(yEnd);end

Precondy=yEnd;Precondy(:)=1;
if usePrecond;Precond=yEnd;end

for s=1:BlSz(1):NT(ndT-2);vS=s:min(s+BlSz(1)-1,NT(ndT-2));
    for t=1:BlSz(2):NT(ndT-1);vT=t:min(t+BlSz(2)-1,NT(ndT-1));   
        if isempty(NZ);yAux=bsxfun(@times,dynInd(y,vS,ndT-2),dynInd(Mk,vT,ndT-1));else;yAux=extractSlabs(dynInd(y,vS,ndT-2),NZ,1,0,5);end      
        yAux=sincRigidTransform(yAux,etInv{s}{t},0,FT,FTH,0);
        if ~isempty(NZ);yAux=extractSlabs(yAux,NZ,0,0,5);end
        yEnd=yEnd+multDimSum(yAux,ndT-2:ndT-1);%else yEnd=yEnd+groupedSincRigidTransform(yAux,etInv{s}{t},0,FT,FTH);end
        if usePrecond
            yAux=bsxfun(@times,Precondy,dynInd(Mk,vT,ndT-1));
            Precond=Precond+multDimSum(sincRigidTransform(yAux,etInv{s}{t},0,FT,FTH,0),ndT-2:ndT-1);
        end            
    end
end
if usePrecond;Precond=1./Precond;else Precond=Precondy;end

y=yEnd;

%computeEnergy(x)
Ap=applyCG(x);
r=y-Ap; 

z=Precond.*r;
p=z; 
rsold=sum(conj(z(:)).*r(:));

%Iterations
n=1;
while 1
    Ap=applyCG(p);
    al=conj(rsold)/sum(conj(p(:)).*Ap(:));
    xup=al*p;
    x=x+xup;
    xup=real(xup.*conj(xup));
    xup=max(xup(:));
    %fprintf('Iteration CG %04d - Error %0.2g \n',n,xup);
    if xup<toler || n>=nX
        %if toler~=0;fprintf('Iteration CG %04d - Error %0.2g \n',n,xup);end
        break
    end
    %computeEnergy(x)
    r=r-al*Ap;
    z=Precond.*r;
    rsnew=sum(conj(z(:)).*r(:));
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew;
    if sqrt(abs(rsnew))<1e-10;break;end
    n=n+1;   
end

%%
function x=applyCG(x) 
    xB=x;
    x=lambda*filtering(x,Hsm);    
    for s=1:BlSz(1):NT(ndT-2);vS=s:min(s+BlSz(1)-1,NT(ndT-2));
        for t=1:BlSz(2):NT(ndT-1);vT=t:min(t+BlSz(2)-1,NT(ndT-1));%NOT SURE ABOUT CONNECTION OF SLABS NON-FULL BLOCK SIZES, NOT IMPORTANT WITH CHOSEN BLSZ(2)
            if ~isempty(NZ);xAux=extractSlabs(xB,NZ,1,1,5);else xAux=xB;end        
            xAux=sincRigidTransform(xAux,etDir{s}{t},1,FT,FTH);
            if isempty(NZ);xAux=bsxfun(@times,dynInd(Mk,vT,ndT-1),xAux);else xAux=extractSlabs(xAux,NZ,0,1,5);xAux=extractSlabs(xAux,NZ,1,0,5);end
            xAux=sincRigidTransform(xAux,etInv{s}{t},0,FT,FTH,0);
            if ~isempty(NZ);xAux=extractSlabs(xAux,NZ,0,0,5);end            
            x=x+multDimSum(xAux,ndT-2:ndT-1);            
        end
    end
end

end
