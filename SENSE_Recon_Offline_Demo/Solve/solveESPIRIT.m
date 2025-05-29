function rec=solveESPIRIT(rec)

%SOLVEESPIRIT   Estimates the sensitivities using the ESPIRIT method based
%on [1] M Uecker, P Lai, MJ Murphy, P Virtue, M Elad, JM Pauly, SS 
%Vasanawala, M Lustig, "ESPIRiT—An eigenvalue approach to autocalibrating
%parallel MRI: where SENSE meets GRAPPA," Magn Reson Med, 71:990-1001, 
%2014, [2] M Uecker, M Lustig, "Estimating absolute-phase maps using 
%ESPIRiT and virtual conjugate coils," Magn Reson Med, 77:1201-1207 (2017),
%[3] M Uecker, P Virtue, SS Vasanawala, M Lustig, "ESPIRiT reconstruction 
%using soft SENSE," ISMRM (2013), [4] M. Buehrer, P. Boesiger, S. Kozerke,
%"Virtual body coil calibration for phased-array imaging," ISMRM (2009).
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and eigenmaps rec.W
%

if useGPU;rec.y=gpuArray(rec.y);if isfield(rec,'x');rec.x=gpuArray(rec.x);end;end


if ~isfield(rec,'Alg');rec.Alg=[];end
if ~isfield(rec.Alg,'parE');rec.Alg.parE=[];end
parE=rec.Alg.parE;

%PARAMETERS FOR ESPIRIT
if ~isfield(parE,'NCV');parE.NCV=1;end%Maximum number of eigenmaps
if ~isfield(parE,'NC');parE.NC=[4 4 4];end%Resolution (mm) of calibration area to compute compression
if ~isfield(parE,'K');parE.K=[100 100 100];end%Size of kernels
if ~isfield(parE,'eigTh');parE.eigTh=0;end%0.02%Threshold for picking singular vectors of the calibration matrix (relative to the largest singular value)
if ~isfield(parE,'absPh');parE.absPh=0;end%1%Flag to compute the absolute phase
if ~isfield(parE,'virCo');parE.virCo=1;end%Flag to use the virtual coil to normalize the maps, virCo1->URef=x, virCo5->URef=None
if ~isfield(parE,'eigSc');parE.eigSc=[0.85 0.3];end%Possible cut-off for the eigenmaps in soft SENSE, not used
if ~isfield(parE,'subSp');parE.subSp=[1 1 1];end%[1 1 1];%Subsampling in image space to accelerate
if ~isfield(parE,'dimLoc');parE.dimLoc=[];end%Dimensions along which to localize to compute virtual coils
if ~isfield(parE,'mirr');parE.mirr=[8 8 8];end%Whether to mirror along a given dimension
if ~isfield(parE,'lambda');parE.lambda=50;end%Regularization factor (not used)
if ~isfield(parE,'Kmin');parE.Kmin=6;end%Minimum K-value before mirroring
if ~isfield(parE,'Ksph');parE.Ksph=200;end%Number of points for spherical calibration area, unless 0, it overrides other K's
if ~isfield(parE,'saveGPUMemory');parE.saveGPUMemory=0;end%To save GPU memory
if ~isfield(parE,'factScreePoint');parE.factScreePoint=1;end%Factor over threshold computed with scree point
if ~isfield(parE,'factorBody');parE.factorBody=1;end%Factor for virtual body coil in ESPIRIT

if ~isfield(rec,'Par') || ~isfield(rec.Par,'Mine') || ~isfield(rec.Par.Mine,'pedsUn');rec.Par.Mine.pedsUn=1;end
if ~isfield(rec,'Enc') || ~isfield(rec.Enc,'AcqVoxelSize');rec.Enc.AcqVoxelSize=ones(1,3);end
if ~isfield(rec.Par,'Labels') || ~isfield(rec.Par.Labels,'SliceGaps');rec.Par.Labels.SliceGaps=0;end
if ~isfield(rec,'Dyn') || ~isfield(rec.Dyn,'Debug');rec.Dyn.Debug=2;end
if ~isfield(rec.Dyn,'Typ2Rec');rec.Dyn.Typ2Rec=[];end


gpu=isa(rec.y,'gpuArray');

tsta=tic;
%WE USE THE MEDIAN OVER THE REPEATS TO ESTIMATE
ND=length(parE.NC);%Number of dimensions
N=size(rec.y);N(end+1:12)=1;
NPE=min(length(rec.Par.Mine.pedsUn),prod(N(5:12)));
if isfield(rec,'x');rec.x=rec.x*parE.factorBody;end
if isfield(rec,'x');xf=permute(rec.x,[1 2 3 5 4]);end
if NPE==1
    %for s=1:N(4);yf=dynInd(yf,s,4,multDimMed(dynInd(rec.y,s,4),5:rec.Plan.NDims));end
    %if isfield(rec,'x');xf=multDimMed(xf,5:rec.Plan.NDims);end
    yf=rec.y(:,:,:,:,1);
    if isfield(rec,'x');xf=xf(:,:,:,:,1);end
else
    yf=rec.y(:,:,:,:,:);
end
if isfield(rec,'x')
    %for m=1:ND;xf=fold(xf,m,size(x,m),N(m));end
    xf=resampling(xf,N(1:ND),3);
end
if parE.saveGPUMemory
    rec.y=gather(rec.y);
    if isfield(rec,'x');rec.x=gather(rec.x);end
end

for n=1:length(parE.mirr)
    if parE.mirr(n)>0
        mirr=parE.mirr(n);
        yext1=dynInd(yf,1:mirr,n);
        yext2=dynInd(yf,N(n)-mirr+1:N(n),n);
        yf=cat(n,flip(yext1,n),yf,flip(yext2,n));
        if isfield(rec,'x')
            xext1=dynInd(xf,1:mirr,n);
            xext2=dynInd(xf,N(n)-mirr+1:N(n),n);
            xf=cat(n,flip(xext1,n),xf,flip(xext2,n));
        end
    end
end;yext1=[];yext2=[];xext1=[];xext2=[];
N=size(yf);N(end+1:12)=1;
yf=gather(yf);
if isfield(rec,'x');xf=gather(xf);end


voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
DeltaK=(1./(N(1:ND).*voxsiz(1:ND)));%1/FOV
DeltaK=DeltaK(1:ND);
parE.NC=(1./parE.NC)./DeltaK;%Resolution of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
%parE.NC=parE.NC.*(1+parE.mirr(1:ND));%Mirroring

if parE.Ksph>0
    parE.K=ones(1,ND)*parE.Ksph^(1/ND);
else
    parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles
    %parE.K=parE.K-mod(parE.K,2)+1;
    parE.K=max(parE.K,parE.Kmin);
    parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution
    %parE.K=parE.K-mod(parE.K,2)+1;
    %parE.K=parE.K.*(1+parE.mirr(1:ND))-parE.mirr(1:ND);
end
NK=round(prod(parE.K));

if rec.Dyn.Debug>=2
    fprintf('Size of calibration area:%s\n',sprintf(' %d',parE.NC));
    fprintf('Approximate size of kernels:%s\n',sprintf(' %d',round(parE.K)));
end

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time preprocessing: %.3f s\n',tend);end

rec.S=[];rec.W=[];rec.Wfull=[];
for pe=1:NPE
    y=dynInd(yf,pe,5);
    yst=y;
    if gpu;y=gpuArray(y);end
    if isfield(rec,'x')
        x=dynInd(xf,pe,5);
        if gpu;x=gpuArray(x);end
    end
    N=size(y);N(end+1:4)=1;
    
    %VIRTUAL BODY COIL
    if parE.virCo>0 && isfield(rec,'x')
        tsta=tic;
        if parE.virCo>=2 && parE.virCo<5
            %COMPRESSION
            dimNoLoc=setdiff(1:3,parE.dimLoc);
            perm=[1 2 3 5 4];    
            yH=permute(conj(y),perm);
            Ainv=bsxfun(@times,y,(sum(abs(y).^2,4)+1e-9).^(-1));
            perm=[4 5 1 2 3];           
            P=permute(multDimSum(bsxfun(@times,Ainv,yH),dimNoLoc),perm);
            yH=[];Ainv=[];
        
            NP=size(P);
            P=resSub(P,3:5);       
            [~,P]=svdm(P);
            P=reshape(P,NP);
        
            if parE.virCo==3%WE ESTIMATE IN ORTHOGONAL SPACE
                y=permute(y,perm);
                y=matfun(@mtimes,matfun(@ctranspose,P),y);
                y=ipermute(y,perm);
                pa=y;
            else%WE ADD A BODY COIL CHANNEL
                pa=permute(y,perm);   
                pa=matfun(@mtimes,matfun(@ctranspose,P),pa);
                pa=ipermute(pa,perm);
            end
            %PHASE CORRECTION
            pa=dynInd(pa,1,4);pb=pa;pb(:)=1;        
            for n=1:length(parE.dimLoc)           
                phDif=signz(multDimSum(dynInd(pa,1:N(parE.dimLoc(n)),parE.dimLoc(n)).*conj(dynInd(pa,[1 1:N(parE.dimLoc(n))-1],parE.dimLoc(n))),setdiff(1:3,parE.dimLoc(n))));                       
                pb=bsxfun(@times,pb,conj(cumprod(phDif,parE.dimLoc(n))));
                pa=bsxfun(@times,pa,pb);
            end
            if parE.virCo==3
                y=bsxfun(@times,pb,y);
            else 
                pa=bsxfun(@times,pb,pa);
                y=cat(4,pa,y);
            end
        else
            y=cat(4,x,y);
        end
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time compressing channels: %.3f s\n',tend);end
    end
    NPm=N;
    
    
    %N(1:ND)=N(1:ND).*(1+parE.mirr);
    
    tsta=tic;
    for n=1:ND
        %mirr=zeros(1,ND);mirr(n)=parE.mirr(n);
        %y=mirroring(y,mirr,1);
        y=fftGPU(y,n)/sqrt(N(n));
        y=fftshift(y,n);   
        %if parE.virCo==-1
        %    x=fftGPU(x,n)/sqrt(N(n));
        %    x=fftshift(x,n);
        %end
    end
    if parE.absPh
        y=cat(4,y,fftflip(y,1:ND));
    end
    y=resampling(y,parE.NC,3);
    if parE.virCo==-1;xu=resampling(x,parE.NC);end
    if parE.absPh;yor=y(:,:,:,1:end/2);end

    nc=size(y,4);
    if ND==2;NZ=size(y,3);else NZ=1;end
    if isempty(parE.NCV);parE.NCV=nc;end
    assert(parE.NCV<=nc,'Not prepared to deal with number of eigenmaps bigger than the number of channels');
    NMaps=max(round(N(1:ND)./parE.subSp(1:ND)),parE.NC);    
    NMapsPm=NMaps;
    S=zeros([nc parE.NCV prod(NMapsPm) NZ],'like',y);
    W=zeros([1 parE.NCV prod(NMapsPm) NZ],'like',y);   
    Wfull=zeros([1 nc prod(NMapsPm) NZ],'like',y);   
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time extracting calibration and booking memory: %.3f s\n',tend);end
        
    %HANKEL INDEXES
    if parE.Ksph>0
        %DeltaK=ones(1,ND);%This would make them spherical, otherwise they consider 1/FOV 
        NGR=ceil(NK.^(1/ND)*ones(1,ND).*(mean(DeltaK)./DeltaK));       
        NGD=2*NGR+1;
        %WEIGHTS OF THE PATCH STRUCTURE
        xx=generateGrid(NGD,gpu,NGD,ceil((NGD+1)/2));
        r=xx{1}(1);r(1)=double(0);
        for n=1:ND;r=bsxfun(@plus,r,abs(xx{n}*DeltaK(n)).^2);end
        [~,ir]=sort(r(:));
        irM=ir(1:NK);
        irs=ind2subV(NGD,irM);
        
        xLim=zeros(ND,2,'like',xx{1});
        for n=1:ND;xx{n}=xx{n}(irs(:,n));xx{n}=xx{n}(:);xLim(n,:)=[min(xx{n}) max(xx{n})];end 
        sw=cat(2,xx{:});
        iX=cell(1,ND);   
        for s=1:ND
            iX{s}=1-xLim(s,1):parE.NC(s)-xLim(s,2);
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    else
        iw=1:NK;
        sw=ind2subV(parE.K,iw);
        iX=cell(1,ND);
        for s=1:ND
            iX{s}=0:parE.NC(s)-parE.K(s);
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    end
    iXi=cell(1,ND);

    NCh=size(y,4);%Number of coils
    for z=1:NZ%Slices
        %COMPUTE CALIBRATION MATRIX
        tsta=tic;      
        if ND==2;yz=dynInd(y,z,3);else yz=y;end
        if parE.virCo==-1
            if ND==2;xz=dynInd(xu,z,3);else xz=xu;end;
        end
        if z==NZ;y=[];end
        A=zeros([NK*NCh NK*NCh],'like',yz);
        NZZ=1;
        if ND==3;NZZ=size(iX{3},1);end
        NSa=size(iX{1},1)*size(iX{2},1);
        BlSz=1;        
        for zz=1:BlSz:NZZ;vZ=zz:min(zz+BlSz-1,NZZ);         
            Ab=zeros([NSa*length(vZ) NK NCh],'like',A);             
            for w=1:NK
                for s=1:2;iXi{s}=iX{s}(:,w);end      
                if ND==3;iXi{3}=iX{3}(vZ,w);end                
                Ab(:,w,:)=reshape(dynInd(yz,iXi,1:ND),[NSa*length(vZ) 1 NCh]);
            end
            Ab=Ab(:,:);
            A=A+Ab'*Ab;
        end;Ab=[];
        A=(A+A')/2;      
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time building Hankel (size %d) slice %d/%d: %.3f s\n',size(A,1),z,NZ,tend);end                
        
        tsta=tic;
        A=gather(A);
        [V,D]=schur(A);
        if gpu;[V,D]=parUnaFun({V,D},@gpuArray);end
        D=diagm(sqrt(abs(D)));
        D=flip(D,2);%We use decreasingly sorted singular values
        V=flip(V,2);
        [D,iS]=sort(D,2,'descend');  
        V=indDim(V,iS,2);
        D(D<1e-9)=1;
        NV=size(V);          
        V=reshape(V,[NK NCh NV(2)]);%Kernel  
        if parE.eigTh<=0
            [~,nv]=screePoint(D.^2);
            nv=find(D>=parE.factScreePoint*D(nv),1,'last');
        else 
            nv=find(D>=D(1)*parE.eigTh,1,'last');
        end%Automatic/tuned detection
        if rec.Dyn.Debug==2;fprintf('Threshold for SV (relative to the largest SV): %.3f\n',D(nv)/D(1));end
        %CROP KERNELS AND COMPUTE EIGEN-VALUE DECOMPOSITION IN IMAGE SPACE TO GET MAPS
        NDV=numDims(V);
        V=dynInd(V,1:nv,NDV);        
        %ROTATE KERNEL TO ORDER BY MAXIMUM VARIANCE
        perm=[1 3 2];
        V=permute(V,perm);
        V=reshape(V,NK*nv,[]);
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing k-space kernels slice %d/%d: %.3f s\n',z,NZ,tend);end                
        if rec.Dyn.Debug>=2;fprintf('Size of k-space kernels:%s\n',sprintf(' %d',size(V)));end        
        
        tsta=tic;
        if size(V,1)<size(V,2);[~,~,v]=svd(V);else [~,~,v]=svd(V,'econ');end     
        V=V*v;
        if parE.Ksph>0
            Vaux=reshape(V,[NK nv nc]);
            K=xLim(:,2)-xLim(:,1)+1;
            V=zeros([K' nv nc],'like',Vaux);
            swi=bsxfun(@minus,sw',xLim(:,1))+1;
            swa=cell(1,ND);
            for w=1:NK
                for s=1:ND;swa{s}=swi(s,w);end
                V=dynInd(V,swa,1:ND,reshape(Vaux(w,:,:),[ones(1,ND) nv nc]));
            end;Vaux=[];            
        else
            V=reshape(V,[parE.K nv nc]);
        end        
        NDV=numDims(V);
        perm=1:NDV;perm([NDV-1 NDV])=[NDV NDV-1];
        V=conj(permute(V,perm));        
        for n=1:ND;V=ifftshift(V,n);end         
        V=V*(prod(NMaps)/sqrt(NK));        
        perm(1:2)=ND+(1:2);perm(3:ND+2)=1:ND;
        V=permute(V,perm);
        if parE.virCo==-1;xz=resampling(xz,NMaps(1:ND));end
        if ND==3
            for n=ND
                NMapsOr=size(V);
                NMapsOr(n+2)=NMaps(n);
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);
                %FH=resampling(FH,[NMapsOr(n+2) size(V,n+2)],1);
                FH=resampling(FH,[NMapsOr(n+2) size(V,n+2)],4);
                V=aplGPU(FH,V,n+2);
            end
        end
        NSp=size(V,5);
        for a=1:NSp
            Va=dynInd(V,a,5);
            for n=1:2
                NMapsOr=size(Va);
                NMapsOr(n+2)=NMaps(n);
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);
                %FH=resampling(FH,[NMapsOr(n+2) size(Va,n+2)],1);
                FH=resampling(FH,[NMapsOr(n+2) size(Va,n+2)],4);
                Va=aplGPU(FH,Va,n+2);
            end    
            if parE.virCo==-1;xza=dynInd(xz,a,3);xza=shiftdim(resSub(xza,1:ND),-2);end
            Va=resSub(Va,3:ND+2);
            NVa=size(Va,3);
            [D,U]=svdm(Va);
            Dfull=D;Dfull=real(Dfull);
            D=D(:,1:parE.NCV,:);U=U(:,1:parE.NCV,:);
            D=real(D);
            if parE.virCo==1 && isfield(rec,'x');Uref=U(1,1,:);
            elseif parE.virCo==-1;Uref=signz(xza).*signz(U(1,1,:));
            elseif parE.virCo~=-2;Uref=signz(U(1,1,:));
            else Uref=1;
            end      
            %min(angle(Uref))
            %max(angle(Uref))      
            %if parE.virCo==1 && isfield(rec,'x');Uref=U(1,:,:);else Uref=signz(U(1,:,:));end           
            U=matfun(@mtimes,v,bsxfun(@rdivide,U,Uref));
            D=reshape(D,[1 parE.NCV NVa]);U=reshape(U,[nc parE.NCV NVa]);Dfull=reshape(Dfull,[1 nc NVa]);
            vA=(1:NVa)+(a-1)*NVa;
            S(:,:,vA,z)=U;W(:,:,vA,z)=D;Wfull(:,:,vA,z)=Dfull;
        end;Va=[];V=[];U=[];D=[];Dfull=[];
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing spatial maps slice %d/%d: %.3f s\n',z,NZ,tend);end     
    end;yz=[];
    
    tsta=tic;
    S=ipermute(S,perm);
    W=ipermute(W,perm);
    Wfull=ipermute(Wfull,perm);

    if ND==3    
        S=reshape(S,[NMapsPm nc parE.NCV NZ]);
        W=reshape(W,[NMapsPm 1 parE.NCV NZ]);    
        Wfull=reshape(Wfull,[NMapsPm 1 nc NZ]);
    else
        S=reshape(S,[NMapsPm NZ nc parE.NCV]);
        W=reshape(W,[NMapsPm NZ 1 parE.NCV]);
        Wfull=reshape(Wfull,[NMapsPm NZ 1 nc]);
    end

    if parE.absPh
        NS=size(S);NS(end+1:5)=1;
        S=resPop(S,4,[NS(4)/2 2]);
        NDS=numDims(S);
        ph=exp(1i*angle(sum(dynInd(S,1,NDS).*dynInd(S,2,NDS),NDS-1))/2);
        S=dynInd(S,1,NDS);
        perm=1:NDS;perm([4 NDS-1])=[NDS-1 4]; 
        S=permute(S,perm);
        S=bsxfun(@times,S,conj(ph));

        yor=resampling(yor,NMaps,3);
        for n=1:ND
            yor=ifftshift(yor,n);
            yor=ifftGPU(yor,n);
        end    
        yor=abs(angle(bsxfun(@times,S,conj(yor))))>pi/2;
        S(yor)=-S(yor);yor=[];
    end
    
    if isfield(parE,'gibbsRinging')
        if parE.gibbsRinging~=0 && parE.gibbsRinging<=1;HY=buildFilter(NMapsPm,'tukeyIso',[],gpu,parE.gibbsRinging);else;HY=buildFilter(NMapsPm,'CubicBSpline',[],gpu);end
        S=filtering(S,HY);
        W=filtering(W,HY);W=max(min(W,1),0);
        Wfull=filtering(Wfull,HY);Wfull=max(min(Wfull,1),0);
    end

    
    NS=size(S);NS(end+1:ND)=1;
    Sou=zeros([NPm(1:ND) NS(ND+1:end)],'single');
    for s=1:NS(4);Sou=dynInd(Sou,s,4,gather(resampling(dynInd(S,s,4),NPm(1:ND))));end
    S=Sou;Sou=[];
    W=abs(resampling(W,NPm(1:ND)));  
    Wfull=abs(resampling(Wfull,NPm(1:ND)));  

    if parE.saveGPUMemory;S=gather(S);W=gather(W);Wfull=gather(Wfull);gpuR=0;elseif gpu;S=gpuArray(S);end
    
    if parE.virCo>0
        if parE.virCo==3
            if parE.saveGPUMemory;pb=gather(pb);P=gather(P);end
            S=bsxfun(@times,S,conj(pb));
            perm=[4 5 1 2 3];
            S=permute(S,perm);   
            S=matfun(@mtimes,P,S);
            S=ipermute(S,perm);
        elseif isfield(rec,'x')
            NS=size(S);
            S=dynInd(S,2:NS(4),4);
        end
    end

    %SOLVE FOR PI UNCERTAINTY
    if isfield(rec,'x')
        if gpu;yst=gpuArray(yst);end
        if parE.saveGPUMemory;yst=gather(yst);x=gather(x);end
        Sr=dynInd(S,1,5);
        xr=sum(bsxfun(@times,conj(Sr),yst),4)./(normm(Sr,[],4)+1e-6);yst=[];
        if ND==3
            if abs(angle(mean(xr(:).*conj(x(:)))))>pi/2;S=-S;end
        else
            for s=1:size(x,3)
                xxr=dynInd(xr,s,3);xx=dynInd(x,s,3);
                if abs(angle(mean(xxr(:).*conj(xx(:)))))>pi/2;S=dynInd(S,s,3,-dynInd(S,s,3));end
            end
        end
    end
    
    if parE.saveGPUMemory && gpu;S=gpuArray(S);W=gpuArray(W);Wfull=gpuArray(Wfull);end
        
    for n=1:length(parE.mirr)
        if parE.mirr(n)>0
            mirr=parE.mirr(n);
            S=dynInd(S,mirr+1:size(S,n)-mirr,n);
            W=dynInd(W,mirr+1:size(W,n)-mirr,n);
            Wfull=dynInd(Wfull,mirr+1:size(Wfull,n)-mirr,n);
        end
    end
    
    if isfield(rec,'x')
        N=size(rec.x);N(end+1:3)=1;N=N(1:3);
        S=resampling(S,N,3);
        W=resampling(W,N,3);
        Wfull=resampling(Wfull,N,3);
    end
    perm=[1:4 6 5];
    S=permute(S,perm);W=permute(W,perm);Wfull=permute(Wfull,perm);
    
    %Fat-shitt
    %if rec.Par.Labels.WFS==0;SH{2}=-round(29.5507);else SH{2}=-round(rec.Par.Labels.WFS);end%Water-Fat shift  
    %Saux=shifting(dynInd(S,1,6),SH,[],1);
    %Saux=shifting(dynInd(S,1,6),SH);
    %for s=1:size(Saux,4);Saux(:,:,:,s)=extrapolating(Saux(:,:,:,s),abs(Saux(:,:,:,s))~=0,'PG-CG');end
    %S=cat(6,S,Saux);
    %W=dynInd(W,size(W,6)+1,6,0);
    %Wfull=dynInd(Wfull,size(Wfull,6)+1,6,0);
    
    if parE.saveGPUMemory;S=gather(S);W=gather(W);Wfull=gather(Wfull);end    
    if isempty(rec.S);rec.S=S;rec.W=W;rec.Wfull=Wfull;else rec.S=cat(5,rec.S,S);rec.W=cat(5,rec.W,W);rec.Wfull=cat(5,rec.Wfull,Wfull);end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing absolute phase and normalized maps: %.3f s\n',tend);end
end

indM=[7 27];
for m=indM
    if ~any(ismember(rec.Dyn.Typ2Rec,m));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,m);end
    rec.Dyn.Typ2Wri(m)=1;
end
