
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

%SET DEFAULT PARAMETERS
%YB TO DO rec=espiritAlgorithm(rec);%YB TO DO
%YB TO DO Print all the parameters used here --> will be logged and can be
%accessed later

if ~isfield(rec.Alg,'writeSnapshots'); rec.Alg.writeSnapshots=0;end
if ~isfield(rec.Dyn,'Debug'); rec.Dyn.Debug=0;end
if ~isfield(rec.Alg.parE,'BlSz');rec.Alg.parE.lSz=1;end%Used in calibration matrix construction
if ~isfield(rec.Dyn,'Log'); rec.Dyn.Log=0;end%Whether to log the information
if ~isfield(rec.Alg.parE,'nIt'); rec.Alg.parE.nIt=1;end%Number of iteration in recursive ESPIRiT
if ~isfield(rec.Alg.parE,'showresult'); rec.Alg.parE.showresult = 1; end % ZN
if rec.Alg.parE.showresult == 0;rec.Alg.writeSnapshots=0;end % ZN
parE=rec.Alg.parE;
gpu=isa(rec.y,'gpuArray');

%FILENAME
folderSnapshots=strcat(rec.Names.pathOu, filesep,'Re-Se_Sn');
[~,fileName]=fileparts(generateNIIFileName(rec));

%LOGGING INFO
if rec.Dyn.Log
    logDir = strcat(rec.Names.pathOu,filesep,'Re-Se_Log', filesep);
    logName =  strcat( logDir, fileName,'_ESPIRIT_n=',num2str(rec.Alg.parE.nIt),rec.Plan.Suff,rec.Plan.SuffOu, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~exist(logDir,'dir'); mkdir(logDir);end
    diary(logName) %eval( sprintf('diary %s ', logName))
end
fprintf('\nRunning ESPIRIT on file: %s \n', strcat( fileName, rec.Plan.Suff, rec.Plan.SuffOu) )
% c = clock; fprintf('Date of reconstruction: %d/%d/%d \n', c(1:3));
tsta=tic;

%WE USE THE MEDIAN OVER THE REPEATS TO ESTIMATE
ND=length(parE.NC);%Number of dimensions
N=size(rec.y);N(end+1:12)=1;%YB: rec.y are the raw coil images
NPE=min(length(rec.Par.Mine.pedsUn),prod(N(5:12)));
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
    xf=resampling(xf,N(1:ND),2);
end
if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1%YB: (To be checked) this takes variables out of the GPU memory that are not used for a while. When they have to be used again, you simply load them again with gpuArray()
    rec.y=gather(rec.y);
    if isfield(rec,'x');rec.x=gather(rec.x);end
end
%%% Mirror image if wanted
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

%%% Set parameters for calibration size
voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
DeltaK=(1./(N(1:ND).*voxsiz(1:ND)));%1/FOV - YB: in units of 1/mm
DeltaK=DeltaK(1:ND);

KMaxCalib = 1./(2*parE.NC); %YB: in units of 1/mm
parE.NC=2*KMaxCalib./DeltaK;%Resolution (array size) of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution %Now parE.NC is in units of voxels
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
%parE.NC=parE.NC.*(1+parE.mirr(1:ND));%Mirroring
%YB: Now NC is in units of voxels to include in AC region and not resolution anymore

%%%Set parameters for kernel size
if parE.Ksph>0
    parE.K=ones(1,ND)*parE.Ksph^(1/ND);
else
    parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles 
    parE.K = round(parE.K);%YB: added
    %parE.K=parE.K-mod(parE.K,2)+1;
    parE.K=max(parE.K,parE.Kmin);
    parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution
    %parE.K=parE.K-mod(parE.K,2)+1;
    %parE.K=parE.K.*(1+parE.mirr(1:ND))-parE.mirr(1:ND);
end
NK=round(prod(parE.K));%Number of samples in kernel

if rec.Dyn.Debug>=1
    fprintf('AC region and kernel information:\n');
    fprintf('Size of calibration area:%s\n',sprintf(' %d',parE.NC));
    fprintf('Approximate size of kernels:%s\n',sprintf(' %d',round(parE.K)));
end

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('\nTime preprocessing: %.3f s\n',tend);end

rec.S=[];rec.W=[];%YB: sensitivities and eigenmaps
for pe=1:NPE
    y=dynInd(yf,pe,5);
    yst=y;%YB: used later for PI uncertainty
    if gpu;y=gpuArray(y);end
    if isfield(rec,'x')%YB: body coil
        x=dynInd(xf,pe,5);
        if gpu;x=gpuArray(x);end
    end
    N=size(y);N(end+1:4)=1;
    
    %VIRTUAL BODY COIL
    if parE.virCo>0 && isfield(rec,'x')%There must also already be a body coil in rec.x
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
                phDif=sign(multDimSum(dynInd(pa,1:N(parE.dimLoc(n)),parE.dimLoc(n)).*conj(dynInd(pa,[1 1:N(parE.dimLoc(n))-1],parE.dimLoc(n))),setdiff(1:3,parE.dimLoc(n))));                       
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
            y=cat(4,x,y);%YB: We add the virual body coil to the channel and ESPIRIT will make all channels relative to the first one (only phase I THINK)
        end
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time compressing channels: %.3f s\n',tend);end
    end
    NPm=N;
    
    %N(1:ND)=N(1:ND).*(1+parE.mirr);
    
    tsta=tic;
    for n=1:ND%YB: move to Fourier domain (kspace)
        %mirr=zeros(1,ND);mirr(n)=parE.mirr(n);
        %y=mirroring(y,mirr,1);
        y=fftGPU(y,n)/sqrt(N(n));%YB: normalised fft
        y=fftshift(y,n);   
        %if parE.virCo==-1
        %    x=fftGPU(x,n)/sqrt(N(n));
        %    x=fftshift(x,n);
        %end
    end
    if parE.absPh%YB: paper on virtual conjugate coils: add complex conjugate channels
        y=cat(4,y,fftflip(y,1:ND));
    end
    y=resampling(y,parE.NC,2);%YB:Extract AC region
    if parE.virCo==-1;xu=resampling(x,parE.NC);end
    if parE.absPh;yor=y(:,:,:,1:end/2);end%YB: Original k-space

    nc=size(y,4);
    if ND==2;NZ=size(y,3);else NZ=1;end%YB: For volumetric scans, do ESPIRiT with 3D kernels
    if isempty(parE.NCV);parE.NCV=nc;end
    assert(parE.NCV<=nc,'Not prepared to deal with number of eigenmaps bigger than the number of channels');
    NMaps=max(round(N(1:ND)./parE.subSp),parE.NC);    
    NMapsPm=NMaps;
    S=zeros([nc parE.NCV prod(NMapsPm) NZ],'like',y);%YB: first dimension are the eigenvectors
    W=zeros([1 parE.NCV prod(NMapsPm) NZ],'like',y); %YB: first dimension are the eigenvectors  
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time extracting calibration and booking memory: %.3f s\n',tend);end
        
    %HANKEL INDEXES (used for calibration matrix contruction) I think put outside z loop since for different slices, the same Hankel indices can be used
    if parE.Ksph>0
        %DeltaK=ones(1,ND);%This would make them spherical, otherwise they consider 1/FOV 
        NGR=ceil(NK.^(1/ND)*ones(1,ND).*(mean(DeltaK)./DeltaK)); %YB: radius?    
        NGD=2*NGR+1;%YB: diameter?
        
        %WEIGHTS OF THE PATCH STRUCTURE
        xx=generateGrid(NGD,gpu,NGD,ceil((NGD+1)/2));
        r=xx{1}(1);r(1)=double(0);
        for n=1:ND;r=bsxfun(@plus,r,abs(xx{n}*DeltaK(n)).^2);end%YB: compute weigth (squared L2) around center - enfores spherical area - By including deltaK, you take into account FOV
        [~,ir]=sort(r(:));%Only select samples around center - diameter too big but have to make sure yo have enough samples
        irM=ir(1:NK);
        irs=ind2subV(NGD,irM);
        
        xLim=zeros(ND,2,'like',xx{1});
        for n=1:ND;xx{n}=xx{n}(irs(:,n));xx{n}=xx{n}(:);xLim(n,:)=[min(xx{n}) max(xx{n})];end 
        sw=cat(2,xx{:});
        iX=cell(1,ND);   
        for s=1:ND
            iX{s}=1-xLim(s,1):parE.NC(s)-xLim(s,2);%YB: xLim(1) is negative, so you don't get negative index!
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    else
        %YB: to understand this part better,look at hankel.m and test with a toy example
        iw=1:NK;
        sw=ind2subV(parE.K,iw);
        iX=cell(1,ND);
        for s=1:ND
            iX{s}=0:floor(parE.NC(s)-parE.K(s));%YB: added floor
            iX{s}=bsxfun(@plus,sw(:,s),iX{s})';
        end
    end
    %YB: at this point the hankel indices are not defined for the coil
    %dimensions, so this should still repeated
    iXi=cell(1,ND);

    NCh=size(y,4);%Number of coils
    fprintf('Activate parpool in advance...\n'); c = parcluster('local'); c.NumWorkers = 10; parpool(c);% ZN: parallel pool activate in advance
    for z=1:NZ%Slices
        
        %CREATE CALIBRATION MATRIX
        tsta=tic;      
        if ND==2;yz=dynInd(y,z,3);else yz=y;end
        if parE.virCo==-1
            if ND==2;xz=dynInd(xu,z,3);else xz=xu;end
        end
        if z==NZ;y=[];end%YB: not needed anymore-remove from memory
        A=zeros([NK*NCh NK*NCh],'like',yz);%Size of symmetric matrix you are gonna construct
        NZZ=1;
        if ND==3;NZZ=size(iX{3},1);end
        NSa=size(iX{1},1)*size(iX{2},1);%First 2 dimensions - if 3D ESPIRiT, you will use BlSz to loop over it
        parE.BlSz=1;if rec.Dyn.Debug>=1;fprintf('BlSz in 3rd dimension for calibration matrix construction: %s\n',num2str(parE.BlSz,'%.0f'));end        
        for zz=1:parE.BlSz:NZZ;vZ=zz:min(zz+parE.BlSz-1,NZZ);%YB: doing this is ok since  A'A with A = [A1 ; A2]  is the same as A1'A1 + A2'A2     
            %YB: Not that choosing a different BlSz affect the stacking order of the calibration matrix Ab --> eigevectors will change
            Ab=zeros([NSa*length(vZ) NK NCh],'like', A);             
            for w=1:NK %Loop over samples in kernel and shift  (see Miki Lustig's code as well)
                for s=1:2;iXi{s}=iX{s}(:,w);end      
                if ND==3;iXi{3}=iX{3}(vZ,w);end                
                Ab(:,w,:)=reshape(dynInd(yz,iXi,1:ND),[NSa*length(vZ) 1 NCh]);%YB: because dynInd is used, hankel indices also valid across coils
            end
            Ab=Ab(:,:);%Move coil dimension to the columns  (see Miki Lustig's code as well)
            A=A+Ab'*Ab;%YB: construct A'A to then do eigenvalue decomposition on it instead of svd on A (actually using the shur). Results in same information
        end;Ab=[];
        %YB: Now A shouldbe symmetric (https://math.stackexchange.com/questions/465799/ata-is-always-a-symmetric-matrix)
        A=(A+A')/2;%YB: This is a useful technical when you have values that are intended to be symmetric but which turn out not to be due to round-off error.  
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time building Hankel (size %d) slice %d/%d: %.3f s\n',size(A,1),z,NZ,tend);end   
        if rec.Dyn.Debug>=1;fprintf('Size of theoretical Hankel (calibration) matrix: [%s %s] \n',num2str(NSa*NZZ) , num2str(NK*NCh));end   
        if rec.Dyn.Debug>=1;fprintf('Size of practical A^H*A calibration matrix to compute svd via Shur decomposition: [%s %s] \n',num2str(size(A,1)),num2str(size(A,2)) );end   
        %YB: note: by making A in the form Ab'*Ab, A is a normal matrix --> Shur decomposition will give eigenvalues and eigenvectors
        
        %%% COMPUTE SVD OF CALIBRATION MATRIX (via Shur decomposition - see comments)
        tsta=tic;
        A=gather(A);%shur.m does not support gpu functionality --> ues eig.m?
        [V,D]=schur(A);
        if gpu;[V,D]=parUnaFun({V,D},@gpuArray);end%YB: have to load into GPU again
        D=diagm(sqrt(abs(D)));%singular values are sqrt of eigenvalues of X*X and T=shur(X*X) also has same eigenvalues as X*X and they are on the diagonal (https://www.youtube.com/watch?v=cTCLCKaFzqw)
        D=flip(D,2);%We use decreasingly sorted singular values
        V=flip(V,2);
        [D,iS]=sort(D,2,'descend');  
        V=indDim(V,iS,2);
        D(D<1e-9)=1;%YB: why setting to ones? - should be not detect in screePoint with D instead of D.^2?
        NV=size(V);          
        V=reshape(V,[NK NCh NV(2)]);%Kernel  
        if parE.eigTh<=0;[~,nv]=screePoint(D.^2);else nv=find(D>=D(1)*parE.eigTh,1,'last');end%Automatic/tuned detection
        
        if rec.Dyn.Debug>=1;fprintf('Threshold for SV (relative to the largest SV): %s\n',num2str(D(nv)/D(1),'%.3f'));end
        if rec.Alg.writeSnapshots; visESPIRIT_ACSVD(V,D,nv, [0 0.03], [], strcat(folderSnapshots,filesep,'ESPIRIT_AutoCalib_SVD'),strcat(rec.Names.Name,rec.Plan.SuffOu,'_n=',num2str(rec.Alg.parE.nIt)));end
        
        %%% CROP KERNELS AND COMPUTE EIGEN-VALUE DECOMPOSITION IN IMAGE SPACE TO GET MAPS
        NDV=numDims(V);
        V=dynInd(V,1:nv,NDV);%Select number of kernels to used based on singular value thresholding  
        
        %Reshape V in kernel size and channels
        perm=[1 3 2];
        V=permute(V,perm);%Now second dimension are the number of kernels
        V=reshape(V,NK*nv,[]);%YB: why flattening this and not going to original k-space kernel size (eg. [5 5 5]) and then use the ifft ==> flatting needed for PCA
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing k-space kernels slice %d/%d: %.3f s\n',z,NZ,tend);end                
        
        %rotate kernel to order by maximum variance - see Miki Lustig's code - PCA
        tsta=tic;
        if size(V,1)<size(V,2);[~,~,v]=svd(V);else [~,~,v]=svd(V,'econ');end %YB: See ML lecture 5: if rows are samples and colums features, the eigenvectors are the V of the SVD
        V=V*v;%YB: See ML lecture 5: projection onto (orthogonal!!!) basis 
        
        %reshape filters
        if parE.Ksph>0
            Vaux=reshape(V,[NK nv nc]);
            K=xLim(:,2)-xLim(:,1)+1;
            V=zeros([K' nv nc],'like',Vaux);%YB:Filters with rectangular size, but will be filled with spherical calibration region
            swi=bsxfun(@minus,sw',xLim(:,1))+1;
            swa=cell(1,ND);
            for w=1:NK%YB: is the extraction here consistent with BlSz used previously?
                for s=1:ND;swa{s}=swi(s,w);end
                V=dynInd(V,swa,1:ND,reshape(Vaux(w,:,:),[ones(1,ND) nv nc]));
            end;Vaux=[];            
        else
            V=reshape(V,[parE.K nv nc]);%Now we go back to kernel size after PCA!
        end 
        if rec.Dyn.Debug>=1;fprintf('Size of k-space kernels:%s\n',sprintf(' %d',size(V)));end   
        
        %go to image space
        NDV=numDims(V);
        perm=1:NDV;perm([NDV-1 NDV])=[NDV NDV-1];%YB: 
        V=conj(permute(V,perm));%YB: conj because the fft will then give flipped result??? But conj is flip AND ONE UNIT DISPLACED
        for n=1:ND;V=ifftshift(V,n);end         
        V=V*(prod(NMaps)/sqrt(NK));     %YB: normalisation??   
        perm(1:2)=ND+(1:2);perm(3:ND+2)=1:ND;
        V=permute(V,perm);
        if parE.virCo==-1;xz=resampling(xz,NMaps(1:ND));end
        if ND==3%YB: ifft to go back to image space
            for n=ND
                NMapsOr=size(V);
                NMapsOr(n+2)=NMaps(n);%+2 because of permute
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);%Non-unitary matrix
                FH=resampling(FH,[NMapsOr(n+2) size(V,n+2)],1);
                V=aplGPU(FH,V,n+2);
            end
        end
        NSp=size(V,5);
        for a=1:NSp%YB: final eigenvalue decomp per voxel, so doing on a slice per slice basis is OK
            Va=dynInd(V,a,5);
            for n=1:2%ifft for 2D slice dimensions - for 3D ifft performed some lines above
                NMapsOr=size(Va);
                NMapsOr(n+2)=NMaps(n);
                [~,FH]=build1DFTM(NMapsOr(n+2),0,gpu);
                FH=resampling(FH,[NMapsOr(n+2) size(Va,n+2)],1);
                Va=aplGPU(FH,Va,n+2);
            end    
            if parE.virCo==-1;xza=dynInd(xz,a,3);xza=shiftdim(resSub(xza,1:ND),-2);end
            Va=resSub(Va,3:ND+2);
            NVa=size(Va,3);
            %YB: note: svdm will take svd of all the matrices (3rd dimension)
            [D,U]=svdm(Va);% Eigenvalue decomposition using the SVD of the non-conjugated matrix - see Miki Lusrig's code to show that this is the same as taking G'G and then eigenvalue decomp
            D=D(:,1:parE.NCV,:);U=U(:,1:parE.NCV,:);%YB: extract only the eigenmaps interested in
            D=real(D);%YB:Should be real already - check also in Miki's code
            if parE.virCo==1 && isfield(rec,'x');Uref=U(1,1,:);%YB: to normalise phase with first coil - otherwise no reference taken! (= different with Miki Lustig's code)
            elseif parE.virCo==-1;Uref=sign(xza).*sign(U(1,1,:));
            elseif parE.virCo~=-2;Uref=sign(U(1,1,:));
            else Uref=1;
            end            
            %if parE.virCo==1 && isfield(rec,'x');Uref=U(1,:,:);else Uref=sign(U(1,:,:));end
            U=bsxfun(@rdivide,U,Uref);%YB: what to reference coils to - phase AND magnitude!! - In Miki's code, only phase of reference was used
            U=matfun(@mtimes,v,U);%YB: v because rotating back after PCA - now applying transformation column in C contains contribution for every 'coil in PCA framwork'. Need to go back to original coils with basis matrix v
            D=reshape(D,[1 parE.NCV NVa]);U=reshape(U,[nc parE.NCV NVa]);
            %YB: vA are the indices of all the voxel of this slice since
            %3rd dimension in S and W is prod(N)
            vA=(1:NVa)+(a-1)*NVa;
            S(:,:,vA,z)=U;%YB: assign eigenvector (sensitivities)
            W(:,:,vA,z)=D;%YB: assifn eigenvalue - largest is first one since singular values were already sorted (see svdm.m)
        end;Va=[];V=[];U=[];D=[];
        tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing spatial maps slice %d/%d: %.3f s\n',z,NZ,tend);end     
    end;%yz=[];
    fprintf('Close parpool after computation.\n'); delete(gcp('nocreate'));  % ZN: delete parpool after computation
    
    tsta=tic;
    S=ipermute(S,perm);
    W=ipermute(W,perm);

    if ND==3    
        S=reshape(S,[NMapsPm nc parE.NCV NZ]);
        W=reshape(W,[NMapsPm 1 parE.NCV NZ]);    
    else
        S=reshape(S,[NMapsPm NZ nc parE.NCV]);
        W=reshape(W,[NMapsPm NZ 1 parE.NCV]);
    end
    
    %plot eigenmaps
    if rec.Alg.writeSnapshots; visESPIRIT_EigenMaps(W, [], 0,strcat(folderSnapshots,filesep,'ESPIRIT_EigenMaps'),strcat(rec.Names.Name,rec.Plan.SuffOu,'_n=',num2str(rec.Alg.parE.nIt)));end

    if parE.absPh%YB: need to remove the virtual conjugate coils in a correct way
        NS=size(S);NS(end+1:5)=1;
        S=resPop(S,4,[NS(4)/2 2]);%YB: split the coils - VCC now in fifth dimension
        NDS=numDims(S);
        %YB: I think this will throw error is parE.NCV >1 - Not true since
        %resPop moves dimensions to the end
        ph=exp(1i*angle(sum(dynInd(S,1,NDS).*dynInd(S,2,NDS),NDS-1))/2);%YB: This is the phase from the paper on VCC in Equation 5 - summing overl coils and then getting phase
        S=dynInd(S,1,NDS);%Only keep the original channels
        perm=1:NDS;perm([4 NDS-1])=[NDS-1 4]; 
        S=permute(S,perm);
        S=bsxfun(@times,S,conj(ph));%Remove arbitraty phase and be left with only image phase

        yor=resampling(yor,NMaps,2);
        for n=1:ND
            yor=ifftshift(yor,n);
            yor=ifftGPU(yor,n);
        end    
        yor=abs(angle(bsxfun(@times,S,conj(yor))))>pi/2;%YB: why pi/2? Since you take abs
        S(yor)=-S(yor);%YB:Allign phase 
        %yor=[];
    end
    
    if isfield(parE,'gibbsRinging')%YB: post-processing filtering
        if parE.gibbsRinging~=0 && parE.gibbsRinging<=1;HY=buildFilter(NMapsPm,'tukeyIso',[],gpu,parE.gibbsRinging);else;HY=buildFilter(NMapsPm,'CubicBSpline',[],gpu);end
        S=filtering(S,HY);
        W=filtering(W,HY);
        W=max(min(W,1),0);%YB: not tricky since eigenvalues bigger than 1 (coming from error), will be set to 1??
    end

    %%%RESAMPLE TO ORIGINAL ARRAY SIZE
    NS=size(S);NS(end+1:ND)=1;
    Sou=zeros([NPm(1:ND) NS(ND+1:end)],'single');
    for s=1:NS(4);Sou=dynInd(Sou,s,4,gather(resampling(dynInd(S,s,4),NPm(1:ND))));end
    S=Sou;Sou=[];
    W=abs(resampling(W,NPm(1:ND)));  

    if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1;S=gather(S);W=gather(W);gpuR=0;elseif gpu;S=gpuArray(S);end
    
    if parE.virCo>0
        if parE.virCo==3
            if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1;pb=gather(pb);P=gather(P);end
            S=bsxfun(@times,S,conj(pb));
            perm=[4 5 1 2 3];
            S=permute(S,perm);   
            S=matfun(@mtimes,P,S);
            S=ipermute(S,perm);
        elseif isfield(rec,'x')
            NS=size(S);
            S=dynInd(S,2:NS(4),4);%YB: remove the virtual body coil out of the list of estimated coils
        end
    end

    %SOLVE FOR PI UNCERTAINTY by looking at phase of body coil vs estiamted
    %image with guess on sensitivities
    if isfield(rec,'x')
        if gpu;yst=gpuArray(yst);end%YB: original y in image space
        if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1;yst=gather(yst);x=gather(x);end
        Sr=dynInd(S,1,5);
        xr=sum(bsxfun(@times,conj(Sr),yst),4)./(normm(Sr,[],4)+1e-6);yst=[];%Estimate of the image with the estimated sens
        if ND==3
            if abs(angle(mean(xr(:).*conj(x(:)))))>pi/2;S=-S;end%YB: allign phase
        else
            for s=1:size(x,3)
                xxr=dynInd(xr,s,3);xx=dynInd(x,s,3);
                if abs(angle(mean(xxr(:).*conj(xx(:)))))>pi/2;S=dynInd(S,s,3,-dynInd(S,s,3));end
            end
        end
    end
    
    if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1 && gpu;S=gpuArray(S);W=gpuArray(W);end
        
    for n=1:length(parE.mirr)
        if parE.mirr(n)>0
            mirr=parE.mirr(n);
            S=dynInd(S,mirr+1:size(S,n)-mirr,n);
            W=dynInd(W,mirr+1:size(W,n)-mirr,n);
        end
    end
    
    if isfield(rec,'x')
        N=size(rec.x);N(end+1:3)=1;N=N(1:3);
        S=resampling(S,N,2);
        W=resampling(W,N,2);
    end
    perm=[1:4 6 5];
    S=permute(S,perm);W=permute(W,perm);
    
    %Fat-shift
    %if rec.Par.Labels.WFS==0;SH{2}=-round(29.5507);else SH{2}=-round(rec.Par.Labels.WFS);end%Water-Fat shift  
    %Saux=shifting(dynInd(S,1,6),SH,[],1);
    %Saux=shifting(dynInd(S,1,6),SH);
    %for s=1:size(Saux,4);Saux(:,:,:,s)=extrapolating(Saux(:,:,:,s),abs(Saux(:,:,:,s))~=0,'PG-CG');end
    %S=cat(6,S,Saux);
    %W=dynInd(W,size(W,6)+1,6,0);
    
    if isfield(parE,'saveGPUMemory') && parE.saveGPUMemory==1;S=gather(S);W=gather(W);end    
    if isempty(rec.S);rec.S=S;rec.W=W;else rec.S=cat(5,rec.S,S);rec.W=cat(5,rec.W,W);end%YB: assign to structure
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing absolute phase and normalized maps: %.3f s\n',tend);end
end

%%% DATA PROJECTION CHECK
if rec.Alg.writeSnapshots; visESPIRIT_DataProjection(yz, S,single(parE.virCo>0), [],[],strcat(folderSnapshots,filesep,'ESPIRIT_SubspaceProj'),strcat(rec.Names.Name,rec.Plan.SuffOu,'_n=',num2str(rec.Alg.parE.nIt)));end
 
%%%  EXPORT INFORMATION
indM=[7 27];%YB: see disorderAlgorithm.m 7=sens 27=eigenmaps 8 =mask and 12=reconstruction
for m=indM
    if ~any(ismember(rec.Dyn.Typ2Rec,m));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,m);end
    rec.Dyn.Typ2Wri(m)=1;
end

if rec.Dyn.Log;diary off;end%Stop logging
