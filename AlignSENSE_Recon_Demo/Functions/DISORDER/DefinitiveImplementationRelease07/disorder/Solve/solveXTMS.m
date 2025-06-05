function rec=solveXTMS(rec)

%SOLVEXTMS   Performs a variable projection based aligned
%SENSE reconstruction following L Cordero-Grande, EJ Hughes, J Hutter, 
%AN Price, JV Hajnal, Three-Dimensional Motion Corrected Sensitivity 
%Encoding Reconstruction for Multi-Shot Multi-Slice MRI: Application to 
%Neonatal Brain Imaging. Magn Reson Med 79:1365-1376, 2018.
%   REC=SOLVEXTMS(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * REC is a reconstruction structure with reconstructed data and
%   surrogate information (rec.(rec.Plan.Types))
%

%GPU TYPES
typ2Rec=rec.Dyn.Typ2Rec;
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end
if gpu
    for n=typ2Rec'
        if n~=5;datTyp=rec.Plan.Types{n};rec.(datTyp)=gpuArray(rec.(datTyp));end
    end
end
on=cell(1,6);for n=1:6;on{n}=ones(1,n);end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:2)=[2 1];
for n=typ2Rec'
    if n~=5;datTyp=rec.Plan.Types{n};rec.(datTyp)=permute(rec.(datTyp),perm);end
end
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end
debug=2;

%USE ONLY THE FIRST AVERAGE
y=dynInd(rec.y,1,12);

%BUILD SHOTS
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
ktraj=unique(rec.Assign.z{2}(rec.Assign.z{8}==0),'stable');

NProfs=numel(ktraj);
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
ktraj=reshape(ktraj,[NEchos NShots]);
ktrajI=ktraj+1;ktrajI(ktrajI<1)=ktrajI(ktrajI<1)+numel(ktrajI);
kRange=rec.Enc.kRange{2};
A=extractShots(ktraj,kRange);%subTree before
if gpuIn;A=gpuArray(A);end

%BUILD PACKAGES
zFullTraj=rec.Assign.z{8}(1:NEchos:end);zFullTrajU=unique(zFullTraj);NZ=length(zFullTrajU);zFullTrajM=nan(size(zFullTraj));for u=1:NZ;zFullTrajM(zFullTraj==zFullTrajU(u))=u;end
kFullTraj=rec.Assign.z{2}(1:NEchos:end);kFullTrajU=unique(kFullTraj);NK=length(kFullTrajU);kFullTrajM=nan(size(kFullTraj));for u=1:NK;kFullTrajM(kFullTraj==kFullTrajU(u))=u;end
indJump=[1 find(diff(kFullTrajM)~=0)+1];%Jumps in k-space
kTrajJump=kFullTrajM(indJump);
NPacka=length(kTrajJump)/seqperiod(kTrajJump);
Packa=ones(1,NZ);
for n=2:NPacka;Packa(unique(zFullTrajM(indJump((n-1)*seqperiod(kTrajJump)+1):end)))=2;end
if rec.Dyn.Debug>=2
    fprintf('Number of packages: %d\n',NPacka);
    fprintf('Number of echoes: %d\n',NEchos);
    fprintf('Number of shots: %d\n',NShots);
end

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
parXT=rec.Alg.parXT;
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,y);
S=dynInd(S,1:eivaS(1),4);y=dynInd(y,1:eivaS(1),4);
NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    if parXT.perc(1)<1;fprintf('Number of compressed coil elements at%s%%: %d\n',sprintf(' %0.2f',parXT.perc*100),NS(4));else fprintf('Number of compressed coil elements: %d\n',NS(4));end
end

%SLICE PROFILE PARAMETERS
parXT.SlTh=rec.Par.Scan.AcqVoxelSize(3);%Slice thickness
parXT.SlOv=rec.Par.Scan.AcqVoxelSize(3)+rec.Par.Scan.SliceGap(1);%Slice overlap

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
NY=size(rec.y);
rec.x=zeros(NX,'like',y);
if ~any(rec.Dyn.Typ2Rec==12);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,12);rec.Dyn.Typ2Wri(12)=1;end
typ2Rec=rec.Dyn.Typ2Rec;

%PARAMETERS
voxSiz=[rec.Par.Scan.AcqVoxelSize(1:2) parXT.SlOv];
nExtern=99;
tolType={'Energy','RelativeResidual2Error'};%For CG--without vs with motion correction
tol=[parXT.tolerSolve 0];%For CG
nIt=[300 1];%For CG

%SLICE PROFILE
if parXT.threeD
    E.Sp=sliceProfile(NX(3),parXT.SlTh,parXT.SlOv,parXT);
    if gpuIn;E.Sp=gpuArray(E.Sp);end
end

%SOLVE WITHOUT MOTION
tsta=tic;
resAni=single(on{3});
[~,indMinRes]=min(voxSiz);        
resAni(indMinRes)=voxSiz(indMinRes);
indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
resAni=voxSiz./resAni;
BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
E.bS=round([2 4].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
E.Je=1;%To use joint encoding/decoding
E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NY(n);E.Uf{n}.NX=NX(n);end%Folding
EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NY(n);EH.Ub{n}.NX=NX(n);end%Unfolding
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)

%SOFT/HARD MASK
[CX,R,discMa]=buildMasking(rec.M,rec.Alg.UseSoftMasking,gpuIn);

%REGULARIZE IN SLICE DIRECTION
if parXT.threeD
    R.Fo.la=parXT.alpha(2);
    R.Fo.Fi=buildFilter([1 1 NX(3)],'2ndFiniteDiscrete',[1 1 1],gpu);
    R.Fo.mi=0;
end

%SENSITIVITIES
ncx=1:eivaS(1);
E.Sf=dynInd(S,ncx,4);E.dS=[1 ncx(end)];
if gpuIn;E.Sf=gpuArray(E.Sf);end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

yX=mean(dynInd(y,ncx,4),5);
if gpuIn;yX=gpuArray(yX);end
gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
if gibbsRing~=0 && gibbsRing<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpuIn,gibbsRing));end
nX=nIt(1);
if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

if nX~=0
    if gpuIn;rec.x=gpuArray(rec.x);end
    [rec.x,EnX]=CGsolver(yX,E,EH,P,CX,R,rec.x,nX,tolType{1},tol(1));
    rec.x=gather(rec.x);
end
if ~isempty(EnX)%Print and store final energy
    EnX=EnX/numel(yX);
    if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
    EnX=[];
end


if parXT.correct==0;estT=0;elseif parXT.correct==1;estT=1;else estT=[1 1];end
if any(estT~=0)
    rec.d=rec.x;
    if gpuIn;rec.d=gpuArray(rec.d);end
    if ~any(rec.Dyn.Typ2Rec==16);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,16);rec.Dyn.Typ2Wri(16)=1;end
    rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file
    
    %INITIALIZATION
    NA=size(A);NX=size(rec.M);NY=size(y);
    outlD=ones([on{2} NX(3) 1 NA(5)],'like',A);
    B=zeros([on{2} NX(3) on{2} NPacka],'like',A);
    for n=1:NPacka;B=dynInd(B,{find(Packa==n),n},[3 6],1);end
    
    L=length(estT);
    for l=1:L
        nReset=0;cont=0;        
        if L==2
            if l==1
                rec.T=single(zeros([on{5} NPacka 6]));
                %rec.T=dynInd(rec.T,[2 3],6:7,1);
            elseif l==2
                T=rec.T;rec.T=single(zeros([on{4} NA(5) NX(3) 6]));              
                for n=1:NPacka;rec.T=dynInd(rec.T,Packa==n,6,repmat(dynInd(T,n,6),[on{4} NA(5) sum(single((Packa==n)))]));end
            end
        else
            rec.T=single(zeros([on{4} NA(5) NX(3) 6]));
            %rec.T=dynInd(rec.T,3,7,1);
        end

        if debug==2;fprintf('Level: %d of %d\n',l,L);tstart=tic;end        

        %MOTION INFO        
        NT=size(rec.T);  
        NStates=prod(NT(5:6));
        cT=zeros(NT(1:6),'single');%Flag to indicate whether a given motion state can be assumed to have converged     
        w=parXT.winit*ones(NT(1:6),'single');%Weights of LM iteration            

        %DATA
        if L==2 && l==1       
            yXX=yX;
            E.Bm=B;%Slice mask
            E.dS(1)=1;E.NMs=1;
            
            %MOTION GRIDS
            [E.rG,E.kG,E.rkG]=generateTransformGrids(NX,gpuIn,NX,ceil((NX+1)/2),1);
            [E.Fof,E.Fob]=buildStandardDFTM(NX(1:3),0,gpuIn);            
        else
            [E.Fms,EH.Fms]=build1DFTM(NY(1),1,gpuIn);
            E.Fms=E.Fms(ktrajI(:),:);EH.Fms=EH.Fms(:,ktrajI(:));
            E.Fms=reshape(E.Fms,[NEchos NShots NProfs]);EH.Fms=reshape(EH.Fms,[NProfs NEchos NShots]);
            E.Fms=resPop(resPop(E.Fms,2,[],5),3,[],2);EH.Fms=resPop(EH.Fms,3,[],5);
            E.dS(1)=NShots;E.NMs=NShots;      
            %E.NRe=E.NSe-E.NMs;%Number of outern segments

            yXX=matfun(@mtimes,E.Fms,yX);

            E.ZSl=7;%For slab extraction
            if isfield(E,'Bm');E=rmfield(E,'Bm');end

            %SLICE PROFILE
            if parXT.threeD
                E.Sp=sliceProfile(E.ZSl,parXT.SlTh,parXT.SlOv,parXT);
                %E.Sp=sliceProfile(NX(3),parXT.SlTh,parXT.SlOv,parXT);
                if gpuIn;E.Sp=gpuArray(E.Sp);end
            end       
            
            %MOTION GRIDS            
            [E.rG,E.kG,E.rkG]=generateTransformGrids(NX,gpuIn,NX,ceil((NX+1)/2),1,[],E.ZSl);
            [E.Fof,E.Fob]=buildStandardDFTM([NX(1:2) E.ZSl],0,gpuIn);            
        end

        %CONTROL OF TRANSFORM CONVERGENCE / CONSTRAIN THE TRANSFORM / ROI IN
        %READOUT FOR MOTION ESTIMATION
        E.w=w;E.cT=cT;E.winit=parXT.winit;
        if parXT.meanT;CT.mV=5:6;else CT=[];end
        %E.tL=parXT.traLimXT*(prod(voxSiz)^(1/3));%Minimum update to keep on doing motion estimation
        E.tL=parXT.traLimXT*parXT.SlTh;
        E.Tr=rec.T;%Assign starting motion
        if isfield(E,'En');E=rmfield(E,'En');end

        %ITERATIVE SOLVER
        for n=1:nExtern
            if rec.Dyn.Debug>=2;fprintf('Iteration: %d\n',n);end
            %RESET CONVERGENCE OF MOTION
            if mod(cont,nReset)==0 || all(E.cT(:)) || parXT.convTransformJoint
                nReset=nReset+1;cont=0;%We increment the number of iterations till next reset
                E.cT(:)=0;%We reset
                E.cT(cT==1)=1;%We fix to 1 the low frequency states in case that for whatever reason cT was set to 1
                if rec.Dyn.Debug>=2;fprintf('Resetting motion convergence\n');end
            end
            if rec.Dyn.Debug>=2;fprintf('Explored motion states: %d of %d\n',NStates-sum(single(E.cT(:))),NStates);end
            cont=cont+1;
            nX=nIt(estT(l)+1);
            if n==1;nX=0;end
            
            if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       
      
            %profile on         
            if nX~=0;[rec.d,EnX]=CGsolver(yXX,E,EH,P,CX,R,rec.d,nX,tolType{estT(l)+1},tol(estT(l)+1));end
            %profile off
            %profsave(profile('info'),'../Profiling/ProfileResultsA')


            if ~isempty(EnX)%Print and store final energy
                EnX=EnX/numel(yXX);
                if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
                EnX=[];
            end
                        
            if estT(l)  
                %LM SOLVER
                if parXT.exploreMemory;E.cT(:)=1;end%To test flow without running the main methods
                E=LMsolver(yXX,E,rec.d,CT);
            else
                break
            end
            
            %RESIDUALS            
            ry=normm(encode(rec.d,E),yXX,[1:2 4])/numel(yXX);%Residuals
            if rec.Dyn.Debug>=2;fprintf('Residuals: %0.6g\n',sum(ry(:)));end
            if l==2
                EH.We=2*ry./(dynInd(ry,[2:NX(3) 1],3)+dynInd(ry,[NX(3) 1:NX(3)-1],3));
                %figure;imshow(squeeze(EH.We),[0 2])
                %pause
                EH.We=single(EH.We<parXT.outlP);
            end
            
            %CONVERGENCE          
            if all(E.cT(:))%Weak convergence
                %if l==L;estT(l)=0;else break;end%If last level we perform another estimation for x           
                estT(l)=0;    
            end
        end
    end
end

%PERMUTE BACK
typ2Rec=rec.Dyn.Typ2Rec;
for n=typ2Rec'
    datTyp=rec.Plan.Types{n};
    if n~=5;rec.(datTyp)=permute(rec.(datTyp),perm);end
    %POSTPROCESSING
    if ismember(n,[12 16])   
        if rec.Alg.MargosianFilter%PARTIAL FOURIER
            Enc=rec.Enc;
            rec.(datTyp)=margosianFilter(rec.(datTyp),Enc);   
        end
       rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);%REMOVE OVERDECODING  
    end
end
    
    
