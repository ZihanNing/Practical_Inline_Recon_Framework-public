function rec=brainTracking(rec,appl)

% BRAINTRACKING performs temporal tracking of volumes
%   REC=BRAINTRACKING(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * REC is a reconstruction structure with information in tracked coordinates
%

if nargin<2 || isempty(appl);appl=0;end

if rec.Fail;return;end

gpu=rec.Dyn.GPU;
if ~appl
    NDims=numDims(rec.w);NDims=min(NDims,3);
    voxsiz=rec.Par.Scan.AcqVoxelSize(1,1:NDims);
    voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
    distan=5*ones(1,3);%For soft masking
    N=size(rec.w);N(end+1:4)=1;

    if ~isempty(regexp(rec.Names.Name,'dhcp','ONCE')) || ~isempty(regexp(rec.Names.Name,'wip','ONCE')) || ~isempty(regexp(rec.Names.Name,'dev','ONCE'))%NOT A FETAL CASE-SAFER TO GO THIS ROUTE
        mirr=zeros(1,3);
        rec.b=rec.w;
        if gpu;rec.b=gpuArray(rec.b);end
        rec.M=refineMask(multDimSum(abs(rec.b),4:5),rec.Alg.parS,voxsiz);
        %SOFT MASK, SIXTH COMPONENT
        Msoft=morphFourier(rec.M,distan,voxsiz,mirr,1);%SOFT MASK USED TO COMPUTE THE ROI
        Msoft(Msoft>1-1e-6)=1;Msoft(Msoft<1e-6)=0;    
        rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,[8;26]);rec.Dyn.Typ2Wri([8 26])=1;

        %ROI COMPUTATION
        rec.Enc.ROI=computeROI(Msoft);
        rec.Enc.ROI(3,:)=[1 N(3) N(3) N(3) 0 0];%To avoid problems later on with MB replication
        fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));

        %ROI EXTRACTION
        typ2Rec=rec.Dyn.Typ2Rec;
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            if ~ismember(n,24);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1:3);end
        end

        %TRANSFORM
        ND=numDims(rec.b);ND=max(ND,4);
        rec.T=single(zeros([ones(1,ND-1) N(ND) 6]));    
    else    
        AdHocArray=rec.Par.Mine.AdHocArray;
        assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
        if AdHocArray(1)==101;MB=AdHocArray(4);else MB=1;end

        NE=length(rec.Par.Labels.TE);
        if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
        if  gpu;rec.w=gpuArray(rec.w);end
        x=abs(rec.w);
        rec.w=gather(rec.w);
        N=size(x);N(end+1:4)=1;
        if NE>1 && numDims(rec.w)<=4;x=reshape(x,[N(1:3) NE N(4)/NE]);end
        if NE==1 && numDims(rec.w)>=5;x=reshape(x,[N(1:3) prod(N(4:end))]);end

        %if numDims(rec.u)>4;fprintf('Motion correction not defined for %d dimensions',numDims(rec.u));return;end
        ND=numDims(x);ND=max(ND,4);
        if ND==5;x=permute(x,[1 2 3 5 4]);end
        N=size(x);N(end+1:4)=1;
        mirr=[0 0 0 1];mirr(end+1:ND)=0;

        %APPROXIMATION TO THE CENTER
        x=gather(x);
        feat=multDimMed(x,4:ND);
        if gpu;feat=gpuArray(feat);end
        %[~,sphcen0,sphrad]=ellipsoidalHoughTransform(feat,voxsiz,[10 25]);%This was before using radious in mmconversion to voxels
        [~,sphcen0,sphrad]=ellipsoidalHoughTransform(feat,voxsiz,[20 60]);
        rec.Par.Mine.sphcen0=sphcen0;
        rec.Par.Mine.sphrad=sphrad;

        %TORSO MASK
        %These are the parameters in use, probably this mask could be made a bit tighter, but haven't explored that yet
        %rec.Alg.parS.nErode=6;rec.Alg.parS.nDilate=16;
        Ml=refineMask(feat,rec.Alg.parS,voxsiz);

        %TRACKING
        if gpu;x=gpuArray(x);end
        xl=x;
        for s=1:N(4)
            %xl=dynInd(xl,s,4,ellipsoidalHoughTransform(dynInd(x,s,4),voxsiz,[10 25],3,[],Ml));%This returns the features
            xl=dynInd(xl,s,4,ellipsoidalHoughTransform(dynInd(x,s,4),voxsiz,[20 60],3,[],Ml));%This returns the features
        end;Ml=[];
        x=gather(x);
        NH=N.*(1+mirr);
        %NF=ones(1,length(N))/8;
        NF=voxsiz/32;%RECENT CHANGE
        %NF(4)=1/4;
        NF(4)=1/8;%RECENT CHANGE
        H=buildFilter(NH,'tukeyIso',NF,0,1,mirr);%It does not fit in the GPU, otherwise it should be done there
        if gpu;H=gpuArray(H);end
        xl=abs(filtering(xl,H,mirr));H=[];

        xl=reshape(xl,[prod(N(1:3)) N(4)]);
        [~,indMaccum]=max(xl,[],1);xl=[];indMaccum=gather(indMaccum);
        sphcen=ind2subV(N(1:3),indMaccum);%Centers after tracking
        sphcen=bsxfun(@minus,sphcen,sphcen0);%Shifts to be applied

        %SHIFTING THE DATA
        NShift=1;%DIRECTIONS FOR SHIFTING, ONLY READOUT
        T=cell(1,NShift);
        perm=1:4;perm([1 4])=[4 1];
        sphcenp=permute(sphcen,perm);
        for c=1:NShift
            T{c}=-dynInd(sphcenp,c,2);
            if gpu;T{c}=gpuArray(T{c});end
        end
        if gpu;x=gpuArray(x);end
        x=abs(shifting(x,T));

        %BRAIN MASK
        [x,feat]=parUnaFun({x,feat},@gather);
        feat=multDimMed(x,4:ND);
        if gpu;feat=gpuArray(feat);end
        K=32;
        S=128;
        rRange=[5 35];%RECENT CHANGE
        like=[0 -1 -1 0];
        pri=[0 0 0 0 5 1];%RECENT CHANGE
        [rec.M,par]=brainSegmentation(feat,sphcen0,rRange,K,S,sphrad,like,pri);feat=[];
        rec.Par.Mine.EllipsoidParameters(1,:)=par;
        rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,8);rec.Dyn.Typ2Wri(8)=1;
        if gpu;x=gpuArray(x);end

        %MASK USED TO COMPUTE THE ROI-FIFTH COMPONENT OF THE MASK IMAGE
        dilate=10;dilate=dilate*[1 1.2 1];
        rec.M=cat(4,rec.M,morphFourier(dynInd(rec.M,3,4),dilate,voxsiz,mirr(1:3)));

        %SOFT MASK-SIXTH COMPONENT OF THE MASK IMAGE
        Msoft=morphFourier(dynInd(rec.M,5,4),distan,voxsiz,mirr(1:3),1);%SOFT MASK USED TO COMPUTE THE ROI
        Msoft(Msoft>1-1e-6)=1;Msoft(Msoft<1e-6)=0;
        rec.M=cat(4,rec.M,Msoft);Msoft=[];

        %ROI COMPUTATION
        rec.Enc.ROI=computeROI(dynInd(rec.M,6,4));
        rec.Enc.ROI(3,:)=[1 N(3) N(3) N(3) 0 0];%To avoid problems later on with MB replication
        fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
        rec.Par.Mine.sphcen0=rec.Par.Mine.sphcen0-rec.Enc.ROI(1:3,1)';
        rec.Par.Mine.EllipsoidParameters(1,1:3)=rec.Par.Mine.EllipsoidParameters(1,1:3)-(rec.Enc.ROI(1:3,1)').*voxsiz;

        %ROI EXTRACTION FOR REGISTRATION
        x=extractROI(x,rec.Enc.ROI,1,1:3);
        M=extractROI(dynInd(rec.M,3,4),rec.Enc.ROI,1,1:3);
        rec.M=gather(rec.M);

        %REGISTRATION
        ND=numDims(x);ND=max(ND,4);
        rec.T=single(zeros([ones(1,ND-1) N(ND) 6]));
        NDT=numDims(rec.T);
        %lev=2;
        %parComp=[1 1 1 0 0 0];
        %[~,rec.T,~]=groupwiseVolumeRegistration(x,M,rec.T,[],lev,rec.Alg.parU.fractionOrder,parComp);x=[];
        kmax=[3 3];
        dk=[2 1];
        lev=[2 2];
        [~,rec.T,~]=integerShiftRegistration(x,M,rec.T,kmax,dk,lev,rec.Alg.parU.fractionOrder);x=[];M=[];

        %ADDING PREVIOUS LOCATION OF THE CENTERS
        perm=1:NDT;perm([1 2 ND NDT])=[ND NDT 1 2];
        sphcen=permute(sphcen,perm);
        rec.T=dynInd(rec.T,1:NShift,NDT,dynInd(rec.T,1:NShift,NDT)+dynInd(sphcen,1:NShift,NDT));
        rec.T=dynInd(rec.T,NShift+1:6,NDT,0);
        rec.T=gather(rec.T);

        %SHIFTING
        for c=1:3
            T{c}=-round(dynInd(rec.T,c,NDT));
            if gpu;T{c}=gpuArray(T{c});end
        end
        if gpu;rec.w=gpuArray(rec.w);end
        
        rec.b=shifting(rec.w,T);
        rec.w=gather(rec.w);
        rec.b=gather(rec.b);
        if isfield(rec,'E')
            if gpu;rec.E=gpuArray(rec.E);end
            rec.E=gather(abs(shifting(rec.E,T)));%KEY INSTRUCTION
        end
        rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,26);rec.Dyn.Typ2Wri(26)=1;

        %ROI EXTRACTION
        typ2Rec=rec.Dyn.Typ2Rec;
        for n=typ2Rec';datTyp=rec.Plan.Types{n};
            if ~ismember(n,24);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1:3);end
        end
        
        %FOR SAVING TRACKING INFORMATION
        rec.t=cellFun(T,@gather);
        rec.Dyn.Typ2Wri(25)=1;
    end
else    
    T=cellFun(rec.t,@gpuArray);
    rec.E=abs(shifting(rec.E,T));
    M=dynInd(rec.M,6,4);rec=rmfield(rec,'M');
    M(:)=1;M([1 end],:,:)=0;M(:,[1 end],:)=0;
    M=mapVolume(M,rec.E,rec.MT{2},rec.MT{1});
    M=single(M>0.5);
    rec.Enc.ROI=computeROI(M,1);    
    fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    rec.E=extractROI(rec.E,rec.Enc.ROI,1,1:3);
    rec.E=gather(rec.E);
end
