function rec=volumeAlignment(rec,appl)

% VOLUMEALIGNMENT performs temporal alignment of volumes
%   REC=VOLUMEALIGNMENT(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   ** REC is a reconstruction structure with information in undistorted 
%   coordinates
%

if nargin<2 || isempty(appl);appl=0;end

if rec.Fail || ~rec.Alg.parU.corrMotion;return;end

AdHocArray=rec.Par.Mine.AdHocArray;
assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
if AdHocArray(1)==101;MB=AdHocArray(4);else MB=1;end

gpu=rec.Dyn.GPU;
if ~appl
    NDims=numDims(rec.M);NDims=min(NDims,3);
    voxsiz=rec.Par.Scan.AcqVoxelSize(1,1:NDims);
    voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
    gpu=rec.Dyn.GPU;

    if ~rec.Alg.parU.useUndist || ~isfield(rec,'u')
        if rec.Par.Mine.Modal==10
            if rec.Alg.SVDRecover;rec.u=abs(rec.r);else rec.u=abs(rec.x);end
        else
            rec.u=abs(rec.b);
        end
    end
    if gpu;rec.u=gpuArray(rec.u);rec.T=gpuArray(rec.T);end
    if isfield(rec,'x');rec.x=gather(rec.x);end
    if isfield(rec,'b');rec.b=gather(rec.b);end
    if isfield(rec,'M') && gpu;rec.M=gpuArray(rec.M);end

    ND=numDims(rec.u);
    mirr=zeros(1,NDims);
    if rec.Par.Mine.Modal==9 && isfield(rec.Par.Mine,'sphcen0')
        %FOR FMRI WITH TRACKING WE SEGMENT AGAIN
        feat=multDimMed(rec.u,4:ND);
        K=32;
        S=128;    
        like=[0 -1 -1 0];    
        pri=[0 0 0 0 10 5];
        %rRange=[5 30];
        rRange=[5 35];%RECENT CHANGE
        %rRange=[max(min(rec.Par.Mine.EllipsoidParameters(1,4:6))-2,6) max(rec.Par.Mine.EllipsoidParameters(1,4:6))+10]/(prod(voxsiz).^(1/3));
        %fprintf('Radious range used for segmentation: [%.2f %.2f]\n',rRange(1),rRange(2));
        [M,par]=brainSegmentation(feat,rec.Par.Mine.sphcen0,rRange,K,S,rec.Par.Mine.sphrad,like,pri);
        %[M,par]=brainSegmentation(feat,rec.Par.Mine.EllipsoidParameters(1,1:3)./voxsiz,rRange,K,S,rec.Par.Mine.sphrad,like,pri);%RECENT CHANGE
        rec.Par.Mine.EllipsoidParameters(2,:)=par;
        rec.M=cat(4,rec.M,M);%7th-10th components    
        Tround=round(rec.T);
    elseif ~isfield(rec.Par.Mine,'sphcen0')
        rec.M=refineMask(multDimSum(abs(rec.u),4:5),rec.Alg.parS,voxsiz);    
        distan=3*ones(1,NDims);
        %rec.M=morphFourier(morphFourier(rec.M,-distan,voxsiz,mirr),distan,voxsiz,mirr,1);%Soft masking
        rec.M=morphFourier(rec.M,distan,voxsiz,mirr,1);%Soft masking
        if isfield(rec,'T');Tround=rec.T;end
    else
        distan=20*ones(1,NDims);
        rec.M=morphFourier(rec.M,distan,voxsiz,mirr,5);%Soft masking
        Tround=rec.T;
    end

    NE=length(rec.Par.Labels.TE);
    if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't
    N=size(rec.u);N(end+1:4)=1;
    if NE>1 && numDims(rec.u)<=4;rec.u=reshape(rec.u,[N(1:3) NE N(4)/NE]);end
    if NE==1 && numDims(rec.u)>=5;rec.u=reshape(rec.u,[N(1:3) prod(N(4:end))]);end

    %if numDims(rec.u)>4;fprintf('Motion correction not defined for %d dimensions',numDims(rec.u));return;end
    ND=numDims(rec.u);ND=max(ND,4);

    %%%%HEREHERE---LAMBDA HAS TO BE 10X FOR DWI, IT IS INSIDE, REMOVE AND
    %%%%CONTROL FROM OUTSIDE. ACCEL HAS TO BE DISABLED, FINAL RECONSTRUCTION
    %%%%HAS TO BE PERFORMED

    if rec.Alg.parU.corrMotion    
        N=size(rec.u);N(end+1:4)=1;
        if ~isfield(rec,'v')
            if rec.Par.Mine.Modal==10 || ~isfield(rec.Par.Mine,'sphcen0')
                if isfield(rec.Par.Mine,'sphcen0');y=updateMask(rec.u);end
                rec.T=single(zeros([ones(1,ND-1) N(ND) 6]));
                M=rec.M;
            else
                rec.T=rec.T-Tround;
                M=softenMask(dynInd(rec.M,9,4),6);
            end        
            if ~rec.Alg.parU.iterMask;lev={[2 1]};else lev={2,[2 1]};end     
            for n=1:1+rec.Alg.parU.iterMask%Iterate between volume registration and mask refinement   
                if ~isfield(rec,'E') || n~=1+rec.Alg.parU.iterMask;[rec.v,rec.T,rec.d]=groupwiseVolumeRegistration(rec.u,M,rec.T,[],lev{n},rec.Alg.parU.fractionOrder);
                else [rec.v,rec.T,rec.d,rec.E]=groupwiseVolumeRegistration(rec.u,M,rec.T,[],lev{n},rec.Alg.parU.fractionOrder,[],rec.E);rec.E=abs(rec.E);
                end
                rec.v=abs(rec.v);
                if rec.Alg.parU.iterMask
                    if rec.Par.Mine.Modal==10 || ~isfield(rec.Par.Mine,'sphcen0')
                        if isfield(rec.Par.Mine,'sphcen0');y=updateMask(rec.v);end
                        M=rec.M;
                    else
                        M=softenMask(dynInd(rec.M,9,4),5);
                    end
                end
            end

            rec.v=gather(rec.v);rec.e=rec.v;
            rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,[14;16]);rec.Dyn.Typ2Wri([14 16])=1;

            if rec.Par.Mine.Modal==10 || ~isfield(rec.Par.Mine,'sphcen0')
                M=rec.M;
                y=rec.u;
            else
                %FINAL SEGMENTATION                   
                feat=multDimMed(rec.v,4:ND);
                K=32;
                S=128;
                %rRange=[5 30];
                rRange=[5 35];%RECENT CHANGE
                %rRange=[max(min(rec.Par.Mine.EllipsoidParameters(2,4:6))-2,6) max(rec.Par.Mine.EllipsoidParameters(2,4:6))+10]/(prod(voxsiz).^(1/3));%RECENT CHANGE
                %fprintf('Radious range used for segmentation: [%.2f %.2f]\n',rRange(1),rRange(2));
                like=[0 -1 -1 0];    
                pri=[0 0 0 0 10 5];
                [M,par]=brainSegmentation(feat,rec.Par.Mine.sphcen0,rRange,K,S,rec.Par.Mine.sphrad,like,pri);
                %[M,par]=brainSegmentation(feat,rec.Par.Mine.rec.Par.Mine.EllipsoidParameters(2,:)./voxsiz,rRange,K,S,rec.Par.Mine.sphrad,like,pri);%RECENT CHANGE
                rec.Par.Mine.EllipsoidParameters(3,:)=par;
                rec.M=cat(4,rec.M,M);%11th-14th components
                M=softenMask(dynInd(rec.M,13,4),4);
                y=rec.u;
            end

            if rec.Alg.parU.corrMotion==2
                rec.T=repmat(rec.T,[ones(1,ND-1) N(3)/MB 1]);           
                rec.T=reshape(rec.T,[ones(1,ND-1) N(ND) N(3)/MB 6]);
                if rec.Par.Mine.Modal==9
                     Tround=repmat(Tround,[ones(1,ND-1) N(3)/MB 1]);
                     Tround=reshape(Tround,[ones(1,ND-1) N(ND) N(3)/MB 6]);
                end         
                %TT=rec.T;
                %parU=rec.Alg.parU;
                %save('/home/lcg13/Work/DataDefinitiveImplementationDebug06/groupSliceYeShearAlt.mat','y','M','TT','parU','MB','Tround');
                %1
                %pause
                if ~isfield(rec,'E');[rec.e,rec.T,rec.d]=groupwiseSliceRegistration(y,M,rec.T,rec.Alg.parU.Lambda,MB,[],rec.Alg.parU.fractionOrder,[],rec.Alg.parU.interpolTime(2));
                else [rec.e,rec.T,rec.d,rec.E]=groupwiseSliceRegistration(y,M,rec.T,rec.Alg.parU.Lambda,MB,[],rec.Alg.parU.fractionOrder,rec.E,rec.Alg.parU.interpolTime(2));rec.E=abs(rec.E);
                end
                [rec.e,rec.d]=parUnaFun({rec.e,rec.d},@abs);
                [rec.e,rec.d]=parUnaFun({rec.e,rec.d},@gather);
                rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,15);rec.Dyn.Typ2Wri(15)=1;
            end
            if rec.Alg.parU.writeVideo && rec.Par.Mine.Modal==9;writeVideo(rec);end
        else%BE CAREFUL THIS WON'T GIVE EXACTLY THE SAME RESULTS AS RUNNING BOTH CORRECTIONS SEQUENTIALLY
            if rec.Par.Mine.Modal==10
                updateMask(rec.v);
                rec.T=single(zeros([ones(1,ND-1) N(ND)*N(3)/MB 6]));
                M=rec.M;
            else
                rec.T=rec.T-Tround;
                rec.T=repmat(rec.T,[ones(1,ND-1) N(3)/MB 1]);
                Tround=repmat(Tround,[ones(1,ND-1) N(3)/MB 1]);
                if isfield(rec.Par.Mine,'sphcen0');M=softenMask(dynInd(rec.M,9,4),4);else M=rec.M;end
            end
            if rec.Alg.parU.corrMotion==2            
                rec.T=reshape(rec.T,[ones(1,ND-1) N(ND) N(3)/MB 6]);
                if rec.Par.Mine.Modal==9
                    Tround=reshape(Tround,[ones(1,ND-1) N(ND) N(3)/MB 6]);
                end
                if ~isfield(rec,'E');[rec.e,rec.T,rec.d]=groupwiseSliceRegistration(rec.v,M,rec.T,rec.Alg.parU.Lambda,MB,[],rec.Alg.parU.fractionOrder,[],rec.Alg.parU.interpolTime(2));
                else [rec.e,rec.T,rec.d,rec.E]=groupwiseSliceRegistration(rec.v,M,rec.T,rec.Alg.parU.Lambda,MB,[],rec.Alg.parU.fractionOrder,rec.E,rec.Alg.parU.interpolTime(2));rec.E=abs(rec.E);
                end
                [rec.e,rec.d]=parUnaFun({rec.e,rec.d},@abs);
                rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,[15;16]);rec.Dyn.Typ2Wri(15:16)=1;
            end
        end
        rec.Dyn.Typ2Wri(17)=1;                
        if rec.Par.Mine.Modal==9
            %SPIN HISTORY
            TR=rec.Par.Labels.RepetitionTime(1)/1000;%in seconds
            FA=convertRotation(rec.Par.Labels.FlipAngle(1),'deg','rad');%in radians
            T1=2.5;%2.75 is mean value%Around lowest values for white matter at 3T in Williams et al., Neonatal Brain: Regional Variability of in Vivo MR Imaging Relaxation Rates at 3.0 T?Initial Experience. Radiology, 2005.
            rec.H=spinHistoryEstimation(y,rec.T,T1,TR,MB,FA,3);
            rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,30);rec.Dyn.Typ2Wri(30)=1;
            rec.T=rec.T+Tround;
        end    
    end
else
    NT=size(rec.T);ND=length(NT);
    Tround=zeros(NT,'like',rec.T);
    for c=1:3;Tround=dynInd(Tround,c,ND,repmat(-round(gather(rec.t{c})),[ones(1,ND-2) NT(ND-1)]));end
    rec.T=rec.T-Tround;
    if gpu;rec.E=gpuArray(rec.E);rec.T=gpuArray(rec.T);end
    [~,~,~,rec.E]=groupwiseSliceRegistration([],[],rec.T,rec.Alg.parU.Lambda,MB,[],rec.Alg.parU.fractionOrder,rec.E,rec.Alg.parU.interpolTime(2));rec.E=abs(rec.E);    
end

function y=updateMask(x)
    feat=bsxfun(@times,rec.M,multDimMed(x,4:ND));
    rec.M=sphericalHoughTransform(feat,[10 25]);
    distan=5*ones(1,NDims);
    rec.M=morphFourier(rec.M,distan,voxsiz,mirr,1);%Soft masking
    rec.M(rec.M>1-1e-6)=1;rec.M(rec.M<1e-6)=0;
    y=rec.u;
end

function M=softenMask(M,sof)   
    dilate=sof*ones(1,3);%It was 4
    M=morphFourier(M,dilate,voxsiz,mirr);
    M(M>1-1e-6)=1;M(M<1e-6)=0; 
    distan=sof*ones(1,3);%It was 4
    M=morphFourier(M,distan,voxsiz,mirr,1);%Soft masking%It was 2
    M(M>1-1e-6)=1;M(M<1e-6)=0;   
end

end