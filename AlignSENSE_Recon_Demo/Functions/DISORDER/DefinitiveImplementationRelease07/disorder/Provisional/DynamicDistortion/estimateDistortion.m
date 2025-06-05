function rec=estimateDistortion(rec)

% ESTIMATEDISTORTION estimates the B0 field from the phase information of a
% gradient echo or hybrid spin echo / gradient echo sequence
%   REC=ESTIMATEDISTORTION(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   ** REC is a reconstruction structure with B0 information (rec.B) 
%

if rec.Fail || ~rec.Alg.parU.useUndist;return;end%This return has been recently introduced to save memory

%Note there are problems in 2014_09_10/SU_18303/ZZ-HE/su_10092014_1426506_14_2_dev1t2stest2mmsenseV4.mat due to 3 echos

AdHocArray=rec.Par.Mine.AdHocArray;
assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
if AdHocArray(1)==101;MB=AdHocArray(4);else MB=1;end

%NUMBER OF ECHOES AND CHECK THAT SAFE METHOD IS BEING USED
NE=length(rec.Par.Labels.TE);
if rec.Par.Labels.TE(end)<rec.Par.Labels.TE(1);NE=1;end%Sometimes it appears there is a second TE while there isn't

%DUAL ECHO DIFFERENCE
if ~rec.Alg.SVDRecover && rec.Par.Mine.Modal==10     
    typ2Rec=rec.Dyn.Typ2Rec;
    for n=typ2Rec';datTyp=rec.Plan.Types{n};%FINISH ROI EXTRACTION
        rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,2);
    end    
    N=size(rec.x);N(end+1:5)=1;
    N(4:5)=[NE N(4)/NE];
    rec.x=reshape(rec.x,N);
end
if NE==1 && rec.Par.Mine.Modal==10;return;end

if ~isfield(rec,'u')
    if rec.Par.Mine.Modal==10
        if ~rec.Alg.SVDRecover;rec.u=abs(rec.x);else rec.u=abs(rec.r);end
    else
        rec.u=abs(rec.b);
    end       
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,13);
    rec.Dyn.Typ2Wri(13)=rec.Alg.parU.useUndist;
end
gpu=rec.Dyn.GPU;
rec.u=gather(rec.u);
if NE>2 && rec.Par.Mine.Modal==10;rec.Par.Mine.Spin2Echo=0;end
if NE>1 && rec.Alg.parU.useUndist && rec.Par.Mine.Modal==10
    if isfield(rec.Par.Mine,'Spin2Echo')
        if ~rec.Alg.SVDRecover
            if ~rec.Par.Mine.Spin2Echo;rec.x=dynInd(rec.x,2,4).*conj(dynInd(rec.x,1,4));else rec.x=dynInd(rec.x,2,4).*dynInd(rec.x,1,4);end
        else
            if ~rec.Par.Mine.Spin2Echo;rec.x=dynInd(rec.r,2,4).*conj(dynInd(rec.r,1,4));else rec.x=dynInd(rec.r,2,4).*dynInd(rec.r,1,4);end
            rec.r=gather(rec.r);
        end
    end
end

%PERMUTE
typ2Rec=rec.Dyn.Typ2Rec;
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    rec.(datTyp)=permute(rec.(datTyp),[1 2 3 5 4]);
end

if rec.Par.Mine.Modal==10;N=size(rec.x);N(end+1:5)=1;else N=size(rec.u);N(end+1:5)=1;end
%ND=numDims(rec.x);
% if rec.Par.Mine.Modal==9
%     %REFINE MASKING USING SECOND ORDER PHASE FINITE DIFFERENCES    
%     q2d=finiteDifferenceComplex(rec.x,3,2)/pi;
%     rad=1;
%     BlSz=100;
%     for s=1:BlSz(1):N(4);vS=s:min(s+BlSz(1)-1,N(4));    
%         q2d=dynInd(q2d,vS,4,(dynInd(q2d,vS,4)-morphFourier(dynInd(q2d,vS,4),rad,[1 1 1],[0 0 0],1)).^2);     
%         q2d=dynInd(q2d,vS,4,sqrt(abs(morphFourier(dynInd(q2d,vS,4),rad,[1 1 1],[0 0 0],1)))); 
%     end
%     q2d=1-mean(abs(q2d),ND+1);
%     q2dv=mean(q2d,4);
% 
%     %if ~rec.Alg.parU.useHough
%         parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
%         parS.maskTh=0.2;%Threshold of the body coil intensities for mask extraction%It was 0.2
%         parS.nErode=6;%12;%16;%Erosion for masking (in mm)
%         parS.nDilate=20;%28*[1 1 1];%Dilation for masking (in mm)
%         parS.conComp=1;%1;
%         parS.Otsu=[0 1];%1;
%         rec.M=refineMask(bsxfun(@times,rec.M,q2dv),parS,rec.Enc.AcqVoxelSize);
%         distan=10*ones(1,ndims(rec.M));
%         rec.M=morphFourier(rec.M,distan,rec.Enc.AcqVoxelSize,[1 1 1],1);%Soft masking
%         rec.M(rec.M>1-1e-6)=1;rec.M(rec.M<1e-6)=0;
%     %end
% end

if rec.Par.Mine.Modal==9 && isfield(rec.Par.Mine,'sphcen0')
    Ms=dynInd(rec.M,6,5);%Full soft mask
else
    Ms=rec.M;
end

voxsiz=rec.Par.Scan.AcqVoxelSize(1,1:3);
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
voxsiz(3)=voxsiz(3)/rec.Alg.parU.weightZ;
ND=max(numDims(rec.x),4);
voxsiz(ND)=voxsiz(1)/rec.Alg.parU.weightT;
voxsiz(voxsiz==0)=1;%We assume there are the echos here

if ~isfield(rec,'B') && rec.Alg.parU.useUndist   
    if rec.Par.Mine.Modal==10
        if strcmp(rec.Alg.parU.UnwrapMeth,'PUMA');rec.B=PUMAUnwrapping(rec.x,Ms,rec.Par.Mine.Modal);
        %elseif strcmp(rec.Alg.parU.UnwrapMeth,'LUCK');rec.B=LUCKUnwrapping(rec.x,rec.M,rec.Par.Mine.Modal);        
        elseif strcmp(rec.Alg.parU.UnwrapMeth,'CNCG')
            rb=rec.x;
            rec.x=gather(rec.x);
            if rec.Alg.parU.interpolTime(1);rb=shear(rb,1,MB);end
            rb=CNCGUnwrapping(rb,rec.Alg.parU.LUnwrap,voxsiz,rec.Alg.parU.Upsampling,rec.Alg.parU.weightFunc);
            if rec.Alg.parU.interpolTime(1);rb=shear(rb,0,MB);end
            rec.B=rb;rb=[];           
        else fprintf('Unknown unwrapping method %s\n',rec.Alg.parU.UnwrapMeth);rec.Fail=1;return;
        end
        rec.x=gather(rec.x);
    else        
        NT=size(rec.b,5);        
        if NT>1000
            NV=round(NT/10);
            wi=circconvmtx(ones(1,2*NV)',NT);          
            wi=wi(:,1:NV:end);
            NW=size(wi,2);
            NB=size(rec.b);
            B=zeros([NB(1:3) 1 2*NV NW],'single');
            for n=1:NW;indw=find(wi(:,n)==1);
                rb=dynInd(rec.b,indw,5);
                if rec.Alg.parU.interpolTime(1);rb=shear(rb,1,MB);end
                rb=CNCGUnwrapping(rb,rec.Alg.parU.LUnwrap,voxsiz,rec.Alg.parU.Upsampling,rec.Alg.parU.weightFunc);
                if rec.Alg.parU.interpolTime(1);rb=shear(rb,0,MB);end
                B=dynInd(B,n,6,gather(rb));
            end
            medB=resPop(B,1:5,prod([NB(1:3) 1 2*NV]),1);
            medB=median(medB,1);
            B=bsxfun(@minus,B,medB)+medB(1);
            rec.B=gather(rec.b);rec.B(:)=0;
            for n=1:NW;indw=find(wi(:,n)==1);
                rec.B=dynInd(rec.B,indw,5,dynInd(B,n,6)+dynInd(rec.B,indw,5));
            end
            B=[];            
            rec.B=bsxfun(@rdivide,rec.B,permute(sum(wi,2),[2 3 4 5 1]));
        else
            %if numel(rec.B)>1e8;rec.B=gather(rec.B);end
            if strcmp(rec.Alg.parU.UnwrapMeth,'PUMA');rec.B=PUMAUnwrapping(rec.b,Ms,rec.Par.Mine.Modal);
            %elseif strcmp(rec.Alg.parU.UnwrapMeth,'LUCK');rec.B=LUCKUnwrapping(rec.b,rec.M,rec.Par.Mine.Modal);
            elseif strcmp(rec.Alg.parU.UnwrapMeth,'CNCG')
                rb=rec.b;
                rec.b=gather(rec.b);
                if rec.Alg.parU.interpolTime(1);rb=shear(rb,1,MB);end
                rb=CNCGUnwrapping(rb,rec.Alg.parU.LUnwrap,voxsiz,rec.Alg.parU.Upsampling,rec.Alg.parU.weightFunc);
                if rec.Alg.parU.interpolTime(1);rb=shear(rb,0,MB);end
                rec.B=rb;rb=[];                
            else fprintf('Unknown unwrapping method %s\n',rec.Alg.parU.UnwrapMeth);rec.Fail=1;return;
            end            
        end
        rec.b=gather(rec.b);
        if gpu;rec.B=gpuArray(rec.B);end
    end
    %DISAMBIGUATE REFERENCE PHASE
    z0=resPop(rec.B,1:3,prod(N(1:3)),1);
%     if rec.Par.Mine.Modal==9--ALTERNATIVE TO DISAMBIGUATE
%         M=dynInd(rec.M,3,4);%Brain mask
%     else
%         M=rec.M;
%     end
    z0=dynInd(z0,Ms~=0,1);    
    %if rec.Par.Mine.Modal==10;z0=median(resPop(z0,1:5,numel(z0),1));else z0=median(z0,1);end
    z0=median(z0(:));%CHANGED
    z0=2*pi*round(z0/(2*pi));%To round to the nearest 2Npi value
    rec.B=bsxfun(@minus,rec.B,z0);
    rec.B=bsxfun(@times,rec.B,Ms);

    %for n=typ2Rec';datTyp=rec.Plan.Types{n};
    %    if n~=12;rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0);end   
    %end
    rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,11);rec.Dyn.Typ2Wri(11)=1;
end
if gpu;rec.u=gpuArray(rec.u);end
if rec.Alg.parU.useUndist==2;rec.B=dynInd(rec.B,1,5);end
%PERMUTE BACK
typ2Rec=rec.Dyn.Typ2Rec;
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    rec.(datTyp)=permute(rec.(datTyp),[1 2 3 5 4]);
end
rec.u=gather(rec.u);
rec.B=gather(rec.B);

end

function x=shear(x,di,MB)%BETTER NOT TO USE DUE TO JUMPS FROM FREQUENCY STABILIZATION PROBABLY
    if MB>1
        %if ~di;x=circshift(x,1,5);end
        
        %NX=size(x);
        %y=dynInd(x,floor(NX(2)/2),2);
        %y=permute(y,[1 3 5 2 4]);
        %y=reshape(y,[NX(1) NX(3) 4 5]);
        %y=permute(y,[1 3 2 4]);
        %y=reshape(y,[NX(1)*4 NX(3)*5]);
        %visReconstruction(y,0)
        NRB=size(x);NRB(end+1:5)=1;        
        x=reshape(x,[NRB(1:2) NRB(3)/MB MB NRB(4:end)]);
        if di;x=shearing(x,1,[6 4]);else x=shearing(x,-1,[6 4]);end%THIS SEEMS TO GIVE BEST RESULTS, BUT I'M NOT SURE WHY, CHECK!
        x=reshape(x,NRB);
        %y=dynInd(x,floor(NX(2)/2),2);
        %y=permute(y,[1 3 5 2 4]);
        %y=reshape(y,[NX(1) NX(3) 4 5]);
        %y=permute(y,[1 3 2 4]);
        %y=reshape(y,[NX(1)*4 NX(3)*5]);
        %visReconstruction(y,1)   
        
         %if di;x=circshift(x,-1,5);end
    end
end
