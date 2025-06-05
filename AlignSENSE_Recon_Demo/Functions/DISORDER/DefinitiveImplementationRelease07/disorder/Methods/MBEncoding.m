function [E,EH,fail]=MBEncoding(rec,pe)

%MBENCODING   Generates the encoding fields necessary to perform a MB 
%reconstruction
%   [E,EH]=MBENCODING(REC,PE)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   * PE is the phase encoding index (for multiple phase encodes)
%   * E is an encoding structure with MB fields
%   * EH is a decoding structure with MB fields
%   * FAIL indicates a failure to reconstruct
%

fail=0;
E=[];EH=[];
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU && rec.Dyn.GPU~=0);
AdHocArray=rec.Par.Mine.AdHocArray;

if isempty(AdHocArray) || (rec.Enc.MB==1 && AdHocArray(7)<=1) || all(AdHocArray==0:127);MBData=0;else MBData=1;end

if MBData && AdHocArray(3)==3%TO DO---USE A GENERATEGRID FOR THIS
    sliceGap=AdHocArray(2)/rec.Enc.SlicDist;
    rampSl=sliceGap*AdHocArray(7);%Number of slices on which 2pi cycle happens.    
    ramp=2*pi/rampSl;%Ramp per slice    
    Isoc=rec.Par.Mine.Isoc(3,1,pe);%-1/3;    
    if isfield(rec.Par.Mine,'ShiftIsoc');shiftIsoc=round(rec.Par.Mine.ShiftIsoc(3,1,pe));end%To avoid ringing better to round
    %THIS HAS BEEN REMOVED RECENTLY
    %if length(rec.Par.Mine.pedsUn)>1;displPhase=single(ramp*(rec.Enc.ExcGrid-(dynInd(Isoc,rec.Par.Mine.Nat==rec.Par.Mine.pedsUn,3)-1)));else displPhase=single(ramp*(rec.Enc.ExcGrid-(dynInd(Isoc,1,3)-1)));end
    %if length(rec.Par.Mine.pedsUn)>1;displPhase=single(ramp*(rec.Enc.ExcGrid-(dynInd(Isoc,rec.Par.Mine.Nat==rec.Par.Mine.pedsUn,3)-1)));else displPhase=single(ramp*(rec.Enc.ExcGrid-(dynInd(Isoc,rec.Par.Mine.pedsUn,3)-1)));end    
    displPhase=single(ramp*(rec.Enc.ExcGrid-(Isoc-1)));
    if isfield(rec.Par.Mine,'ShiftIsoc') && isfield(rec.Par.Mine,'splitInfo')
        displPhase=cat(12,displPhase,single(ramp*(rec.Enc.ExcGrid-(Isoc+rec.Par.Mine.ShiftIsoc(3,1,pe)-1))));
    end
else
    displPhase=single(zeros(1,rec.Enc.FOVSize(3)));
end
perm=1:12;perm([2 3])=[3 2];
displPhase=permute(displPhase,perm);

%Z-BLIP ENCODING PROFILES
if abs(AdHocArray(7))<=1 || ~MBData;maux=0;else maux=0:abs(AdHocArray(7))-1;end%Shift pattern
E.Bl=length(maux);
if AdHocArray(8)~=0 || ~MBData;maux=maux-maux(end)/2;end%Balancing
if AdHocArray(8)==11;maux=circshift(flip(circconvmtx(maux,E.Bl)',1),[1 0]);end%Self-calibration
NgPE=rec.Enc.AcqSize(2);
kFullGrid=generateGrid(NgPE,gpu,NgPE,ceil((NgPE+1)/2));


%assert(size(unique(rec.Enc.kRange{2},'rows'),1)==1,'The encoding is different for different echoes, which is not supported by MB reconstructor');

%ATTEMPT TO ADD +mod(NgPE,2) INSTEAD OF 1 HAS NOT WORKED WELL SO HAD TO REVERT...
kGrid=(rec.Enc.kRange{2}(1,1):rec.Enc.kRange{2}(1,2))-kFullGrid{1}(1)+1;
Nblips=length(kGrid);
kGridor=kGrid(1);
kGrid=horzcat(1:kGridor-1,kGrid);

if MBData || rec.Alg.GhosCorRec
    NTotalblips=length(kGrid);
    repm=ones(1,12);repm([2:3 12])=[NgPE rec.Enc.FOVSize(3) size(displPhase,12)];
    B=single(ones(repm));repm(:)=1;
    if gpu;B=gpuArray(B);end
    if AdHocArray(8)==11;B=repmat(B,[1 1 1 1 1 size(rec.y,6)]);end
    for b=1:size(B,6)
        for m=1:E.Bl
            displPhaseFin=exp(1i*maux(mod(b-1,E.Bl)+1,m)*displPhase);
            %vect=m:length(maux):Nblips;
            vect=mod((m+kGridor-1-1),E.Bl)+1:size(maux,2):NTotalblips;
            repm(2)=length(vect);
            displPhaseFin=repmat(displPhaseFin,repm);
            B=dynInd(B,{kGrid(vect),b},[2 6],single(displPhaseFin));    
        end
    end

    %Y-ENCODING TRANSFORM
    [~,dkGrid]=sincKernel(kFullGrid{1}(:)/max(rec.Enc.SamRate(2),1),kFullGrid{1}(:)/max(rec.Enc.SamRate(2),1),gpu);
    [E.Ef,EH.Eb]=buildDFTM(rec.Enc.rGrid{pe}{2}(:),kFullGrid{1}(:),[],dkGrid(:));
end

if MBData
    if AdHocArray(3)==1;B=dynInd(B,1:rec.Enc.FOVSize(3)/2,3,1);end

    if gpuIn;[E.Ef,EH.Eb,B]=parUnaFun({E.Ef,EH.Eb,B},@gpuArray);end
    NB=size(B);NB(end+1:rec.Plan.NDims)=1;    
    B=reshape(B,[NB(1:2) NB(3)/rec.Enc.MB rec.Enc.MB 1 NB(6:rec.Plan.NDims)]);
    B=applyMBPhases(B);
    B=reshape(B,NB);
    E.Bf=B;EH.Bb=conj(B);
end

%GHOSTING INFORMATION
if rec.Alg.GhosCorRef==3 || rec.Alg.GhosCorRec
    rec.Corr.z{5}=dynInd(rec.Corr.z{5},1,7);%We assume same alternation: TODO-> generalize (wouldn't be a problem but for back-compatibility with first SAFE tests...)
    if size(rec.Corr.z{5},rec.Plan.NDims+2)>1         
        perm=1:rec.Plan.NDims+2;perm([2 rec.Plan.NDims+2])=[rec.Plan.NDims+2 2];
        rec.Corr.z{5}=permute(rec.Corr.z{5},perm);
        rec.Corr.z{5}=dynInd(rec.Corr.z{5},1:Nblips,2);
    end
    E.Gh.Mf=rec.Corr.z{5};
    E.Gh.nD=unique(rec.Corr.z{5}(:))';    
    Nadd=NgPE-size(E.Gh.Mf,2);    
    E.Gh.Mf=horzcat(E.Gh.Mf(2:Nadd+1),E.Gh.Mf);%We assume partial Fourier: TODO-> generalize    
    EH.Gh.Mb=E.Gh.Mf;    
    EH.Gh.nD=E.Gh.nD;
    if ~isfield(rec.Par.Mine,'ShiftIsoc')
        if gpuIn;rec.Corr.P{2}=gpuArray(rec.Corr.P{2});end
        x=exp(1i*angle(dynInd(rec.Corr.P{2},pe,11)));
    else
        if gpuIn;rec.Corr.PExtr{2}=gpuArray(rec.Corr.PExtr{2});end
        x=exp(1i*angle(dynInd(rec.Corr.PExtr{2},pe,11)));
        H{8}=shiftIsoc;
        if isfield(rec.Par.Mine,'splitInfo');x=cat(12,x,shifting(x,H));else x=shifting(x,H);end
    end 
    perm=1:max(ndims(x),8);perm(8)=3;perm(3)=8;
    x=permute(x,perm);
    N=size(x);
    if rec.Enc.FOVSize(3)>size(x,3);fprintf('z-FOV of MB data (%d slices) is larger than for SB data (%d slices)\n',rec.Enc.FOVSize(3),size(x,3));fail=1;end
    N(3)=rec.Enc.FOVSize(3);   
    x=ifftshift(x,3);
    x=resampling(x,N,1);
    x=fftshift(x,3);
    E.Gh.Pf=x;
    if rec.Alg.GhosCorRef~=3;E.Gh.Pf(:)=1;end
    EH.Gh.Pb=conj(x);
end

function B=applyMBPhases(B)
    if rec.Enc.MB>=4
        mbPhases=[  0         0         0         0         0         0         0         0         0         0         0         0;  
                    0         0         0         0         0         0         0         0         0         0         0         0;
                    0    3.1416         0         0         0         0         0         0         0         0         0         0;
                    0    3.1416    3.1416         0         0         0         0         0         0         0         0         0;
                    0         0    3.1416         0         0         0         0         0         0         0         0         0;
               1.6910    2.8120    1.1570   -1.1570   -2.8120   -1.6910         0         0         0         0         0         0;
               2.5820   -0.5620    0.1020         0   -0.1020    0.5620   -2.5820         0         0         0         0         0;
               2.1120    0.2200    1.4640    1.9920   -1.9920   -1.4640   -0.2200   -2.1120         0         0         0         0;
               0.4790   -2.6670   -0.6460   -0.4190         0    0.4190    0.6460    2.6670   -0.4790         0         0         0;
               1.6830   -2.3950    2.9130    0.3040    0.7370   -0.7370   -0.3040   -2.9130    2.3950   -1.6830         0         0;
               1.4050    0.8870   -1.8540    0.0700   -1.4940         0    1.4940   -0.0700    1.8540   -0.8870   -1.4050         0;
               1.7290    0.4440    0.7220    2.1900   -2.1960    0.9840   -0.9840    2.1960   -2.1900   -0.7220   -0.4440   -1.7290]; 
           mbPhases=mbPhases(rec.Enc.MB,1:rec.Enc.MB);
           mbPhases=permute(mbPhases,[1 3 4 2]);
           mbPhases=single(exp(1i*mbPhases));
           if gpu;mbPhases=gpuArray(mbPhases);end
           B=bsxfun(@times,B,mbPhases);
    end
end

end