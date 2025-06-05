function Par=computeRAF(MR,Par,typ,red)

%COMPUTERAF   Obtains the transform from acquired grid to volume FOV as
%well as the isocenter in acquired grid
%   PAR=COMPUTERAF(MR,PAR,{TYP},{RED})
%   * MR is a MRecon object where to obtain the transform
%   * PAR contains current parameter information which is to be filled with
%   the matrix returning the coordinate transform in homogeneous
%   coordinates (MINE.APHIACQ), the isocenter (MINE.ISOC), the native
%   orientation (MINE.NAT), the offcentres in pixels (SCAN.MPSOFFCENTRES)
%   --not strictly necessary-- and mm (SCAN.MPSOFFCENTRESMM), the range in 
%   the readout (ENCODING.XRANGE) --not strictly necessary--, the range in the 
%   PE (ENCODING.YRANGE) --not strictly necessary--, the matrix for 
%   conversion of reconstructed FOVs (MINE.SIGNS), and the matrix returning
%   the coordinate transform in homogeneous coordinates for writing the
%   data (MINE.APHIREC)
%   * {TYP} is the type of transform to be obtained ('acq' for the
%   transform at acquired resolution, default, and 'rec' for the transform
%   at prescribed reconstructed resolution
%   * {RED} is a flag that indicates whether to use the minimum of acquired
%   resolutions in the PE direction, defaults to 0
%   ** PAR is the filled specific parameter information
%   ** FAIL indicates failures
%
%TODO: THIS CODE COMES FROM RELEASE04 BUT IT SHOULD BE REMOVED IN THE
%FUTURE. IN THIS RELEASE WE ONLY USE THE TRANSFORM FROM ACQUIRED GRID TO
%ACQUIRED FOV, NOT THE TRANSFORM TO PRESCRIBED RECONSTRUCTION RESOLUTION
%

if nargin<3 || isempty(typ);typ='acq';end
if nargin<4 || isempty(red);red=0;end

fail=0;
%if isfield(MR.Parameter.Labels,'AdHocArray');AdHocArray=MR.Parameter.Labels.AdHocArray;elseif MR.Parameter.IsParameter('MP_RECFRAME_ad_hoc');AdHocArray=MR.Parameter.GetValue('MP_RECFRAME_ad_hoc');else AdHocArray=[];end
AdHocArray=Par.Mine.AdHocArray;
if ~isempty(AdHocArray)
    if AdHocArray(1)==102 && AdHocArray(22)==0;SIR=2;else SIR=1;end
    if AdHocArray(1)==101 || AdHocArray(1)==102;MultiBand=AdHocArray(4);else MultiBand=1;end
else
    MultiBand=1;SIR=1;
end
AccelSlice=MultiBand*SIR;

if size(MR.Parameter.Scan.RecVoxelSize,2)==3;RecVoxelSize=MR.Parameter.Scan.RecVoxelSize(1,:);elseif size(MR.Parameter.Scan.RecVoxelSize,2)==6;RecVoxelSize=MR.Parameter.Scan.RecVoxelSize(1,1:2:end);end
AcqVoxelSize=MR.Parameter.Scan.AcqVoxelSize(1,:);
AcqVoxelSizeRec=AcqVoxelSize;


Par.Mine.UseReadRes=0;%To use the maximum resolution in plane
if Par.Mine.Modal==10 && Par.Labels.PartialFourierFactors(1,2)~=1;Par.Mine.UseReadRes=1;end
Par.Mine.AddedSamples=zeros(1,2);

%%%HEREHEREHERE---ADD PAR TO CALL THE FUNCTION

if ~isempty(Par.Mine.diInfo)
    if strcmp(MR.Parameter.Scan.MPS,'PA RL HF');Par.Mine.Nat=1;%LR
    elseif strcmp(MR.Parameter.Scan.MPS,'AP LR HF');Par.Mine.Nat=2;%RL
    elseif strcmp(MR.Parameter.Scan.MPS,'LR PA HF');Par.Mine.Nat=3;%AP
    elseif strcmp(MR.Parameter.Scan.MPS,'RL AP HF');Par.Mine.Nat=4;end%PA
    Par.Mine.pedsUn=unique(Par.Mine.diInfo(:,5));
    
    %%%%HEREHEREHERE---I'VE MODIFIED THIS BUT NEVER TESTED ON ADULT DATA
    %if numel(Par.Mine.pedsUn)==1 %&& strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') && strcmp(rec.Par.Mine.curStud.IssuId,'DiscardDWI');
    %    Par.Mine.pedsUn=Par.Mine.Nat;Par.Mine.diInfo(:,5)=Par.Mine.Nat;
    %end
    %if numel(Par.Mine.pedsUn)<=2;Par.Mine.UseReadRes=0;end
    if Par.Mine.UseReadRes
        AcqVoxelSizeRec(1:2)=min(AcqVoxelSize(1:2))*ones(1:2);
        modifyRes;
        if fail;Par=[];return;end
    end
    
    %INFORMATION FOR FAKED PARAMETER GENERATION
    iFoldOverDir={'RL','RL','AP','AP'};
    iFatShiftDir={'L','R','A','P'};
    ijkMPS{1}={'PA','RL','HF'};ijkMPS{2}={'AP','LR','HF'};ijkMPS{3}={'LR','PA','HF'};ijkMPS{4}={'RL','AP','HF'};
    scanFOV=MR.Parameter.Scan.FOV;
    xyzOffcentres=MR.Parameter.Scan.xyzOffcentres;
    xx=MR.Transform([0 0 0],'MPS','ijk')';
    NoX=round(xx(1)-0.5)*2;
    NoY=round((xx(2)-0.5)*2/MR.Parameter.Labels.SENSEFactor(2));%This has given small problems; hope should not need to be modified
    A=MR.Transform('xyz','rec')';A=A(1:2,1:2);A(abs(A)<1e-6)=0;A=sign(A);%This may present problems in case of oblique acquisition
    B=MR.Transform('xyz','MPS');B=B(1:2,1:2);B(abs(B)<1e-6)=0;B=sign(B);%This may present problems in case of oblique acquisition
    XRange=[-floor(NoX/2) ceil(NoX/2)-1];
    YRange=[-floor(NoY/2) ceil(NoY/2)-1];
    c=1;
    %fprintf('Native: %d\n',Par.Mine.Nat);
    %fprintf('PEs:%s\n',sprintf(' %d',Par.Mine.pedsUn));
    %fprintf('Added samples:%s\n',sprintf(' %d',Par.Mine.AddedSamples));    
    
    %%%%HEREHEREHERE---PERHAPS MORE INTERESTING TO COMPUTE THESE OVER ALL
    %%%%PE'S ALWAYS
    for v=1:4%Par.Mine.pedsUn'
        MR.Parameter.Scan.FoldOverDir=iFoldOverDir{v};%This instruction changes the xyzOffcentres             
        MR.Parameter.Scan.FatShiftDir=iFatShiftDir{v};
        MR.Parameter.Scan.MPS(1:2)=ijkMPS{v}{1};
        MR.Parameter.Scan.MPS(4:5)=ijkMPS{v}{2};
        MR.Parameter.Scan.ijk(1:2)=ijkMPS{v}{1};
        MR.Parameter.Scan.ijk(4:5)=ijkMPS{v}{2};
        for s=1:size(MR.Parameter.Labels.FoldOverDir,1)
            MR.Parameter.Labels.FoldOverDir(s,:)=iFoldOverDir{v};%This instruction changes the xyzOffcentres             
            MR.Parameter.Labels.FatShiftDir(s,:)=iFatShiftDir{v};
        end
        
        %This changes XYZ coordinate system and changes when changing the FoldOver and FatShift so it has to be recomputed
        MR.Parameter.Scan.xyzOffcentres=xyzOffcentres;
        %This changes REC coordinate system but does not change (therefore it does not need to be recomputed)
        %FakedParameters{v}.Scan.Offcentre=Offcentre;

        Par.Scan.MPSOffcentresMM(:,:,c)=MR.Transform([0 0 0],'xyz','MPS')';%This is MPSOffcentresMM        
        %MPSOffcentres=round(2*MRAux.Transform([0 0 0],'xyz','MPSpix')')/2;%This was MPSOffcentres
        Par.Scan.MPSOffcentres(:,:,c)=round(MR.Transform([0 0 0],'xyz','MPSpix')');%This is MPSOffcentres%Without rounding it is not fully integer  
        shifts=A*Par.Scan.MPSOffcentres(1,1:2,c)';
        Par.Encoding.XRange(:,:,c)=XRange-shifts(1);Par.Encoding.YRange(:,:,c)=YRange-shifts(2);        
        
        %FakedParameters{v}.Scan.MPSOffcentres(1,1:2)=centerSorted(v,:);
        MR.Parameter.Scan.MPSOffcentres=Par.Scan.MPSOffcentres(:,:,c);
        MR.Parameter.Scan.MPSOffcentresMM=Par.Scan.MPSOffcentresMM(:,:,c);
        for s=1:size(MR.Parameter.Labels.MPSOffcentres,1)
            MR.Parameter.Labels.MPSOffcentres(s,1:2)=MR.Parameter.Scan.MPSOffcentres(1,1:2);
            MR.Parameter.Labels.MPSOffcentresMM(s,1:2)=MR.Parameter.Scan.MPSOffcentresMM(1,1:2);
        end                                    
        MR.Parameter.Encoding.XRange=Par.Encoding.XRange(:,:,c);MR.Parameter.Encoding.YRange=Par.Encoding.YRange(:,:,c);
        MR.Parameter.Labels.XRange=MR.Parameter.Encoding.XRange;MR.Parameter.Labels.YRange=MR.Parameter.Encoding.YRange;
        if ceil(v/2)==ceil(Par.Mine.Nat/2);MR.Parameter.Scan.FOV=scanFOV;else MR.Parameter.Scan.FOV=flip(scanFOV);end
        C=MR.Transform('xyz','MPS');C=C(1:2,1:2);C(abs(C)<1e-6)=0;C=sign(C);%This may present problems in case of oblique acquisition     
        Par.Mine.Signs(:,:,c)=B/C;%This would convert from MPSOffcentres in current PE to MPSOffcentres in native PE
        Par.Mine=geometryFOV(MR,Par.Mine,c);                       
        %fprintf('Current PE: %d\n',v);
        %fprintf('MPSOffcentres:%s\n',sprintf(' %d',Par.Scan.MPSOffcentres(:,:,c)));
        %fprintf('MPSOffcentresMM:%s\n',sprintf(' %.2f',Par.Scan.MPSOffcentresMM(:,:,c)));
        %fprintf('XRange:%s\n',sprintf(' %d',Par.Encoding.XRange(:,:,c)));
        %fprintf('YRange:%s\n',sprintf(' %d',Par.Encoding.YRange(:,:,c)));
        %fprintf('Geometry rotations:\n%s',sprintf(' %.2f %.2f\n',Par.Mine.Signs(:,:,c)'));
        %fprintf('Isocenter:%s\n',sprintf(' %.2f',Par.Mine.Isoc(:,:,c)));
        %fprintf('Geometry FOV:\n%s',sprintf(' %.2f %.2f %.2f %.2f\n',Par.Mine.APhiAcq(:,:,c)'));
        %fprintf('Geometry FOV:\n%s',sprintf(' %.2f %.2f %.2f %.2f\n',Par.Mine.APhiRec(:,:,c)'));
        c=c+1;
    end
else    
    %We may use this in the future but not yet activated
    %if Par.Mine.UseReadRes
    %    AcqVoxelSizeRec(1:2)=AcqVoxelSize(1)*ones(1:2);
    %    modifyRes;    
    %end    
    Par.Mine.Nat=1;Par.Mine.pedsUn=1;Par.Mine.Signs=[];
    Par.Mine=geometryFOV(MR,Par.Mine,1);
    %fprintf('Isocenter:%s\n',sprintf(' %.2f',Par.Mine.Isoc));
    %fprintf('Geometry FOV:\n%s',sprintf(' %.2f %.2f %.2f %.2f\n',Par.Mine.APhiAcq'));
    %fprintf('Geometry FOV:\n%s',sprintf(' %.2f %.2f %.2f %.2f\n',Par.Mine.APhiRec'));
end

function Par=geometryFOV(MR,Par,c)
    if strcmp(typ,'acq')
        if isempty(MR.Parameter.Encoding.XRes);return;end
        MatrixSize=[MR.Parameter.Encoding.XRes(1) MR.Parameter.Encoding.YRes(1)];       
        if ~isempty(MR.Parameter.Encoding.KzOversampling);MatrixSize(3)=MR.Parameter.Encoding.ZRes(1);else MatrixSize(3)=size(MR.Parameter.Labels.Thicknesses,1)*AccelSlice;end    
        FOV=MatrixSize.*RecVoxelSize;        
        if isempty(MR.Parameter.Encoding.KzOversampling);FOV(3)=(RecVoxelSize(3)+MR.Parameter.Scan.SliceGap(1))*MatrixSize(3)-MR.Parameter.Scan.SliceGap(1);end        
        MatrixSizeRec=round(MatrixSize.*RecVoxelSize./AcqVoxelSizeRec);
        MatrixSize=round(MatrixSize.*RecVoxelSize./AcqVoxelSize);        
        if red;MatrixSize(1:2)=min(MatrixSize(1:2));end
        curFOV=MR.Parameter.Scan.curFOV;
    else%rec
        MatrixSize=[MR.Parameter.Encoding.XReconRes(1) MR.Parameter.Encoding.YReconRes(1)];
        MatrixSizeRec=MatrixSize;
        if ~isempty(MR.Parameter.Encoding.KzOversampling);MatrixSize(3)=MR.Parameter.Encoding.ZReconRes;else MatrixSize(3)=size(MR.Parameter.Labels.Thicknesses,1)*AccelSlice;end   
        FOV=MatrixSize.*RecVoxelSize;
        if isempty(MR.Parameter.Encoding.KzOversampling);FOV(3)=(RecVoxelSize(3)+MR.Parameter.Scan.SliceGap(1))*MatrixSize(3)-MR.Parameter.Scan.SliceGap(1);end
        curijk=MR.Parameter.Scan.ijk;
        curFOV=MR.Parameter.Scan.curFOV;
        MR.Parameter.Scan.ijk=MR.Parameter.Scan.REC;
    end
    MR.Parameter.Scan.curFOV=repmat(FOV,[size(MR.Parameter.Scan.curFOV,1) 1]);
    APhi=MR.Transform('ijk', 'RAF','MatrixSize',MatrixSize);
    APhiRec=MR.Transform('ijk', 'RAF','MatrixSize',MatrixSizeRec);
    AIso=MR.Transform('xyz', 'ijk','MatrixSize',MatrixSize);

    if strcmp(typ,'rec');MR.Parameter.Scan.ijk=curijk;end
    MR.Parameter.Scan.curFOV=curFOV;
    Par.Isoc(:,:,c)=dynInd(AIso,1,3)*[0;0;0;1];
    for n=1:size(APhi,3);APhi=dynInd(APhi,{1:3,4,n},1:3,dynInd(dynInd(APhi,n,3)*[1;1;1;1],1:3,1));end
    for n=1:size(APhiRec,3);APhiRec=dynInd(APhiRec,{1:3,4,n},1:3,dynInd(dynInd(APhiRec,n,3)*[1;1;1;1],1:3,1));end
    Par.APhiAcq(:,:,c)=dynInd(APhi,1,3);
    Par.APhiRec(:,:,c)=dynInd(APhiRec,1,3);
    Par.AIso(:,:,c)=dynInd(AIso,1,3);
end

function modifyRes        
    kRange{1}=Par.Encoding.KxRange;kRange{2}=Par.Encoding.KyRange;%Sampled ranges
    AcqSizeOrig=zeros(1,2);AcqSize=AcqSizeOrig;
    for m=1:2        
        if ~ismember(size(diff(kRange{m},1,2)+1,1),[size(unique(Par.Labels.PartialFourierFactors(:,m)),1) 1]) || ~ismember(size(unique(Par.Labels.PartialFourierFactors(:,m)),1),[size(diff(kRange{m},1,2)+1,1) 1])
            fprintf('Inconsistent encoding for file %s (perhaps due to interaction between direction file and multiple echoes)\n',Par.Filename.Data);fail=1;return;
        end
        AcqSizeOrig(m)=round(max(bsxfun(@times,(diff(kRange{m},1,2)+1),1./unique(Par.Labels.PartialFourierFactors(:,m))),[],1));%This tells the number of points to include in the oversampled/undersampled spectral grid
        AcqSize(m)=round((AcqVoxelSize(m)/AcqVoxelSizeRec(m))*max(bsxfun(@times,(diff(kRange{m},1,2)+1),1./unique(Par.Labels.PartialFourierFactors(:,m))),[],1));%This tells the number of points to include in the oversampled/undersampled spectral grid    
        AcqVoxelSize(m)=AcqVoxelSizeRec(m);
        Par.Scan.AcqVoxelSize(m)=AcqVoxelSize(m);
        Par.Labels.PartialFourierFactors(:,m)=repmat(min((diff(kRange{m},1,2)+1),[],1)/AcqSize(m),[size(Par.Labels.PartialFourierFactors,1) 1]);
        Par.Mine.AddedSamples(m)=AcqSize(m)-AcqSizeOrig(m);
    end
end

end
