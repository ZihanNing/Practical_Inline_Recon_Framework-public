function [rec,perm,fl]=mps2Rec(rec,di,minRes,typ)

%MPS2REC   Transforms data from MPS to REC coordinates
%   [REC,PERM,FL]=MPS2REC(REC,{DI},{MINRES},{TYP})
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), the
%   .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)), the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * {DI} is a the direction of transform (1, default, for mps2Rec and 0 for rec2Mps)
%   * {MINRES} is a flag to denote whether the data has to be resampled to the minimum of the acquired resolution, it defaults to 0
%   * {TYP} is the type of operation, with 1 a full transformation of the data is performed. With 0 only the permutation and flips are returned and
%   the data is left untouched
%   * REC is a reconstruction structure with transformed (or not, depending on TYP) coordinates
%   * PERM is the set of permutations to transform the data from mps2Rec or rec2Mps
%   * FL is the set of flips to transform the data from mps2Rec or rec2Mps
%

if ~exist('di','var') || isempty(di);di=1;end
if ~exist('minRes','var') || isempty(minRes);minRes=0;end
if ~exist('typ','var') || isempty(typ);typ=1;end

typ2Rec=rec.Dyn.Typ2Rec;name=rec.Names.Name;

%PRECOMPUTATIONS FOR COORDINATE TRANSFORMS
[Orientation,FoldOverDir,FatShiftDir]=parUnaFun({rec.Par.Scan.Orientation,rec.Par.Scan.FoldOverDir,rec.Par.Scan.FatShiftDir},@unique,'rows');
assert(size(Orientation,1)==1,'Different orientations %s in different stacks for file %s',permute(Orientation,[2 1]),name);
assert(size(FoldOverDir,1)==1,'Different foldover dirs %s in different stacks for file %s',permute(FoldOverDir,[2 1]),name);
assert(size(Orientation,1)==1,'Different fat shift dirs %s in different stacks for file %s',permute(FatShiftDir,[2 1]),name);
if strcmp(rec.Par.Scan.AcqMode,'Cartesian') || strcmp(rec.Par.Scan.AcqMode,'EPI')
    if strcmp(Orientation,'SAG')
        if strcmp(FoldOverDir,'AP')
            perm=[1 2 3];
            if strcmp(FatShiftDir,'F') || strcmp(FatShiftDir,'A');fl=[0 1 1];else fl=[1 0 1];end%H-P
        else%FH
            perm=[2 1 3];
            if strcmp(FatShiftDir,'A') || strcmp(FatShiftDir,'H');fl=[1 1 1];else fl=[0 0 1];end%P-F
        end
    elseif strcmp(Orientation,'TRA')
        if strcmp(FoldOverDir,'AP')
            perm=[2 1 3];
            if strcmp(FatShiftDir,'L') || strcmp(FatShiftDir,'P');fl=[0 0 1];else fl=[1 1 1];end%R-A
        else%RL
            perm=[1 2 3];
            if strcmp(FatShiftDir,'A') || strcmp(FatShiftDir,'L');fl=[1 0 1];else fl=[0 1 1];end%P-R
        end
    else%COR
        if strcmp(FoldOverDir,'RL')
            perm=[1 2 3];
            if strcmp(FatShiftDir,'F') || strcmp(FatShiftDir,'R');fl=[0 1 1];else fl=[1 0 1];end%H-L
        else%FH
            perm=[2 1 3];
            if strcmp(FatShiftDir,'L') || strcmp(FatShiftDir,'F');fl=[0 0 1];else fl=[1 1 1];end%R-H
        end
    end
elseif strcmp(rec.Par.Scan.AcqMode,'Radial');perm=[1 2 3];fl=[1 0 1];
elseif strcmp(rec.Par.Scan.AcqMode,'Kooshball');perm=[1 2 3];fl=[0 0 1];
elseif strcmp(rec.Par.Scan.AcqMode,'Spiral');perm=[2 1 3];fl=[0 0 1];
end
%This perhaps to be introduced because for GT convention the z axis is not flipped...
fl(3)=0;
%Actually, in the transforms given by transformRAF, for data at acquired resolution mps2Rec does not seem necessary, review if this function needs to 
%be called

if ~typ;return;end

%CHANGE OF COORDINATES
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if n>=6 && n<=13
        permin=perm;permin(4:ndims(rec.(datTyp)))=4:ndims(rec.(datTyp));
        if di
            rec.(datTyp)=permute(rec.(datTyp),permin);
            for m=1:3
                if fl(m);rec.(datTyp)=flip(rec.(datTyp),m);end
            end
        else
            for m=1:3
                if fl(m);rec.(datTyp)=flip(rec.(datTyp),m);end
            end
            rec.(datTyp)=ipermute(rec.(datTyp),permin);
        end

        %TODO, THIS MAY BE DONE OUTSIDE THIS FUNCTION
        if minRes                        
            M=size(rec.(datTyp));
            N=min(M(1:2));
            if n==8;rec.(datTyp)=single(abs(rec.(datTyp))>1e-12);end                
            rec.(datTyp)=resampling(rec.(datTyp),[N N]);
            if n==8;rec.(datTyp)=single(abs(rec.(datTyp))>0.5);end
        end
    end
end

%MASKING IN THE NEW COORDINATE SYSTEM
for n=[9 11:13];datTyp=rec.Plan.Types{n};
    if minRes && isfield(rec,datTyp) && any(typ2Rec==8);rec.(datTyp)=bsxfun(@times,rec.(datTyp),rec.W);end
end
