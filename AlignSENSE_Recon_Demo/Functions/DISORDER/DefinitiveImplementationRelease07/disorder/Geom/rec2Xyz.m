function [rec,perm,fl]=rec2Xyz(rec,di,typ)

%REC2XYZ   Transforms data from REC to XYZ coordinates
%   [REC,PERM,FL]=REC2XYZ(REC,{DI},{TYP})
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), the
%   .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)), the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * {DI} is a the direction of transform (1, default, for rec2Xyz and 0 for xyz2Rec)
%   * {TYP} is the type of operation, with 1 a full transformation of the data is performed. With 0 only the permutation and flips are returned and
%   the data is left untouched
%   * REC is a reconstruction structure with transformed coordinates
%   * PERM is the set of permutations to transform the data from rec2Xyz or xyz2Rec
%   * FL is the set of flips to transform the data from rec2Xyz or xyz2Rec
%

if ~exist('di','var') || isempty(di);di=1;end
if ~exist('typ','var') || isempty(typ);typ=1;end

typ2Rec=rec.Dyn.Typ2Rec;name=rec.Names.Name;

%PRECOMPUTATIONS FOR COORDINATE TRANSFORMS
[Orientation,PatientPosition,PatientOrientation]=parUnaFun({rec.Par.Scan.Orientation,rec.Par.Scan.PatientPosition,rec.Par.Scan.PatientOrientation},@unique,'rows');
assert(size(Orientation,1)==1,'Different orientations %s in different stacks for file %s',permute(Orientation,[2 1]),name);
assert(size(PatientPosition,1)==1,'Different patient positions %s in different stacks for file %s',permute(PatientPosition,[2 1]),name);
assert(size(PatientOrientation,1)==1,'Different patient orientations %s in different stacks for file %s',permute(PatientOrientation,[2 1]),name);

if strcmp(Orientation,'SAG')
    if strcmp(PatientOrientation,'Supine')
        perm=[2 3 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 1 1];else fl=[1 0 0];end%FeetFirst
    elseif strcmp(PatientOrientation,'Prone')
        perm=[2 3 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 0 1];else fl=[0 1 0];end%FeetFirst
    elseif strcmp(PatientOrientation,'Left')
        perm=[3 2 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 1 1];else fl=[0 0 0];end%FeetFirst
    else%Right
        perm=[3 2 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 0 1];else fl=[1 1 0];end%FeetFirst
    end
elseif strcmp(Orientation,'TRA')
    if strcmp(PatientOrientation,'Supine')
        perm=[1 2 3];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 0 0];else fl=[1 1 1];end%FeetFirst
    elseif strcmp(PatientOrientation,'Prone')
        perm=[1 2 3];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 1 0];else fl=[0 0 1];end%FeetFirst
    elseif strcmp(PatientOrientation,'Left')
        perm=[2 1 3];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 1 0];else fl=[1 0 1];end%FeetFirst
    else%Right
        perm=[2 1 3];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 0 0];else fl=[0 1 1];end%FeetFirst
    end
else%COR
    if strcmp(PatientOrientation,'Supine')
        perm=[3 2 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 0 1];else fl=[1 1 0];end%FeetFirst
    elseif strcmp(PatientOrientation,'Prone')
        perm=[3 2 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 1 1];else fl=[0 0 0];end%FeetFirst
    elseif strcmp(PatientOrientation,'Left')
        perm=[2 3 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[1 1 1];else fl=[1 0 0];end%FeetFirst
    else%Right
        perm=[2 3 1];
        if strcmp(PatientPosition,'HeadFirst');fl=[0 0 1];else fl=[0 1 0];end%FeetFirst
    end
end

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
    end
end
