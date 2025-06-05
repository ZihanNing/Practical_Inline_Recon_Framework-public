function rec=maskFromShim(rec)

% MASKFROMSHIM builds a brain mask from shimming box isocenter
%   REC=MASKFROMSHIM(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * REC is a reconstruction structure with mask information (rec.M) 
%

name=rec.Names.Name;
%MASKING EITHER BASED ON SHIM CENTER, ON G-FACTOR MAPS OR ON THE DATA ITSELF
if rec.Alg.parU.useShimBox
    if rec.Dyn.Debug>=1 && ~isfield(rec.Par.Mine,'ShimBoxCenter');fprintf('Shim box center not found for file %s\n',name);rec.Fail=1;return;end
    sbCenter=rec.Par.Mine.ShimBoxCenter;sbCenter(4)=1;
    %sbCenter=dynInd(rec.Par.Mine.APhiRec,rec.Par.Mine.Nat==rec.Par.Mine.pedsUn,3)\sbCenter(:);
    sbCenter=dynInd(rec.Par.Mine.APhiRec,rec.Par.Mine.Nat,3)\sbCenter(:);
    sbCenter=round(sbCenter(1:3)+1);
end
rec.M=dynInd(rec.x,1,4);
N=size(rec.M);NDims=numDims(rec.M);NDims=min(NDims,3);
voxsiz=rec.Par.Scan.AcqVoxelSize(1,1:NDims);
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
mirr=ones(1,NDims);
NPE=length(rec.Par.Mine.pedsUn);
%NPE=size(rec.Par.Mine.APhiRec,3);
if isfield(rec,'G');rec.G=reshape(rec.G,[N(1:3) NPE 2]);end
if rec.Alg.parU.useShimBox    
    rec.M=single(rec.M>1e-3);
    rec.M=repmat(rec.M,[1 1 1 1 2]);
    rec.M=dynInd(rec.M,2,5,0);
    assert(all(sbCenter>0.5) && all(sbCenter'<N+0.5),'Shim box center (%s) outside the FOV (%s) for file %s',sprintf(' %d',sbCenter),sprintf(' %d',N),name);
    rec.M=dynInd(rec.M,[sbCenter' 2],[1:3 5],1);           
    dilate=60;dilate=dilate*[1 1.2 1];dilate(NDims+1:end)=[];
    rec.M=dynInd(rec.M,2,5,morphFourier(dynInd(rec.M,2,5),dilate,voxsiz,mirr));
    if isfield(rec,'G');rec.M=dynInd(rec.M,1,5,dynInd(rec.G,2,5));end
    rec.M=dynInd(rec.M,1,5,single(dynInd(rec.M,1,5)>1e-3));
elseif isfield(rec,'G')   
    rec.M=dynInd(rec.G,2,5);
    rec.M=single(rec.M>1e-3);
    rec.M=cat(5,rec.M,refineMask(rec.M,rec.Alg.parS,voxsiz));
else
    rec.M=single(rec.M>1e-3);
    rec.M=cat(5,rec.M,refineMask(rec.M,rec.Alg.parS,voxsiz));
end
rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,8);rec.Dyn.Typ2Wri(8)=1;

%ROI COMPUTATION
rec.Enc.ROI=computeROI(dynInd(rec.M,2,5));
rec.Enc.ROI(3,:)=[1 N(3) N(3) N(3) 0 0];%To avoid problems later on with MB replication
if rec.Dyn.Debug>=1;fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));end

if NPE>2;rec.Enc.ROI(1,:)=[1 N(1) N(1) N(1) 0 0];end%To avoid problems with partial Fourier

%ROI EXTRACTION
typ2Rec=rec.Dyn.Typ2Rec;
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,[1 3]);
end
rec.Enc.FOVSize([1 3])=rec.Enc.ROI([1 3],4)';