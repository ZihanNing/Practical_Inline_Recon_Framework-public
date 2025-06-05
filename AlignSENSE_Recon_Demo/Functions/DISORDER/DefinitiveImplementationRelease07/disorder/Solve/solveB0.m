function rec=solveB0(rec)

%SOLVEB0   Estimates the B0 field
%   REC=SOLVEB0(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   * REC is a reconstruction structure with estimated sensitivities rec.S
%   and mask rec.M
%

gpu=isa(rec.x,'gpuArray');

%COMPILATION OF THE GC CODE---THIS SHOULD BE PRECOMPILED IN THE FUTURE
%GCFold=fullfile(fileparts(mfilename('fullpath')),'../PhaseUnwrapping','to_compile_mf2');
%mex('-silent','-outdir',GCFold,fullfile(GCFold,'mf2.cpp'),fullfile(GCFold,'graph.cpp'),fullfile(GCFold,'maxflow.cpp'));
rec.B=dynInd(rec.x,2,4).*conj(dynInd(rec.x,1,4));

%MASKING
voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);

parS=rec.Alg.parS;
parS.nErode=0;parS.nDilate=0;
W=refineMask(sqrt(abs(rec.B)),parS,voxsiz);
distan=30*ones(1,numDims(W));
W=morphFourier(W,distan,rec.Enc.AcqVoxelSize(1:numDims(W)),ones(1,numDims(W)),1);%Soft masking
W(W>1-1e-6)=1;W(W<1e-6)=0;

%UNWRAPPING PARAMETERS
N=size(rec.B);ND=numDims(rec.B);
if strcmp(rec.Alg.parU.UnwrapMeth,'PUMA')
    lnorm=2;
    w=[10;10;1;1];w=w(1:ND);
    %fprintf('Used weights:%s\n',sprintf(' %.2f',w));
    c=eye(ND,ND); 
    q=repmat(1-W,[ones(1,ND) ND]);
    rec.B=puma_hoND(angle(rec.B),lnorm,'c',c,'q',q,'w',w);
    %DISAMBIGUATE REFERENCE PHASE
    z0=resPop(rec.B,1:ND,prod(N(1:ND)),1);
    z0=dynInd(z0,W~=0,1);    
    z0=median(z0);
    z0=2*pi*round(z0/(2*pi));%To round to the nearest 2Npi value
    rec.B=bsxfun(@times,bsxfun(@minus,rec.B,z0),W);    
elseif strcmp(rec.Alg.parU.UnwrapMeth,'CNCG')
    [rec.B,D]=CNCGUnwrapping(rec.B,1,rec.Enc.AcqVoxelSize);
    z0=resPop(D,1:ND,prod(N(1:ND)),1);
    z0=dynInd(z0,W~=0,1);    
    z0=median(z0);
    rec.B=bsxfun(@times,bsxfun(@minus,rec.B,z0),W);     
else
    fprintf('Unknown unwrapping method %s\n',rec.Alg.parU.UnwrapMeth);rec.Fail=1;return;
end

%RADIANS TO HERTZ
TE=2.3;
if rec.Dyn.Debug>=2;fprintf('TE has been hardcoded to %.02f due to inconsistencies in ReconFrame information\n',TE);end
rec.B=convertB0Field(rec.B,TE,[],'rad','Hz');

%FILTERING
filt=max(voxsiz);%Maximum frequency (in mm)
H=buildFilter(2*N,'tukeyIso',voxsiz/filt,gpu,1,1);
x=filtering(rec.B,H,1);
rec.B=cat(4,x,rec.B);

rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,11);
