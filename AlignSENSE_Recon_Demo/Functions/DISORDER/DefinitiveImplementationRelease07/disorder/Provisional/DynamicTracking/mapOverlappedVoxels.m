function rec=mapOverlappedVoxels(rec,sh)

% MAPOVERLAPPEDVOXELS shifts among overlapped voxels
%   REC=MAPOVERLAPPEDVOXELS(REC,SH)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%   * {SH} is the amount to shift (from 1 to MB-1)
%   ** REC is a reconstruction structure with information in undistorted 
%   coordinates (rec.u) 
%

AdHocArray=rec.Par.Mine.AdHocArray;
assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
if AdHocArray(1)==101;MB=AdHocArray(4);else MB=1;end

gpu=rec.Dyn.GPU;
if gpu;rec.E=gpuArray(rec.E);end
N=size(rec.E);M=N;
M(2)=rec.Enc.AcqSize(2);
rec.E=resampling(rec.E,M,2);
if AdHocArray(7)~=0;Sh{2}=-sh*rec.Enc.AcqSize(2)/abs(AdHocArray(7));end
Sh{3}=-sh*size(rec.E,3)/MB;
rec.E=shifting(rec.E,Sh);
rec.E=gather(resampling(rec.E,N,2));

    
    
    
