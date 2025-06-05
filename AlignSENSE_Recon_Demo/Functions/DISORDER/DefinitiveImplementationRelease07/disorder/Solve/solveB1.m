function rec=solveB1(rec)

%SOLVEB1   Estimates the B1 field using either DREAM ([1] K Nehrke and P
%Bornert, "DREAM—a novel approach for robust, ultrafast, multislice B1 
%mapping," Magn Reson Med, 68(5):1517-1526, Nov 2012) or AFI ([2] VL 
%Yarnykh, "Actual flip‐angle imaging in the pulsed steady state: A method 
%for rapid three‐dimensional mapping of the transmitted radiofrequency 
%field," Magn Reson Med, 57(1):192-200, Jan 2007)
%   REC=SOLVEB1(REC)
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

voxsiz=rec.Par.Scan.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);

%HARDCODED PARAMETER
dynR=0.1;%Minimum FA ratio

if strcmp(rec.Par.Scan.FastImgMode,'TFE')%We assume it is DREAM
    rec.B=convertRotation(atan(2*abs(dynInd(rec.x,2,4))./(abs(dynInd(rec.x,1,4))+eps)),'rad','deg');%FA from DREAM paper    
    rec.B=rec.B/rec.Par.Labels.FlipAngle(2);%FA normalized
    rec.B(rec.B==0)=1;
else%We assume it is AFI
    rec.B=abs(dynInd(rec.x,2,4))./(abs(dynInd(rec.x,1,4))+eps);
    r=rec.Par.Labels.RepetitionTime(2)/rec.Par.Labels.RepetitionTime(1);
    M=rec.B==0;
    rec.B=convertRotation(acos(min(max((1-rec.B*r)./(rec.B-r),0),1)),'rad','deg'); %min( ,1) because acos  
    rec.B=rec.B/rec.Par.Labels.FlipAngle(1);%FA normalized
    rec.B(M)=1;
end
N=size(rec.B);

%FILTERING
filt=max(voxsiz);%Maximum frequency (in mm)
H=buildFilter(2*N,'tukeyIso',voxsiz/filt,gpu,1,1);
x=filtering(rec.B,H,1);
x=max(x,dynR);
rec.B=cat(4,x,rec.B);
    
rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,11);
