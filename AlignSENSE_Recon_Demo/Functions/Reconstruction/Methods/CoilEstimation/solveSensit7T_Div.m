function rec=solveSensit7T_Div(rec, recB)

%SOLVESENSIT7T_DIV Solves for the sensitivities in 7T data by dividing out
%a (virtual) body coil
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and eigenmaps rec.W
%

if nargin < 2; recB=[];end

%HARDCODED GENERAL PARAMETERS
resolDesired=0;
resol=max(resolDesired,rec.Enc.AcqVoxelSize);%Resolution of sensitivities
NY=size(rec.y);
RSOSAtFullResol = 1; 

if ~isempty(recB)
    %Make body coil
    B = dynInd(recB.y,1,4) - 1i*dynInd(recB.y,2,4);%Quadrature combination
    phU = CNCGUnwrapping(B,recB.Enc.AcqVoxelSize, 'Magnitude','LS');
    M = abs(B)>.2*multDimMax(abs(B));
    ph = phU - multDimMea(phU(M));
    xDiv=abs(B).*exp(1i*ph);
%     xDiv=RSOS(recB.y);
    %Map on FOV of reference scan
    xDiv = mapVolume(xDiv, rec.y, recB.Par.Mine.APhiRec,rec.Par.Mine.APhiRec,[],[],'spline');
else
    xDiv = RSOS(rec.y);
end

if RSOSAtFullResol
%     S = (abs(rec.y).*exp(-1i*angle(rec.y)))./xDiv;
    S = rec.y./xDiv;
    x = SWCC(rec.y, S);
end

%%% RESAMPLE TO REQUIRED RESOLUTION
NResol=round(NY(1:3).*rec.Enc.AcqVoxelSize(1:3)./resol);
%rec.y=filtering(rec.y,buildFilter(NY(1:ND),'tukeyIso',1,gpu,0.05));

if RSOSAtFullResol
    rec.S = resampling(S,NResol);
    rec.x = resampling(x,NResol);
else
    rec.y=resampling(rec.y,NResol);
    xDiv=resampling(xDiv,NResol);
end

%%% CHANGE HEADER
[rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec] = mapNIIGeom( rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec, 'resampling',[], NY(1:3), NResol(1:3));
rec.Enc.FOVSize = NResol(1:3);
%rec.y=filtering(rec.y,buildFilter(NResol,'tukeyIso',1,gpu,0.05),1);

%%% COIL ESTIMATION
if ~RSOSAtFullResol
    rec.S = rec.y./xDiv;
    rec.x = SWCC(rec.y, rec.S);
end

%%% REPORT
rec = gatherStruct(rec);

end