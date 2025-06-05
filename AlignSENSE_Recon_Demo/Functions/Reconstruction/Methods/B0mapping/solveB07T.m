
function rec=solveB07T(rec)

%SOLVEB07T   Estimates the B0 field coming from 7T data
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
rec.B0=dynInd(rec.x,2,4).*conj(dynInd(rec.x,1,4));

%MASKING
voxsiz=rec.Enc.AcqVoxelSize;
voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);

parS=rec.Alg.parS;
parS.nErode=0;parS.nDilate=0;
W=refineMask(sqrt(abs(rec.B0)),parS,voxsiz);
distan=30*ones(1,numDims(W));
W=morphFourier(W,distan,rec.Enc.AcqVoxelSize(1:numDims(W)),ones(1,numDims(W)),1);%Soft masking
W(W>1-1e-6)=1;W(W<1e-6)=0;

%UNWRAPPING PARAMETERS
N=size(rec.B0);ND=numDims(rec.B0);
if strcmp(rec.Alg.parU.UnwrapMeth,'PUMA')
    lnorm=2;
    w=[10;10;1;1];w=w(1:ND);
    %fprintf('Used weights:%s\n',sprintf(' %.2f',w));
    c=eye(ND,ND); 
    q=repmat(1-W,[ones(1,ND) ND]);
    rec.B0=puma_hoND(angle(rec.B0),lnorm,'c',c,'q',q,'w',w);
    %DISAMBIGUATE REFERENCE PHASE
    z0=resPop(rec.B0,1:ND,prod(N(1:ND)),1);
    z0=dynInd(z0,W~=0,1);    
    z0=median(z0);
    z0=2*pi*round(z0/(2*pi));%To round to the nearest 2Npi value
    rec.B0=bsxfun(@times,bsxfun(@minus,rec.B0,z0),W);    
elseif strcmp(rec.Alg.parU.UnwrapMeth,'CNCG')
    phaseOrig = angle(rec.B0);
    [rec.B0,D]=CNCGUnwrapping(rec.B0, rec.Enc.AcqVoxelSize, 'Magnitude','LSIt');%MagnitudeGradient4
    plotND([],cat(4, phaseOrig,rescaleND(rec.B0,[-pi,pi]), angle(exp(1i*rec.B0))),[],[],0,{[],2},rec.Par.Mine.APhiRec,{'Original phase';'Unwrapped';'Unwrapped prediced phase'});
    z0=resPop(D,1:ND,prod(N(1:ND)),1);
    z0=dynInd(z0,W~=0,1);    
    z0=median(z0);
    rec.B0=bsxfun(@times,bsxfun(@minus,rec.B0,z0),ones(size(W)));%YB: disables softmasking application since it compromises the B0 maps          
else
    fprintf('Unknown unwrapping method %s\n',rec.Alg.parU.UnwrapMeth);rec.Fail=1;return;
end
% 
% if ~isempty(folderName) && ~isempty(fileName)
%     if ~exist(folderName,'dir');mkdir(folderName);end
%     export_fig(strcat(folderName,filesep,fileName,'.png'));
%     close(h);
% end

%RADIANS TO HERTZ
TE=rec.Par.Labels.TE(2) - rec.Par.Labels.TE(1);%Delta TE
rec.B0=convertB0Field(rec.B0,TE,[],'rad','Hz');

%INCLUDE POSSIBLE THIRD ECHO TO INCREASE SNR
if length(rec.Par.Labels.TE)>2 && rec.Alg.useTE3 && strcmp(rec.Alg.parU.UnwrapMeth,'CNCG')
    TE=rec.Par.Labels.TE(3) - rec.Par.Labels.TE(1);
    P = convertB0Field(rec.B0,TE,[],'Hz','rad');%Phase increment based on current B0 estimate
    PP = dynInd(rec.x,3,4) .* conj( dynInd(rec.x,1,4)) .* conj(exp(1i*P)); %with P phase away to disambiguate
    [PP,D]=CNCGUnwrapping(PP,rec.Enc.AcqVoxelSize,'Magnitude','LSIt');%Unwrapped phase relative to scan 3
    z0=resPop(D,1:ND,prod(N(1:ND)),1);
    z0=dynInd(z0,W~=0,1);z0=median(z0);  
    B0TE3=bsxfun(@times,bsxfun(@minus,PP + P,z0),ones(size(W)));%YB: disables softmasking application since it compromises the B0 maps     
    B0TE3=convertB0Field(B0TE3,TE,[],'rad','Hz');%Need to combine phase first and second estimate
end

%FILTERING
filt=max(voxsiz);%Maximum frequency (in mm)
H=buildFilter(2*N,'tukeyIso',voxsiz/filt,gpu,1,1);
x=filtering(rec.B0,H,1);
rec.B0=cat(4,x,rec.B0);
if rec.Alg.useTE3%Store the B0 map with with using the 3rd echo first in rec.B0
    x = filtering(B0TE3,H,1);
    rec.B0 = cat(4, x, B0TE3 , rec.B0);
end
 rec.B0 = cat(4, rec.B0, abs(W));%Export mask for inspection


