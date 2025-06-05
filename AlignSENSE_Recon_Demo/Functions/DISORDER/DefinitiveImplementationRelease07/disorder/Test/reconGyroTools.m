function recGT=reconGyroTools(rec,step)

%RECONGYROTOOLS   Performs a GyroTools reconstruction to compare the methods against
%   REC=RECONGYROTOOLS(REC,{STEP})
%   * REC is a reconstruction structure. At this stage it should contain the naming information (.Names) and the plan information (.Plan)
%   * STEP is the step at which to output the reconstructions, empty, default, to write the recons to NIFTI
%   * RECGT is a GyroTools reconstruction object
%

if ~exist('step','var');step=[];end

%GENERIC RECONFRAME RECONSTRUCTION

%FIRST WE DETECT THE REFERENCE DATA
indSENSE=find(rec.Names.prot.A_Modal(1:rec.Names.ind)==2);
targetFile=strtrim(rec.Names.prot.B_FileName(indSENSE(end),:));
senseRawFile=fullfile(rec.Names.pathIn,sprintf('%s.raw',targetFile));
rawFile=rec.Names.rawFile;
recGTS= MRsense(senseRawFile,rawFile);
recGTS.Perform;
if strcmp(step,'Sensitivities');recGT=recGTS;return;end

%NOW WE PERFORM THE RECONSTRUCTION
recGT=MRecon(rawFile);if strcmp(step,'Parse');return;end
recAux.Par.Mine.diInfo=[];
recAux.Par=computeRAF(MR,recAux.Par,'rec');
recGT.Parameter.Recon.Sensitivities=recGTS;
if ~strcmp(recGT.Parameter.Labels.FastImagingMode,'EPI') && ~strcmp(recGT.Parameter.Labels.FastImagingMode,'GRASE');recGT.Parameter.Parameter2Read.typ = 1;else recGT.Parameter.Parameter2Read.typ=double(rec.Dyn.Typ2Rec);end%Otherwise for anatomicals Oversampling Removal in image space throws an error...
for n=2:rec.Plan.NDims;recGT.Parameter.Parameter2Read.(rec.Plan.Dims{n})=double(rec.Dyn.Ind2Rec{n});end
recGT.ReadData;if strcmp(step,'Read');return;end
recGT.RandomPhaseCorrection;if strcmp(step,'RandomPhase');return;end
recGT.RemoveOversampling;if strcmp(step,'RemoveMOversampling');end%This is not in our code, remove for comparison
recGT.PDACorrection;if strcmp(step,'PDA');return;end
recGT.DcOffsetCorrection;if strcmp(step,'DCOffset');return;end
recGT.MeasPhaseCorrection;if strcmp(step,'MeasPhase');return;end
recGT.SortData;if strcmp(step,'Sort');return;end
recGT.GridData;if strcmp(step,'Grid');return;end
recGT.RingingFilter;if strcmp(step,'RingingFilter');return;end%This is not in our code, remove for comparison
recGT.ZeroFill;if strcmp(step,'ZeroFillK');return;end%This is not in our code, remove for comparison
recGT.K2IM;if strcmp(step,'K2IM');return;end
recGT.EPIPhaseCorrection;if strcmp(step,'EPIPhase');return;end
%recGT.Parameter.Recon.RemovePOversampling='No';%In our pipeline we would do something like this
%recGT.RemoveOversampling;if strcmp(step,'RemoveMOversampling');return;end%In our pipeline we would do something like this
recGT.K2IP;if strcmp(step,'K2IP');return;end
recGT.GridderNormalization;if strcmp(step,'GridNormalization');return;end
%recGT.RemoveOversampling;if strcmp(step,'RemovePOversampling');return;end%In our pipeline we would do something like this (not exactly but
%reconstructing to the final FOV directly)
recGT.SENSEUnfold;if strcmp(step,'SENSEUnfold');return;end
recGT.PartialFourier;if strcmp(step,'PartialFourier');return;end
recGT.Average;if strcmp(step,'Average');return;end
recGT.GeometryCorrection;if strcmp(step,'GeometryCorrection');return;end
recGT.RemoveOversampling;if strcmp(step,'RemovePOversampling');return;end%This is not in our code, remove for comparison
recGT.ZeroFill;if strcmp(step,'ZeroFillI');return;end%This is not in our code, remove for comparison
recGT.FlowPhaseCorrection;if strcmp(step,'FlowPhaseCorrection');return;end%This is not in our code, remove for comparison
recGT.RotateImage;if strcmp(step,'Rotate');return;end

%ASSIGN INFORMATION OF INTEREST TO AN AUXILIARY REC STRUCTURE AND WRITE THE FILE
recAux.Fail=rec.Fail;
recAux.Par.Mine.Modal=12;
recAux.Names=rec.Names;
recAux.Dyn.Typ2Rec=12;
recAux.Dyn.Typ2Wri=rec.Dyn.Typ2Wri;
recAux.Plan.Types=rec.Plan.Types;
recAux.Plan.TypeNames=rec.Plan.TypeNames;
recAux.Plan.Suff=rec.Plan.Suff;

recAux.x=recGT.Data;
if size(recAux.x,8)~=1;recAux.x=permute(recAux.x,[1 2 8 3 4 5 6 7 9 10 11 12]);end
writeData(recAux);
