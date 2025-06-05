
function rec = prepareRec(fileAcq, fileRef, fileB1, fileB0, writeNIIFlag, pathOu, resRec, supportReadout)

%PREPAREREC prepares a reconstruction structure that can be called by a reconstruction algorithm. It combines acquistition data with the corresponding reference data (Sensitivities/B0/B1). 
%   REC = PREPAREREC(FILEACQ,FILEREF,{FILEB1},{FILEB0},{WRITENIIFLAG},{PATHOU},{RESREC})
%   * FILEACQ is the name of the raw data (.dat) of the acquisition.
%   * FILEREF is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * {FILEB1} is the name of the raw data (.dat) of the reference acquisition used for B1 mapping
%   * {FILEB0} is the name of the raw data (.dat) of the reference acquisition used for B0 mapping
%   * {WRITENIIFLAG} is a flag to write a simple Sensitivity Weighted Coil Combined recon to a NIFTI files. Defaults to 1.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the files are stored.
%   * {RESREC} is the resolution at which to downsample the acquisition.
%   ** REC is the reconstruction structure with all the data arranged.
%   

if nargin<2 || isempty(fileRef);error('prepareRec:: Sensitivities must be provided.'); end
if ~(nargin<3) && ~isempty(fileB0);[~,fileb0]=fileparts(fileB0);else; fileb0=''; end
if ~(nargin<4) && ~isempty(fileB1);[~,fileb1]=fileparts(fileB1);else; fileb1=''; end


fprintf('<strong>Creating reconstruction structure:</strong>\n');
if ~iscell(fileAcq)%Single acquisition
    [pathOuTemp,file]=fileparts(fileAcq);
    fprintf('  Acquisition:   %s\n',file);
else%Multiple acquisitions
    for i=1:length(fileAcq)
        [pathOuTemp,file]=fileparts(fileAcq{i});
        fprintf('  Acquisition:   %s\n',file);
    end
end

assert(~iscell(fileRef),'prepareRec:: The reference file fileRef cannot contain multiple scans. It must be a singe filename.')
[~,fileref]=fileparts(fileRef);
fprintf('  Sensitivities:   %s\n',fileref);
fprintf('  B0 map:   %s\n',fileb0);
fprintf('  B1 map:   %s\n',fileb1);

if nargin<2 || isempty(fileRef); error('matchRec:: Reference data needs to be provided.');end
if nargin<3 || isempty(fileB1); includeB1 =0;else includeB1 =1;end
if nargin<4 || isempty(fileB0); includeB0 =0;else includeB0 =1;end
if nargin<5 || isempty(writeNIIFlag); writeNIIFlag=1;end
if nargin<6 || isempty(pathOu); pathOu=[];end
if nargin<7 || isempty(resRec); resRec=[];end
if nargin<8 || isempty(supportReadout); supportReadout=[];end

if isempty(pathOu); pathOu = pathOuTemp;end

%%% ASSING DATA  
if ~iscell(fileAcq)%Single acquisition
    if exist(strcat(fileAcq,'.mat'),'file')~=2; dat2Rec(fileAcq);end
    ss= load(strcat(fileAcq,'.mat'));rec=ss.rec;ss=[]; 
    %Resample to specific resolution
    if ~isempty(resRec) || ~isempty(supportReadout); rec=resampleRec(rec,resRec,supportReadout);end
else%Combine multiple scans
    for n=1:length(fileAcq)%Run over acquistitions
        if exist(strcat(fileAcq{n},'.mat'),'file')~=2; dat2Rec(fileAcq{n});end
        ss=load(strcat(fileAcq{n},'.mat'));rec=ss.rec;ss=[];  
        %Resample to specific resolution
        if ~isempty(resRec) || ~isempty(supportReadout); rec=resampleRec(rec,resRec,supportReadout);end
        if n==1
           y=rec.y;
           z{2}=rec.Assign.z{2};
           z{3}=rec.Assign.z{3};
           Name=rec.Names.Name;
           if isfield(rec,'PT'); pSliceImage=rec.PT.pSliceImage;end
           if isfield(rec,'PT'); pTimeTest=rec.PT.pTimeTest;end
        else
           y=cat(5,y,rec.y);
           z{2}=cat(2,z{2},rec.Assign.z{2});
           z{3}=cat(3,z{3},rec.Assign.z{3});
           Name=strcat(Name,rec.Names.Name);
           if isfield(rec,'PT'); pSliceImage=cat(5,pSliceImage,rec.PT.pSliceImage);end
           if isfield(rec,'PT'); pTimeTest=cat(2,pTimeTest,rec.PT.pTimeTest);end
        end
    end
    rec.y=y;y=[];
    rec.Assign.z{2}=z{2};z{2}=[];
    rec.Assign.z{3}=z{3};z{3}=[];
    rec.Names.Name=Name;
    if isfield(rec,'PT'); rec.PT.pSliceImage=pSliceImage;end
    if isfield(rec,'PT'); rec.PT.pTimeTest=pTimeTest;end
end

%%% ASSIGN COIL SENSITIVITIES
if exist(strcat(fileRef,'.mat'),'file')~=2; dat2Se(fileRef);end%Check if data has been parsed
load(strcat(fileRef,'.mat'), 'recS');if ~exist('recS','var'); recS = dat2Se(fileRef);end %Check if sensitivities have been estimated
NS=size(recS.S);NS=NS(1:3);
NY=size(rec.y);NY=NY(1:3);

%%% FILTER COIL SENITIVITIES
gibbsRinging=1;
if gibbsRinging~=0 && gibbsRinging<=1;HY=buildFilter(2*NS,'tukeyIso',0.5,0,gibbsRinging,1);elseif gibbsRinging==-1;HY=buildFilter(2*NS,'CubicBSpline',[],0,1);else HY=[];end
if ~isempty(HY)
    rec.S=filtering(recS.S,HY,1);
    HY=buildFilter(NS,'tukeyIso',1,0,gibbsRinging);
    rec.M=abs(filtering(abs(recS.x),HY));%Mask from reference is more likely to be motion-free (shorter scan) & fully sampled          
else            
    rec.S=recS.S;rec.W=recS.W;rec.M=recS.x;
end

%%% MAP FOV OF ACQUISITION DATA AND REFERENCE DATA
MTS = recS.Par.Mine.APhiRec;MTy = rec.Par.Mine.APhiRec;
rec.S = mapVolume (rec.S, ones(NY(1:3)), MTS, MTy,[],[],'spline');
rec.M = gather(rec.M);
rec.M = mapVolume (rec.M, ones(NY(1:3)), MTS, MTy,[],[],'spline');

%%% CHECK IF GEOMTRIES ARE CORRECT
xTest = mapVolume (recS.x, ones(NY(1:3)), MTS, MTy,[],[]);
plotND([], cat(4,xTest,SWCC(multDimMea(rec.y,5:16),rec.S)),[0 15],[],0,[],MTy,{sprintf('Reference scan: %s',recS.Names.Name);sprintf('RSOS of acquisition: %s',rec.Names.Name)},single(abs(xTest)>2),{2},1000);
title('Test to see if geometries are consistent')
if ~exist(fullfile( rec.Names.pathOu ,'An-Ve_Sn','GeomTesting'),'dir'); mkdir(fullfile( rec.Names.pathOu ,'An-Ve_Sn','GeomTesting'));end
%saveFig(fullfile( rec.Names.pathOu ,'An-Ve_Sn','GeomTesting',rec.Names.Name)) ;
recS=[];xTest=[];

%%% STANDARDISE COILS - I don't think it is necessary nice you have de-correlated both k-space of reference and acquisition data
%if isfield(rec,'N');rec.S=standardizeCoils(rec.S,rec.N);end
    
%     %%% ASSIGN B0 map
%     if ~exist(strcat(fileB0,'.mat'),'file'); dat2B0(fileB0);end
%     ss = load(strcat(fileB0,'.mat'));recB0 =ss.recB0; ss=[];
%     MTB0 = recB0.Par.Mine.APhiRec;
%     rec.B0 = mapVolume (dynInd(recB0.B0,1,4), ones(N(1:3)) , MTB0, MTy);recB0=[]; %Only extract one B0 map   
%     
%     %%% ASSIGN B1 map
%     if ~exist(strcat(fileB1,'.mat'),'file'); dat2B1(fileB1);end
%     ss = load(strcat(fileB1,'.mat'));recB1 =ss.recB1; ss=[];
%     MTB1 = recB1.Par.Mine.APhiRec;
%     rec.B = mapVolume (dynInd(recB1.B,1,4), ones(N(1:3)) , MTB1, MTy);recB1=[];%Only extract one B1 map   
 

%%% RECONSTRUCT IMAGE AND WRITE TO INSPECT SIMPLEST PARALLEL IMAGING RECON
if writeNIIFlag
    recW = rec;
    x = SWCC( multDimMea(recW.y,5), recW.S);%Sensitivity Weighted Coil Combination where repeats are averaged
    recW.Par.Mine.Modal=7;%Writing data to An-Ve = Volumetric encoding (3D)
    recW.Plan.Suff = '_SWCC'; recW.Plan.SuffOu='';recW.Names.pathOu = pathOu;
    recW.Dyn.Typ2Wri=zeros(1,50);
    recW.x=x;recW.Plan.Types{10}='x';recW.Plan.TypeNames{10}='Aq';recW.Dyn.Typ2Wri(10)=1;recW.Dyn.Typ2Rec=vertcat([],10);
    recW.Alg.OnlyJSON = 0;recW.Fail=0;
    recW.Par.Mine.Proce=1;recW.Par.Mine.Nat=1;recW.Alg.OverDec=ones(1,3);recW.Par.Mine.Signs=[];
    writeData(recW);
    fprintf('Sensitivity Weighted Coil Combination reconstruction saved to NIFTI file.\n')
end

