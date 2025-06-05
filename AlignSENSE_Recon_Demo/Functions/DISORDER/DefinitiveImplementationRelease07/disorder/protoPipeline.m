function fRec=protoPipeline(path,modal,rootOu,rootIn,series,specific,repeated)

%PROTOPIPELINE   Runs the reconstruction pipeline for a given study
%   FREC=PROTOPIPELINE(PATH,{MODAL},{ROOTOU},{ROOTIN},{SERIES},{SPECIFIC},{REPEATED})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%   * {ROOTIN} is the folder where the raw data has been mounted. It
%   defaults to Data/rawSource (relative to the user's home folder)
%   * {SERIES} restricts the reconstructions to a specific set of series
%   * {SPECIFIC} indicates to use a specific configuration of parameters as
%   stated in reconSpecific.m
%   * {REPEATED} looks for joint reconstructions of repeated scans
%   ** FREC returns the rec structures where the method failed
%

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('\nReconstructing study %s\n',path);
%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode,~,pathRF]=versRecCode;
%ADD THE RECONFRAME CODE TO THE PATH (REQUIRED ONLY FOR TESTING)
if exist(pathRF,'dir') && ~isempty(pathRF);addpath(genpath(pathRF));end

if nargin<2;modal=[];end
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);
if nargin<4 || isempty(rootIn);rootIn=fullfile(filesep,'home',user,'Data','rawSource');end
%DETECTING THE RAW FOLDER
pathIn=rawFolderDetection(path,rootIn);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end
%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');parseLabFolder(path,[],rootOu,rootIn);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

if nargout>0;fRec=[];end
contF=1;
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin>=5 && ~isempty(series);nV=nV(ismember(nV,series));end
if nargin<6;specific=[];end
if nargin<7 || isempty(repeated);repeated=0;end
interveneNext=0;
rep=1;
for n=nV
    if interveneNext~=2        
        rec.Rec=1;
        rec.Names.Name=strtrim(prot.B_FileName(n,:));
        matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));rawFile=fullfile(pathIn,sprintf('%s.raw',rec.Names.Name));
        %if ~isempty(regexp(matFile,'t2mbz','once')) 
            rec.Names.rawFile=rawFile;rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.pathIn=pathIn;rec.Names.path=path;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;rec.Names.Specific=specific;
            if exist(rawFile,'file')                
                load(matFile);                        
                rec.Par=Par;Par=[];
                pathModa=fullfile(rec.Names.pathOu,numbe2Modal(rec.Par.Mine.Modal));
                if ~exist(pathModa,'dir');mkdir(pathModa);end%We create the modality folder at the beginning to detect that an execution is in process                
                if isempty(gcp('nocreate'));evalc('parpool');end
                if interveneNext==1
                    rec.Par.Mine.InitVol=initVol;
                    rec.Par.Mine.ShiftIsoc=shiftIsoc;
                end
                interveneNext=0;
                if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') && rec.Par.Mine.Modal==10 && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1%MB DWI neonates                              
                    NVS=size(rec.Par.Mine.diInfo,1);
                    if rec.Par.Parameter2Read.extr2(end)~=NVS-1%Incomplete Split Scan, we read the following ones
                        Alg=reconAlgorithm;                    
                        curV=rec.Par.Parameter2Read.extr2(end);
                        if Alg.JointSplitScans
                            s=n;                        
                            rec.Par.Mine.splitInfo=zeros(NVS,1);                                        
                            rec.recV{s-n+1}=rec;
                            rec.recV{s-n+1}.Par.Mine.InitVol=0;
                            rec.Par.Mine.splitInfo(rec.Par.Parameter2Read.extr2+1)=s-n+1;
                            fprintf('Series %s: %d vols\n',rec.Names.Name,rec.Par.Parameter2Read.extr2(end)+1);
                            s=s+1;
                            if s<=nV(end)
                                recAux.Names.Name=strtrim(prot.B_FileName(s,:));
                                matFile=fullfile(headFolder,sprintf('%s.mat',recAux.Names.Name));rawFile=fullfile(pathIn,sprintf('%s.raw',recAux.Names.Name));
                                recAux.Names.rawFile=rawFile;recAux.Names.matFile=matFile;recAux.Names.pathOu=pathOu;recAux.Names.pathIn=pathIn;recAux.Names.path=path;recAux.Names.prot=prot;recAux.Names.ind=n;recAux.Names.headFolder=headFolder;recAux.Names.versCode=versCode;recAux.Names.Specific=specific;
                                if exist(rawFile,'file')
                                    load(matFile);
                                    recAux.Par=Par;Par=[];
                                    if recAux.Par.Mine.Modal==10 && recAux.Par.Mine.AdHocArray(1)==101 && recAux.Par.Mine.AdHocArray(4)~=1%MB DWI neonates
                                        if recAux.Par.Parameter2Read.extr2(end)>=NVS-1%The second one is full, so we discard this one --- note this is a sensitive change, it was like it is below!!!
                                        %if recAux.Par.Parameter2Read.extr2(end)==NVS-1
                                            rec.Rec=0;%We discard this one
                                        elseif recAux.Par.Parameter2Read.extr2(end)+1+curV+1>=NVS%They complete a study 
                                            rec.recV{s-n+1}=recAux;
                                            rec.Par.Mine.splitInfo(NVS-recAux.Par.Parameter2Read.extr2(end):NVS)=s-n+1;
                                            rec.recV{s-n+1}.Par.Mine.InitVol=NVS-recAux.Par.Parameter2Read.extr2(end)-1;
                                            shiftIsoc=recAux.Par.Mine.Isoc-rec.Par.Mine.Isoc;
                                            rec.Par.Mine.ShiftIsoc=shiftIsoc;
                                            fprintf('Series %s: %d vols\n',recAux.Names.Name,recAux.Par.Parameter2Read.extr2(end)+1);
                                            interveneNext=2;
                                        end
                                    end
                                    recAux=[];
                                end
                            end
                        else
                            s=n+1;
                            if s<=nV(end)
                                recAux.Names.Name=strtrim(prot.B_FileName(s,:));
                                matFile=fullfile(headFolder,sprintf('%s.mat',recAux.Names.Name));rawFile=fullfile(pathIn,sprintf('%s.raw',recAux.Names.Name));
                                recAux.Names.rawFile=rawFile;recAux.Names.matFile=matFile;recAux.Names.pathOu=pathOu;recAux.Names.pathIn=pathIn;recAux.Names.path=path;recAux.Names.prot=prot;recAux.Names.ind=n;recAux.Names.headFolder=headFolder;recAux.Names.versCode=versCode;recAux.Names.Specific=specific;
                                if exist(rawFile,'file')
                                    load(matFile);
                                    recAux.Par=Par;Par=[];
                                    if recAux.Par.Mine.Modal==10 && recAux.Par.Mine.AdHocArray(1)==101 && recAux.Par.Mine.AdHocArray(4)~=1%MB DWI neonates
                                        if recAux.Par.Parameter2Read.extr2(end)+1+curV+1>=NVS && recAux.Par.Parameter2Read.extr2(end)+1<NVS%They complete a study and second one is not full 
                                            interveneNext=1;
                                            initVol=NVS-recAux.Par.Parameter2Read.extr2(end)-1;
                                            shiftIsoc=recAux.Par.Mine.Isoc-rec.Par.Mine.Isoc;
                                        end
                                    end
                                end
                                recAux=[];
                            end
                        end
                    end
                end     
                if rec.Par.Mine.Modal==7;recRep{rep}=rec;end                                
                if rec.Rec && repeated<2;rec=reconPipeline(rec);end                  
                if rec.Par.Mine.Modal==7 
                    if isfield(rec,'Fail') && rec.Fail;recRep=recRep(1:rep-1);else rep=rep+1;end
                end
                if repeated<2
                    if (isfield(rec,'Fail') && rec.Fail) || (isfield(rec,'Plan') && any(rec.Plan.ReturnHeader))
                        if nargout>0;fRec{contF}=finalizeReconstruction(rec);end
                        contF=contF+1;
                        if rec.Alg.DetectAnomalies<=2 && rec.Anomaly.Detected(1);return;end              
                    end
                end
            end
        %end
    else
        interveneNext=0;
    end
    rec=[];
end

if repeated
    lNotUse=[];
    for n=rep-1:-1:2
        repStr=n;
        for l=n-1:-1:1
            if ~ismember(l,lNotUse) && checkRepeatedScans(recRep{n},recRep{l})
                repStr=[repStr l];
                lNotUse=[lNotUse l];
            end
        end
        if length(repStr)>1
            fprintf('Repeated scans:\n');
            for s=1:length(repStr)
                fprintf('%s\n',recRep{repStr(s)}.Names.Name);
            end            
            %KEY VALUES
            for s=2:length(repStr);recRep{n}.addRec{s-1}=recRep{repStr(s)};end
            rec=reconPipeline(recRep{n});
            if (isfield(rec,'Fail') && rec.Fail) || any(rec.Plan.ReturnHeader)
                if nargout>0;fRec{contF}=finalizeReconstruction(rec);end
                contF=contF+1;
            end
        end
    end
end

function rec=finalizeReconstruction(rec)

if isfield(rec,'Dyn')
    for m=setdiff(union(single(rec.Dyn.Typ2Rec'),rec.Plan.Typ2Rec'),rec.Plan.ReturnHeader)
        if isfield(rec,rec.Plan.Types{m});rec=rmfield(rec,rec.Plan.Types{m});end
    end
    if any(rec.Plan.ReturnHeader)
        for m=rec.Plan.ReturnHeader
            if isfield(rec,rec.Plan.Types{m});rec.(rec.Plan.Types{m})=gather(rec.(rec.Plan.Types{m}));end
        end
    end
end

if isfield(rec,'Corr');rec=rmfield(rec,'Corr');end%The problem with these is that they may be gpuArrays, it would be good to gather to avoid deleting
if isfield(rec,'Enc')
    if isfield(rec.Enc,'rGrid');rec.Enc=rmfield(rec.Enc,'rGrid');end%This is problematic due to similar reasons
    if isfield(rec.Enc,'rGridAcq');rec.Enc=rmfield(rec.Enc,'rGridAcq');end
    if isfield(rec.Enc,'kGrid');rec.Enc=rmfield(rec.Enc,'kGrid');end
    if isfield(rec.Enc,'DFTM');rec.Enc=rmfield(rec.Enc,'DFTM');end
    if isfield(rec.Enc,'DFTMH');rec.Enc=rmfield(rec.Enc,'DFTMH');end
    if isfield(rec.Enc,'DFTMAcq');rec.Enc=rmfield(rec.Enc,'DFTMAcq');end
    if isfield(rec.Enc,'DFTMHAcq');rec.Enc=rmfield(rec.Enc,'DFTMHAcq');end   
end

end

function isrep=checkRepeatedScans(recA,recB)
    isrep=1;
    scanValsStr={'ScanType','ScanMode','AcqMode','Technique','FastImgMode','PartialEcho','UTE','Kooshball','Orientation','FoldOverDir','FatShiftDir','ijk', ...
                 'MPS','xyz','REC','PatientPosition','PatientOrientation','Kt','FlowComp','Multivenc','PCAcqType','SPIR','MTC','Diffusion'};
    for ss=1:length(scanValsStr)
        if isfield(recA.Par.Scan,scanValsStr{ss}) && isfield(recB.Par.Scan,scanValsStr{ss}) && ~strcmp(recA.Par.Scan.(scanValsStr{ss}),recB.Par.Scan.(scanValsStr{ss}))
            isrep=0;break;
        end
        if xor(isfield(recA.Par.Scan,scanValsStr{ss}),isfield(recB.Par.Scan,scanValsStr{ss}))            
            isrep=0;break;
        end
    end
    scanValsVec={'TFEFactor','EPIFactor','HalfScanFactors','TR','TE','FlipAngle','FOV','curFOV','AcqVoxelSize','SliceGap','Samples','Stacks', ...
                 'ImagesPerStack','Offcentre','Angulation','MPSOffcentres','MPSOffcentresMM','xyzOffcentres','WaterFatShiftPix','SENSEFactor','KtFactor','Venc','kv','ASLNoLabelTypes'};
    for ss=1:length(scanValsVec)
        if isfield(recA.Par.Scan,scanValsVec{ss}) && isfield(recB.Par.Scan,scanValsVec{ss}) && (any(size(recA.Par.Scan.(scanValsVec{ss}))~=size(recB.Par.Scan.(scanValsVec{ss}))) || ~all(abs(recA.Par.Scan.(scanValsVec{ss})-recB.Par.Scan.(scanValsVec{ss}))<1e-3))
            isrep=0;break;
        end
        if xor(isfield(recA.Par.Scan,scanValsVec{ss}),isfield(recB.Par.Scan,scanValsVec{ss}))
            isrep=0;break;
        end
    end
    mineValsVec={'PhaseIncrement'};
    for ss=1:length(mineValsVec)
        if isfield(recA.Par.Mine,mineValsVec{ss}) && isfield(recB.Par.Mine,mineValsVec{ss}) && (any(size(recA.Par.Mine.(mineValsVec{ss}))~=size(recB.Par.Mine.(mineValsVec{ss}))) || ~all(recA.Par.Mine.(mineValsVec{ss})==recB.Par.Mine.(mineValsVec{ss})))
            isrep=0;break;
        end
        if xor(isfield(recA.Par.Mine,mineValsVec{ss}),isfield(recB.Par.Mine,mineValsVec{ss}))
            isrep=0;break;
        end
    end        
end

end