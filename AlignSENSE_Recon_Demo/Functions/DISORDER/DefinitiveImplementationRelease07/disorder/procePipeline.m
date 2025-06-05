function fPro=procePipeline(path,modal,rootOu,series,specific)

%PROCEPIPELINE   Runs a preprocessing pipeline for a given modality
%   FPRO=PROCEPIPELINE(PATH,{MODAL},{ROOTOU},{SERIES},{SPECIFIC})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%   * {SERIES} restricts the reconstructions to a specific set of series
%   * {SPECIFIC} indicates to use a specific configuration of parameters as
%   stated in reconSpecific.m
%   * FPRO returns the rec structures where the method failed
%

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('\nPreprocessing study %s\n',path);
%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;

if nargin<2;modal=[];end
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end

%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');error('Protocol file %s not found',protFile);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

if nargout>0;fPro=[];end
contF=1;
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin<4;series=[];end
if ~isempty(series)
    if ~iscell(series);nV=nV(ismember(nV,series));else nV=nV(ismember(nV,series{1}));end
end
if nargin<5;specific=[];end

discardNext=0;
for n=nV%Traverse through series
    if iscell(series) && n~=nV(1);break;end
    if ~discardNext
        rec.Names.Name=strtrim(prot.B_FileName(n,:));
        matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));
        rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;rec.Names.Specific=specific;
        if exist(matFile,'file')
            load(matFile);
            rec.Par=Par;Par=[];rec.Par.Mine.Proce=1;
            modal=rec.Par.Mine.Modal;
            rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);rec.Dyn.Typ2Wri(:)=0;  
            if modal==10 && ~isfield(rec.Par.Mine,'StrFactorMax');rec.Alg.parU.useShimBox=0;rec.Alg.parU.corrMotion=0;end%Now we correct distortions
            suff=[];
            suff{1}='Aq';
            suff{2}='Ch';
            if ~rec.Fail
                if rec.Alg.SVDRecover;suff{3}='No';cont=4;                    
                else cont=3;
                end
                %For fMRI testing
                %suff{cont}='B0';suff{cont+1}='Un';%suff{cont+2}='Vo';
                %For DWI testing
                %suff{cont}='Re';
                existFile=0;            
                %MB DWI NEONATES (POTENTIAL SPLIT SCANS)
                if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') && modal==10 && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1 && rec.Par.Parameter2Read.extr2(end)~=size(rec.Par.Mine.diInfo,1)-1%Incomplete Split Scan, we read the following ones
                    Alg=reconAlgorithm;
                    NVS=size(rec.Par.Mine.diInfo,1);
                    curV=rec.Par.Parameter2Read.extr2(end);
                    [rec,fileRECON,existFile]=readDirectives(rec);
                    if Alg.JointSplitScans                                                                                                                                           
                        s=n+1;
                        if s<=nV(end)
                            recAux.Names.Name=strtrim(prot.B_FileName(s,:));
                            matFile=fullfile(headFolder,sprintf('%s.mat',recAux.Names.Name));
                            recAux.Names.matFile=matFile;recAux.Names.pathOu=pathOu;recAux.Names.prot=prot;recAux.Names.ind=n;recAux.Names.headFolder=headFolder;recAux.Names.versCode=versCode;recAux.Specific=specific;
                            if exist(matFile,'file')
                                load(matFile);
                                recAux.Par=Par;Par=[];
                                if recAux.Par.Mine.Modal==10 && recAux.Par.Mine.AdHocArray(1)==101 && recAux.Par.Mine.AdHocArray(4)~=1 && recAux.Par.Parameter2Read.extr2(end)+1+curV+1>=NVS && recAux.Par.Parameter2Read.extr2(end)+1<NVS%MB DWI neonates. They complete a study and second one is not full 
                                    discardNext=1;                             
                                end
                            end
                            recAux=[];
                        end
                    else
                        s=n+1;
                        if s<=nV(end)
                            recAux.Names.Name=strtrim(prot.B_FileName(s,:));
                            matFile=fullfile(headFolder,sprintf('%s.mat',recAux.Names.Name));
                            recAux.Names.matFile=matFile;recAux.Names.pathOu=pathOu;recAux.Names.prot=prot;recAux.Names.ind=n;recAux.Names.headFolder=headFolder;recAux.Names.versCode=versCode;recAux.Specific=specific;
                            if exist(matFile,'file')
                                load(matFile);
                                recAux.Par=Par;Par=[];recAux.Par.Mine.Proce=1;
                                recAux.Fail=0;recAux=reconPlanner(recAux);recAux=reconAlgorithm(recAux);recAux=reconSpecific(recAux);recAux.Dyn.Typ2Wri(:)=0;
                                if recAux.Par.Mine.Modal==10 && recAux.Par.Mine.AdHocArray(1)==101 && recAux.Par.Mine.AdHocArray(4)~=1 && recAux.Par.Parameter2Read.extr2(end)+1+curV+1>=NVS && recAux.Par.Parameter2Read.extr2(end)+1<NVS%MB DWI neonates. They complete a study and second one is not full                                                                         
                                    [recAux,~,isFile]=readDirectives(recAux);
                                    if isFile
                                        rec.Names.Name=strcat(rec.Names.Name,recAux.Names.Name);
                                        initVol=NVS-recAux.Par.Parameter2Read.extr2(end)-1;
                                        rec.x=cat(4,dynInd(rec.x,1:initVol,4),recAux.x);recAux.x=[];                                    
                                        if isfield(rec,'E');rec.E=cat(4,dynInd(rec.E,1:initVol,4),recAux.E);recAux.E=[];end                                    
                                        if rec.Alg.SVDRecover;rec.G=(rec.G*initVol+(NVS-initVol)*recAux.G)/NVS;recAux.G=[];end
                                        rec.Dyn.Typ2Wri([9:10 12])=1;
                                        if rec.Dyn.Debug>=1;fprintf('Writing %s\n',rec.Names.Name);tsta=tic;end
                                        writeData(rec);
                                        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time writing: %.3f s\n\n',tend);end
                                        rec.Dyn.Typ2Wri([9:10 12])=0;
                                        discardNext=1;                                  
                                    end
                                end
                                recAux=[];
                            end
                        end
                    end
                elseif (modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || modal==10%For other modalities it also run (see some examples in case 2017_12_18/LA_25710) but GPU problems appeared. TODO: free it to run for other modalities
                %elseif modal==9 || modal==10%For other modalities it also run (see some examples in case 2017_12_18/LA_25710) but GPU problems appeared. TODO: free it to run for other modalities
                    if iscell(series)
                        for ss=1:length(series{1})
                            rec.Names.Name=strtrim(prot.B_FileName(series{1}(ss),:));
                            [rec,fileRECON,existFile]=readDirectives(rec,ss);
                        end
                        rec.Names.Name=strcat(rec.Names.Name,'Concatenated');
                    else
                        [rec,fileRECON,existFile]=readDirectives(rec);
                    end
                end
            else
                fileRECON=fullfile(rec.Names.pathOu,numbe2Modal(modal),rec.Names.Name,'');
            end
            if existFile && (~strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') || size(rec.x,4)==size(rec.Par.Mine.diInfo,1) || modal~=10)     
                fprintf('Series %s\n',rec.Names.Name);
                if isfield(rec,'u');rec.Alg.SVDRecover=0;end
                if rec.Alg.SVDRecover && isempty(gcp('nocreate'));evalc('parpool');end%evalc is simply used to suppress warning messages
                
                    if rec.Plan.Quick%For testing                        
                        ND=size(rec.Enc.PartFourier,1)*ceil(60/size(rec.Enc.PartFourier,1));
                        ND=20;
                        rec.x=dynInd(rec.x,1:min(ND,size(rec.x,4)),4);
                        if isfield(rec,'E');rec.E=dynInd(rec.E,1:min(ND,size(rec.E,4)),4);end
                    end
                    
                    %THIS FUNCTION IS OVERRIDDEN-IT ONLY RETURNS ZEROS
                    if rec.Dyn.Debug>=1;fprintf('Stabilizing frequency %s\n',rec.Names.matFile);tsta=tic;end
                    rec=frequencyStabilization(rec);
                    if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time stabilizing frequency: %.3f s\n\n',tend);end
                    
                    if modal==10%DWI
                        if rec.Dyn.Debug>=1;fprintf('Masking %s\n',rec.Names.matFile);tsta=tic;end
                        rec=maskFromShim(rec);                
                        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time masking: %.3f s\n\n',tend);end
                        
                        if rec.Dyn.Debug>=1;fprintf('Filtering %s\n',rec.Names.matFile);tsta=tic;end
                        rec=SVDRecovery(rec);
                        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time filtering: %.3f s\n\n',tend);end
                    else%fMRI                        
                        %if rec.Dyn.Debug>=1;fprintf('Filtering %s\n',rec.Names.matFile);tsta=tic;end
                        %rec=SVDRecovery(rec);
                        %if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time filtering: %.3f s\n\n',tend);end
                        
                        if rec.Dyn.Debug>=1;fprintf('Tracking and detecting %s\n',rec.Names.matFile);tsta=tic;end
                        rec=brainTracking(rec);                    
                        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time tracking and detecting: %.3f s\n\n',tend);end
                    end 
                    
                    if rec.Dyn.Debug>=1;fprintf('Estimating distortion %s\n',rec.Names.matFile);tsta=tic;end
                    rec=estimateDistortion(rec);
                    if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time estimating distortion: %.3f s\n\n',tend);end
                    
                    if rec.Dyn.Debug>=1;fprintf('Correcting distortion %s\n',rec.Names.matFile);tsta=tic;end
                    rec=reverseDistortion(rec);
                    if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time correcting distortion: %.3f s\n\n',tend);end
                    
                    if rec.Dyn.Debug>=1;fprintf('Aligning %s\n',rec.Names.matFile);tsta=tic;end
                    rec=volumeAlignment(rec);
                    if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time aligning: %.3f s\n\n',tend);end

                    if ~rec.Fail
                        if rec.Dyn.Debug>=1;fprintf('Writing %s\n',rec.Names.matFile);tsta=tic;end
                        rec=writeData(rec);
                        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time writing: %.3f s\n\n',tend);end
                    else
                        fprintf('PROCESSING FAILED. LAST MESSAGE SHOULD PROVIDE THE REASON\n\n');  
                        if isfield(rec,'Dyn')
                            for m=rec.Dyn.Typ2Rec';rec.(rec.Plan.Types{m})=[];end            
                        end
                        %%TODO: CHECK IF THERE IS ANYTHING YET IN GPU
                        %%MEMORY BEFORE DOING THIS ASSIGNMENT                        
                        if nargout>0;fPro{contF}=rec;end
                        contF=contF+1;
                    end                    
            elseif ismember(modal,9:10)
                if rec.Dyn.Debug>=1;fprintf('Reconstruction data file %s not found/used\n',strcat(fileRECON,'_',suff{1},'.nii'));end
            end
        end
        rec=[];
    else
        discardNext=0;
    end
end

function [rec,fileRECONIn,existFile]=readDirectives(rec,ss)
    if nargin<2;ss=[];end    
    fileRECONIn=fullfile(rec.Names.pathOu,numbe2Modal(modal),rec.Names.Name,'');    
    existFile=1;
    for e=1:length(suff)
        suffIn{e}=strcat(suff{e},rec.Plan.Suff);
        if ~exist(strcat(fileRECONIn,'_',suffIn{e},'.nii'),'file');existFile=0;break;end
    end     
    rec.Dyn.Typ2Rec=[12;28];
    if (rec.Alg.SVDRecover && cont==3) || (~rec.Alg.SVDRecover && cont==2);rec.Dyn.Typ2Rec=12;end
    typ2Rec=rec.Dyn.Typ2Rec;
    if existFile               
        if rec.Dyn.Debug>=1;fprintf('Reading %s\n',rec.Names.Name);tsta=tic;end
        [y,MS,MT]=readNII(fileRECONIn,suffIn,0);        
        if isempty(ss) || ss==1
            recAuxInv{1}=rec;rec=invertData(recAuxInv,1);recAuxInv=[];
            for d=1:length(typ2Rec);rec.(rec.Plan.Types{typ2Rec(d)})=y{d};end         
            if rec.Alg.SVDRecover;rec.G=y{3};rec.Dyn.Typ2Rec;rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,9);end
            if length(suffIn)>=cont && modal==9;rec.B=y{cont};rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,11);end
            if length(suffIn)>=cont && modal==10;rec.r=y{cont};rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,22);end
            if length(suffIn)>=cont+1;rec.u=y{cont+1};rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,13);end
            if length(suffIn)>=cont+2;rec.v=y{cont+2};rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,14);end           
        else
            for d=1:length(typ2Rec);rec.(rec.Plan.Types{typ2Rec(d)})=cat(4,rec.(rec.Plan.Types{typ2Rec(d)}),y{d});end
        end
        y=[];
        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time reading: %.3f s\n\n',tend);end            
    end
end

end
