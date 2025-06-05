function fPro=procePipelineTransform(path,modal,rootOu,series,suffTr,shiftOverlap,specific)

%PROCEPIPELINETRANSFORM   Applies the transforms of a preprocessing
%pipeline to a different piece of data
%   FPRO=PROCEPIPELINETRANSFORM(PATH,{MODAL},{ROOTOU},{SERIES},{SPECIFIC})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions)
%   * {SERIES} restricts the reconstructions to a specific set of series
%   * {SUFFTR} is the suffix to transform
%   * {SPECIFIC} indicates to use a specific configuration of parameters as
%   stated in reconSpecific.m
%   * FPRO returns the rec structures where the method failed
%

addpath(genpath(fileparts(mfilename('fullpath'))));

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

fPro=[];
contF=1;
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin<4;series=[];end
if ~isempty(series)
    if ~iscell(series);nV=nV(ismember(nV,series));else nV=nV(ismember(nV,series{1}));end
end
if nargin<5;suffTr='Aq';end%Suffix to transform
if nargin<6;shiftOverlap=[];end
if nargin<7;specific=[];end

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
            suff{1}=suffTr;
            suff{2}='Ma';
            if ~rec.Fail
                existFile=0;            
                %MB DWI NEONATES (POTENTIAL SPLIT SCANS)---THIS IS
                %PROBLEMATIC AT THE MOMENT, JUST LEFT FOR FUTURE
                %EXTENSION...
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
                                        rec.E=cat(4,dynInd(rec.E,1:initVol,4),recAux.E);recAux.E=[];                                    
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
                    if modal==10 && ~isfield(rec.Par.Mine,'StrFactorMax');rec.Alg.parU.useShimBox=0;rec.Alg.parU.corrMotion=0;end%Now we correct distortions
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
            if existFile && (~strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') || size(rec.E,4)==size(rec.Par.Mine.diInfo,1) || modal~=10)     
                if isfield(rec,'u');rec.Alg.SVDRecover=0;end
                if rec.Alg.SVDRecover && isempty(gcp('nocreate'));parpool;end                
                    if rec.Plan.Quick%For testing                        
                        ND=size(rec.Enc.PartFourier,1)*ceil(60/size(rec.Enc.PartFourier,1));
                        ND=10;
                        rec.E=dynInd(rec.E,1:min(ND,size(rec.E,4)),4);
                    end
                    
                    if modal==9%FOR DWI NOT IMPLEMENTED                       
                        if ~isempty(shiftOverlap)
                            rec=mapOverlappedVoxels(rec,shiftOverlap);
                            suffTrO=strcat(suffTr,num2str(shiftOverlap));                     
                        else
                            suffTrO=suffTr;
                        end
                        
                        suffFil={'Fr','Fo','B0','Tr'};%Frequency-Tracking-Distortion-Motion
                        suffTyp={'F','t','B','T'};%Idem
                        suffFor={'.mat','.mat','.nii','.mat'};%Idem
                        suffStr={'Stabilizing frequency','Tracking and detecting','Correcting distortion','Aligning'};                                                
                        
                        for t=1:length(suffFil) 
                            fileFreq=fullfile(rec.Names.pathOu,numbe2Modal(modal),sprintf('%s_%s%s',rec.Names.Name,suffFil{t},suffFor{t}));
                            if ~exist(fileFreq,'file');fprintf('%s (%s) not found\n',suffStr{t},fileFreq);
                            else
                                if strcmp(suffFor{t},'.mat');recAux=load(fileFreq,suffTyp{t});rec.(suffTyp{t})=recAux.(suffTyp{t});recAux=[];
                                else y=readNII(fileRECON,suffFil(t),0);rec.(suffTyp{t})=y{1};
                                end
                                
                                if rec.Dyn.Debug>=1;fprintf('%s %s\n',suffStr{t},rec.Names.matFile);tsta=tic;end
                                rec
                                if t==1
                                    %rec.F(:)=0;
                                    rec=frequencyStabilization(rec,1);
                                elseif t==2
                                    %for
                                    %k=1:length(rec.t);rec.t{k}(:)=0;end%To check consistency
                                    rec=brainTracking(rec,1);
                                elseif t==3
                                    %rec.B(:)=0;%To check consistency
                                    if strcmp(suffTr,'Sh');norm=1;else norm=0;end
                                    rec=reverseDistortion(rec,rec.Par.Labels.TE,norm,[],1);
                                elseif t==4
                                    %rec.T(:)=0;%To check consistency
                                    rec=volumeAlignment(rec,1);
                                end
                                if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time %s: %.3f s\n\n',suffStr{t},tend);end                                             
                            end
                            rec
                        end
                        if ~rec.Fail
                            for t=1:length(rec.Plan.TypeNames);rec.Plan.TypeNames{t}=strcat(rec.Plan.TypeNames{t},suffTrO);end                            
                            rec.Dyn.Typ2Wri(28)=1;
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
                            fPro{contF}=rec;contF=contF+1;
                        end       
                    end
            elseif ismember(modal,9:10)
                fprintf('Reconstruction data file %s not found\n',strcat(fileRECON,'_',suff{1},'.nii'));
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
    rec.Dyn.Typ2Rec=[28;8];
    typ2Rec=rec.Dyn.Typ2Rec;
    if existFile
        if rec.Dyn.Debug>=1;fprintf('Reading %s\n',rec.Names.Name);tsta=tic;end
        [y,MS,MT]=readNII(fileRECONIn,suffIn,0);        
        if isempty(ss) || ss==1
            recAuxInv{1}=rec;rec=invertData(recAuxInv,1);recAuxInv=[];
            for d=1:length(typ2Rec);rec.(rec.Plan.Types{typ2Rec(d)})=y{d};end     
        else
            for d=1:length(typ2Rec);rec.(rec.Plan.Types{typ2Rec(d)})=cat(4,rec.(rec.Plan.Types{typ2Rec(d)}),y{d});end
        end
        rec.MS=MS;
        rec.MT=MT;        
        y=[];
        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time reading: %.3f s\n\n',tend);end            
    end
end

end
