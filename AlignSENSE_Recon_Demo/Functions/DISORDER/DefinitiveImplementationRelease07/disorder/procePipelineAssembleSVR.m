function [fPro,svr]=procePipelineAssembleSVR(path,modal,rootOu,series,specific)

%PROCEPIPELINE   Runs a preprocessing pipeline for a given modality
%   [FPRO,SVR]=PROCEPIPELINE(PATH,{MODAL},{ROOTOU},{SERIES},{SPECIFIC})
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
%   ** FPRO returns the rec structures where the method failed
%   ** SVR returns the SVR structure
%

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('\nAssembling study %s\n',path);
%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;
if nargin<2 || isempty(modal);modal=5:6;end%WE USE A SINGLE MODALITY, THE T1 CODE IS DEEMED AS EXPERIMENTAL
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end
modal(~ismember(modal,5:6))=[];

%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');error('Protocol file %s not found',protFile);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

%PATH OF SVR CODE
svrPath=strcat(fileparts(mfilename('fullpath')),'/../../Software');

%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
if nargout>0;fPro=[];end
for moda=modal
    svr=[];
    csvr=1;

    nV=find(ismember(prot.A_Modal,moda))';
    if nargin>=4 && ~isempty(series);nV=nV(ismember(nV,series));end
    if nargin<5;specific=[];end

    for n=nV
        rec.Names.Name=strtrim(prot.B_FileName(n,:));
        matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));
        rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;rec.Names.Specific=specific;
        if exist(matFile,'file')
            load(matFile);
            rec.Par=Par;Par=[];rec.Par.Mine.Proce=1;
            if strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra')
                rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);
                if ~rec.Fail                    
                    rec.Dyn.Typ2Wri(:)=0;modaEff=rec.Par.Mine.Modal;
                    if modaEff==5 || modaEff==6
                        fprintf('Series %s\n',rec.Names.Name);
                        fileRECON=fullfile(rec.Names.pathOu,numbe2Modal(modaEff),rec.Names.Name);

                        suff=[];
                        suff{1}='Re';
                        if ~isempty(rec.Alg.parXT.outputResol);suff{1}=strcat(suff{1},sprintf('%.2f',rec.Alg.parXT.outputResol));end
                        suffAux='SurrVoluReco';%FWHM4';
                        if ~isempty(rec.Alg.parXT.outputResol);suffAux=strcat(suffAux,sprintf('%.2f',rec.Alg.parXT.outputResol));end
                        existFile=1;
                        for e=1:length(suff)
                            suff{e}=strcat(suff{e},rec.Plan.Suff);
                            if ~exist(strcat(fileRECON,'_',suff{e},'.nii'),'file');existFile=0;break;end
                        end                       
                        if existFile                       
                            [~,MS,MT]=readNII(fileRECON,suff,rec.Dyn.GPU);
                            if rec.Dyn.Debug>=2
                                fprintf('Series: %s\n',rec.Names.Name);
                                fprintf('Resolution:%s\n',sprintf(' %.2f',MS{1}));
                                fprintf('Acquired resolution:%s\n',sprintf(' %.2f',rec.Par.Scan.AcqVoxelSize));
                                fprintf('Slice gap:%s\n',sprintf(' %.2f',rec.Par.Scan.SliceGap(1)));
                                fprintf('Scan duration:%s\n',sprintf(' %.2f',rec.Par.Labels.ScanDuration(1)));
                                fprintf('Scan orientation: %s\n',rec.Par.Scan.Orientation);
                            end
                            svr.fileRECON{csvr}=sprintf('%s_%s.nii',fileRECON,suff{1});
                            svr.rec{csvr}=rec;
                            csvr=csvr+1;
                        else
                            fprintf('Reconstruction data file %s not found\n',strcat(fileRECON,'_',suff{1},'.nii'));
                        end
                    end   
                end
            end
        end
        rec=[];
    end    

    if ~isfield(svr,'rec');NV=0;else NV=length(svr.rec);end
    if NV>=2   
        %SAGITTAL AS PRIMARY VIEW (FOR LARGEST FOV)
        perm=1:NV;
        nSag=[];
        for n=NV:-1:1
            if strcmp(svr.rec{n}.Par.Scan.Orientation,'SAG');nSag=n;break;end
        end
        if ~isempty(nSag);perm([1 nSag])=[nSag 1];end
        svr.fileRECON=svr.fileRECON(perm);
        svr.rec=svr.rec(perm);
        
        resol=svr.rec{1}.Alg.parXT.SVRoutputResol;
        
        %PARAMETERS OF CALL
        fileOut=strsplit(svr.fileRECON{1},'_');
        fileOut=strcat(strjoin(fileOut(1:end-1),'_'),sprintf('_Ro%.2f',resol));
        if ~isempty(svr.rec{1}.Alg.parXT.outputResol);fileOut=strcat(fileOut,sprintf('%.2f',svr.rec{1}.Alg.parXT.outputResol));end
        fileOut=strcat(fileOut,'.nii');
        
        recInp='';
        pack='';
        for n=1:NV
            recInp=[recInp svr.fileRECON{n} ' '];         
            pack=[pack ' 2'];
        end
        
        %SYSTEM CALL
        %Not sure if to use nohup here, activate in case of troubles
        outFol=[svr.rec{1}.Names.pathOu '/' numbe2Modal(modaEff) '/' suffAux '/'];
        curFol=pwd;
        if ~exist(outFol,'dir');mkdir(outFol);end
        cd(outFol);        
        if svr.rec{1}.Dyn.Debug>=2
            systCall=sprintf('time %s/irtk-build-release-multith-notbb/bin/reconstructionNeo %s %d %s-packages%s -resolution %.2f -remove_black_background -no_robust_statistics -iterations 3 -log_prefix %s',svrPath,fileOut,NV,recInp,pack,resol,outFol);
            fprintf('System call: %s\n',systCall);
        else
            systCall=sprintf('%s/irtk-build-release-multith-notbb/bin/reconstructionNeo %s %d %s-packages%s -resolution %.2f -remove_black_background -no_robust_statistics -iterations 3 -log_prefix %s > /dev/null',svrPath,fileOut,NV,recInp,pack,resol,outFol);
        end            
        system(systCall);
        outFol=[svr.rec{1}.Names.pathOu '/' numbe2Modal(modaEff) '/'];
        cd (outFol);
        rmdir(suffAux,'s');
        cd(curFol)
    end
end
