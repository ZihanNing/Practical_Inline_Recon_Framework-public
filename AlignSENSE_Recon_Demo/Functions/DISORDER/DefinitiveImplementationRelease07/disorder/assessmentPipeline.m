function fAss=assessmentPipeline(path,modal,rootOu,series,specific)

%ASSESSMENTPIPELINE   Runs an assessment pipeline for a given modality
%   FPRO=ASSESSMENTPIPELINE(PATH,{MODAL},{ROOTOU},{SERIES},{SPECIFIC})
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
%   ** FASS returns the rec structures where the method failed
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

fAss=[];
contF=1;
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if nargin>=4 && ~isempty(series);nV=nV(ismember(nV,series));end
if nargin<5;specific=[];end
for n=nV
    rec.Names.Name=strtrim(prot.B_FileName(n,:));
    matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));
    rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;rec.Names.Specific=specific;
    if exist(matFile,'file')
        load(matFile);
        rec.Par=Par;Par=[];rec.Par.Mine.Proce=1;
        rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);
        if ~rec.Fail
            rec.Dyn.Typ2Wri(:)=0;modal=rec.Par.Mine.Modal;
            if modal==9 || modal==10%For other modalities it also run (see some examples in case 2017_12_18/LA_25710) but GPU problems appeared. TODO: free it to run for other modalities                
                fileRECON=fullfile(rec.Names.pathOu,numbe2Modal(modal),rec.Names.Name,'');
                suff=[];
                %THIS SHOULD PROBABLY BE ENCAPSULATED INTO A FUNCTION FOR
                %EACH ASSESSMENT TO BE PERFORMED---HERE, I THINK IT IS AN
                %EXAMPLE OF LOOKING AT THE FREQUENCY FIELDS
                suff{1}='Fr';
                existFile=1;
                for e=1:length(suff)
                    suff{e}=strcat(suff{e},rec.Plan.Suff);
                    if ~exist(strcat(fileRECON,'_',suff{e},'.mat'),'file');existFile=0;break;end
                end              
                if existFile               
                    fprintf('Reading %s\n',rec.Names.matFile);tsta=tic;
                    load(sprintf('%s_%s',fileRECON,suff{1}),'F');                    
                    recAux{1}=rec;rec=invertData(recAux,1);recAux=[];
                    rec.Dyn.Typ2Rec=23;rec.F=F;                   
                    tend=toc(tsta);fprintf('Time reading: %.3f s\n\n',tend);
                    if size(rec.F,4)==350
                        figure(1)
                        numplots=numel(get(gca,'Children'));
                        hold on
                        grid on
                        %if numplots<6
                        %    plot(rec.F(:),':','DisplayName',path)
                        if numplots<10
                            plot(rec.F(:),'--','DisplayName',path)
                        else
                            plot(rec.F(:),'DisplayName',path)
                        end
                        xlabel('Volume')
                        ylabel('Hz')
                        l=legend('-DynamicLegend');
                        set(l,'Interpreter', 'none');
                    end                    
                    if rec.Fail
                        fprintf('PROCESSING FAILED. LAST MESSAGE SHOULD PROVIDE THE REASON\n\n');  
                        if isfield(rec,'Dyn')
                            for m=rec.Dyn.Typ2Rec';rec.(rec.Plan.Types{m})=[];end            
                        end                    
                        fAss{contF}=rec;contF=contF+1;
                    end
                else
                    fprintf('Reconstruction data file %s not found\n',strcat(fileRECON,'_',suff{1},'.mat'));
                end
                %%if modal==10;break;end%The heuristic says that only the first one can be processed, as Anthony acquired the doubly reversed thing
            end    
        end
    end
    rec=[];
end
