path='2017_01_18/ge_193200';
modal=5;
rootOu='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06';
rootIn='/home/lcg13/Data/pnrawOrRelease03/archive-rawdata/archive-nnu';
series=[];

pathRF=fullfile(filesep,'home','lcg13','Work','MRecon-3.0.515');


addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(pathRF));

%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;

if nargin<2;modal=5;end%WE USE A SINGLE MODALITY, THE T1 CODE IS DEEMED AS EXPERIMENTAL
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end

%LOAD / GENERATE PROTOCOL INFO
protFile=fullfile(pathOu,'prot.txt');headFolder=fullfile(pathOu,'ZZ-HE');
if ~exist(protFile,'file');error('Protocol file %s not found',protFile);end
warning('off','MATLAB:namelengthmaxexceeded');prot=tdfread(protFile);warning('on','MATLAB:namelengthmaxexceeded');

fPro=[];
%TRAVERSE THROUGH THE DIFFERENT FILES IN THE SERIES TO BE RECONSTRUCTED
nV=find(ismember(prot.A_Modal,modal))';
if exist('series','var') && ~isempty(series);nV=nV(ismember(nV,series));end
for n=nV
    rec.Names.Name=strtrim(prot.B_FileName(n,:));
    matFile=fullfile(headFolder,sprintf('%s.mat',rec.Names.Name));
    rec.Names.matFile=matFile;rec.Names.pathOu=pathOu;rec.Names.prot=prot;rec.Names.ind=n;rec.Names.headFolder=headFolder;rec.Names.versCode=versCode;
    if exist(matFile,'file')
        load(matFile);
        rec.Par=Par;Par=[];rec.Par.Mine.Proce=1;
        rec.Fail=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);
        if ~rec.Fail
            rec.Dyn.Typ2Wri(:)=0;modal=rec.Par.Mine.Modal;
            if modal==5 || modal==6                                
                fileRECON=fullfile(rec.Names.pathOu,numbe2Modal(modal),rec.Names.Name,'');
                                       
                suff=[];
                suff{1}='Aq';
                cont=2;
                existFile=1;
                for e=1:length(suff)
                    suff{e}=strcat(suff{e},rec.Plan.Suff);
                    if ~exist(strcat(fileRECON,'_',suff{e},'.nii'),'file');existFile=0;break;end
                end                                        
                if existFile                                                    
                    fprintf('Reading %s\n',rec.Names.matFile);tsta=tic;                        
                    [y,MS,MT]=readNII(fileRECON,suff,rec.Dyn.GPU);
                    rawFile=fullfile(rootIn,path,rec.Names.Name);
                    
                    MR=MRecon(sprintf('%s.raw',rawFile));                    
                    rec.x=abs(y{1});
                    %MSw=MS{1};
                    %MTw=MT{1};
                    tend=toc(tsta);fprintf('Time reading: %.3f s\n\n',tend);                        
                    rec=geometryCorrection(rec,MR);
                    y{1}=rec.x;
                    writeNII(fileRECON,{'Gc'},y,MS,MT);
                    return
                else
                    fprintf('Reconstruction data file %s not found\n',strcat(fileRECON,'_',suff{1},'.nii'));
                end

                %%if modal==10;break;end%The heuristic says that only the first one can be processed, as Anthony acquired the doubly reversed thing
            end    
        end
    end
    rec=[];
end

gpu=rec.Dyn.GPU;
GeoCorrPars=rec.Par.Labels.GeoCorrPars;
NG=length(GeoCorrPars);
Spl=[2 (NG-18)/2 (NG-18)/2 16];
Spl=cumsum(Spl);

%r0=GeoCorrPars(1:Spl(1));assert(sum(r0~=0)==1,'Not 1 (%d) non-zero radious parameter',sum(r0~=0));r0=r0(r0~=0);
%for n=1:2;Gl{n}=GeoCorrPars(Spl(n)+1:Spl(n+1));assert(sum(G{n}~=0)~=25,'Not 25 (%d) non-zero Gl{%d} parameters',sum(Gl{n}~=0),n);Gl{n}=Gl{n}(Gl{n}~=0);end
%Gt=GeoCorrPars(Spl(3)+1:Spl(4));assert(sum(Gt~=0)==8,'Not 8 (%d) non-zero Gt parameters',sum(Gt~=0));Gt=Gt(Gt~=0);
%%Orders of the Legendre Polinomials (as deducted from observed parameters...)
%M=9;N=15;
%L=size(rec.x);

%for n=1:2:N
%    Gx

    
    