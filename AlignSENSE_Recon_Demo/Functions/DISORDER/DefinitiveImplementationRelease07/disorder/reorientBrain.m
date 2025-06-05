function reorientBrain(path,modal,rootOu)

%REORIENTBRAIN   Reorients the brain using R Wright, B Khanal, A Gomez, E 
%Skelton, J Matthew, JV Hajnal, D Rueckert, JA Schnabel. LSTM Spatial 
%Co-transformer Networks for Registration of 3D Fetal US and MR Brain. 
%PIPPI/DATRA MICCAI. 2018. LNCS 11076:149-159.
%   REORIENTBRAIN(PATH,{MODAL},{ROOTOU})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%

addpath(genpath(fileparts(mfilename('fullpath'))));

%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;

if nargin<2;modal=5;end%WE USE A SINGLE MODALITY, THE T1 CODE IS DEEMED AS EXPERIMENTAL
if nargin<3 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);

if ~exist('gpu','var') || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end;

if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end

suff={'_Ex.nii','_Ma.nii','_RoAux.nii','_TrAux.mat'};
suffRead={'Ex','Ma','RoAux'};
suffWrite={'Ro','Mo'};

%suff={'_ExFull.nii','_MaFull.nii','_RoAuxFull.nii','_TrAuxFull.mat'};
%suffRead={'ExFull','MaFull','RoAuxFull'};
%suffWrite={'RoFull','MoFull'};

cnnPath=strcat(fileparts(mfilename('fullpath')),'/../../Software');
comm=sprintf('python %s/fetal-brain-pose-correction/fetal-brain-pose-correction.py',cnnPath);
PATH = getenv('PATH');
setenv('PATH', [sprintf('%s/miniconda2/bin:',cnnPath) PATH]);

%comm='python /home/lcg13/Software/fetal-brain-pose-correction/fetal-brain-pose-correction.py';
%PATH = getenv('PATH');
%setenv('PATH', ['/home/lcg13/Software/miniconda2/bin:' PATH]);

%WE LOOK FOR AN EX FILE
for m=1:length(modal)
    searchFolder=fullfile(rootOu,path,numbe2Modal(modal(m)));
    structId=dir(searchFolder);
    structId={structId.name};
    structId=structId(~cellfun(@isempty,regexp(structId,suff{1})));
    for s=1:length(structId)
        fileRecon=structId{s};
        fullFileRecon=fullfile(searchFolder,fileRecon);        
        system(sprintf('%s -i %s -m %s -o %s -t %s',comm,fullFileRecon,strrep(fullFileRecon,suff{1},suff{2}),strrep(fullFileRecon,suff{1},suff{3}),strrep(fullFileRecon,suff{1},suff{4})));     
        a=strsplit(fullFileRecon,'_','CollapseDelimiters',false);       
        fileRecon=strjoin(a(1:end-1),'_');
        [y,MS,MT]=readNII(fileRecon,suffRead,gpu);
        R=MT{3}(1:3,1:3);%NOTE THIS IS ONLY POSSIBLE BECAUSE THE RESOLUTION IS 1MM
        Tr=SpinCalc('DCMtoEA321',R,1e-3,0);
        Tr=-convertRotation(Tr,'deg','rad');
        Tr=reshape(Tr,[ones(1,3) 1 3]);
        Tr=repmat(Tr,[ones(1,4) 2]);
        Tr=dynInd(Tr,1:3,5,0);
        y(3)=[];MS(3)=[];MT(3)=[];
        for n=1:length(y)
            gpu=isa(y{n},'gpuArray');
            N=size(y{n});
            Nup=16*round(N(1:3)*sqrt(2)/16);           
            if n==2
               dist=4;%THIS IS A HARDCODED PARAMETER       
               y{n}=morphFourier(y{n},dist*ones(1,3),MS{2},ones(1,3),1);
            end              
            y{n}=resampling(y{n},Nup,2);        
            [~,kGrid,rkGrid,~,cGrid]=generateTransformGrids(Nup,gpu,Nup);
            et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,Tr,[],[],[],[],cGrid);
            y{n}=sincRigidTransform(y{n},et);
            y{n}=real(resampling(y{n},N(1:3),2));
            if n==1;y{n}=max(y{n},0);end
            if n==2;y{n}=single(y{n}>0.5);end
        end
        writeNII(fileRecon,suffWrite,y,MS,MT);
        delete(strrep(fullFileRecon,suff{1},suff{4}));
        delete(strrep(fullFileRecon,suff{1},suff{3}));
    end
end
