
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

paths=studiesDetection(2);

suffIn{1}='Aq';
gpu=1;

%14->2018_03_23/WI_4221->Slight problem
%49->2018_05_18/MA_12831->Failed!!!
%50->2018_05_21/PA_13130->Failed!!!
%51->2018_05_23/HI_14030->Slight problem
%75->2018_06_25/MC_22331->Slight problem
%85->2018_07_11/HU_26231->Problem
%paths{1}='2018_05_18/MA_12831';%49
%paths{1}='2018_05_21/PA_13130';%50
%paths{1}='2018_03_23/WI_4221';%14
%paths{1}='2018_05_23/HI_14030';%51
%paths{1}='2018_06_25/MC_22331';%75
%paths{2}='2018_07_11/HU_26231';%85

%Check for rings:2018_07_04/SA_24430(82)
%paths{1}='2018_07_04/SA_24430';

for p=114:length(paths)%Till 112 ready
    pathRIn=strcat('/home/lcg13/Data/pnrawDe/ReconstructionsRelease06/',paths{p},'/Dy-Di/');
    pathROu=strcat('/home/lcg13/Data/pnrawDe/ReconstructionsRelease06SVR/',paths{p},'/Dy-Di/');
    structPath = dir(pathRIn);%Returns all the files and folders in the directory
    structPath(ismember( {structPath.name}, {'.', '..'}))=[];
    structFil=structPath(~[structPath.isdir]);

    c=0;
    for n=1:length(structFil)        
        if strcmp(structFil(n).name(end-6:end),'_Aq.nii')
            c=c+1;
            candF{c}=structFil(n).name;            
        end
    end
    if c>0        
        fileRECONIn=strcat(pathRIn,candF{c}(1:end-7));
        fileRECONOu=strcat(pathROu,candF{c}(1:end-7));
        fprintf(sprintf('%s\n',fileRECONIn));
        tic
        [y,MS,MT]=readNII(fileRECONIn,suffIn,gpu);
        toc
        
        if size(y{1},4)==282       

            %FEATURES
            tic
            direcs=load('/home/lcg13/Work/DefinitiveImplementationDebug06/Data/DynamicFiles/dhcp-fetal-final.txt');
            b=direcs(:,4);
            bun=unique(b);
            bmaxind=find(b==max(bun));
            y=y{1};
            y=cat(4,dynInd(y,2*bmaxind-1,4),2*dynInd(y,2*bmaxind,4));
            
            x=ridgeDetection(y,1:2,1);
            y=bsxfun(@times,y,conj(x));
            NY=size(y);NY(end+1:4)=1;
            y=reshape(y,[NY(1:3) NY(4)/2 2]);
            y=mean(y,4);
            y=abs(y);
            y=mean(y,5);
            z{1}=y;
            suff{1}='AveBMax';
            toc

            %ROUGH MASK
            tic
            Alg.parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
            Alg.parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
            Alg.parS.Otsu=[0 1];%Binary vector of components for multilevel extraction (it picks those components with 1)
            Alg.parS.nErode=6;%Erosion for masking (in mm)
            Alg.parS.nDilate=16;%Dilation for masking (in mm)
            Alg.parS.conComp=1;%Whether to get the largest connected component after erosion
            voxsiz=MS{1};
            Mla=refineMask(y,Alg.parS,voxsiz);
            z{2}=Mla;
            suff{2}='AveBMaxMaskFirst';
            toc

            %ELLIPSOID FIT
            tic
            [par,Mlb]=ellipsoidFromImage(Mla,[],0);
            %[Mlb,sphcen0,sphrad]=ellipsoidalHoughTransform(y,voxsiz,[20 60],[],[],Mla);
            z{3}=Mlb;
            suff{3}='AveBMaxMaskSecond';
            toc

            tic
            NS=length(suff);
            MSV=cell(1,NS);
            MTV=cell(1,NS);
            for n=1:NS;MSV{n}=MS{1};MTV{n}=MT{1};end
            writeNII(fileRECONOu,suff,z,MSV,MTV);
            toc
        end
    end
end



