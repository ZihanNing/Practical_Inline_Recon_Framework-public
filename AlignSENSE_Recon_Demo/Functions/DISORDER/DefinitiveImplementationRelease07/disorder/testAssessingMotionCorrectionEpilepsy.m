
addpath(genpath(fileparts(mfilename('fullpath'))));

inpPath='/home/lcg13/Data/pnrawDe';
inpFol='EpilepsyRelease07';
outFol='ES';

%WE GENERATE THE EPILEPSY CASES
pathor=studiesDetection(5);%Epilepsy cases
path=pathor;path(8:11)=[];NCase=length(path);
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
modaNam={'mprage','t232','flair'};NModa=length(modaNam);
recoNam={'Aq','Di','Re'};NReco=length(recoNam);
%fileList=cell(NCase,NModa,NReco);
fileList=[];
for n=1:NCase
    pathVe=strcat(inpPath,filesep,inpFol,filesep,path{n},filesep,'An-Ve');    
    for m=1:NModa
        for p=1:NReco
            fileReco=strcat(pathVe,filesep,'*',modaNam{m},'*',recoNam{p},'.nii');
            structReco=dir(fileReco);%Returns all the files and folders in the directory
            nameReco={structReco.name};
            if isempty(nameReco);error('File %s does not exist',fileReco);end
            if strcmp(path{n},'2019_06_18/JA_106830') && m==1;sc=2;else sc=length(nameReco);end                  
            nameReco=nameReco{sc};
            fileList{n}{m}{p}=nameReco;
        end
    end
end
pathOu=strcat(inpPath,filesep,outFol);

generateLinks=0;
if generateLinks
    NObse=4;
    N=NReco*NModa*NCase;
    randObseI=zeros(4,N);
    randObseV=zeros([4 N 3]);%Case/Modality/Reconstruction
    for n=1:NObse
        randObseI(n,:)=randperm(N);
        randObseV(n,:,:)=shiftdim(ind2subV([NCase NModa NReco],randObseI(n,:)),-1);
        outFolObse=strcat(inpPath,filesep,outFol,filesep,sprintf('O%d',n));
        if ~exist(outFolObse,'dir');mkdir(outFolObse);end
        delete(strcat(outFolObse,filesep,'*'));
        for r=1:N      
            outFolObseLink=strcat(outFolObse,filesep,sprintf('%04d.nii',r));
            indObse=randObseV(n,r,:);indObse=indObse(:);
            pathVe=strcat(inpPath,filesep,inpFol,filesep,path{indObse(1)},filesep,'An-Ve');            
            inpFil=strcat(pathVe,filesep,fileList{indObse(1)}{indObse(2)}{indObse(3)});
            system(['cp ' inpFil ' ' outFolObseLink]);
            %system(['ln -s ' inpFil ' ' outFolObseLink]);
        end
    end
    save(strcat(pathOu,filesep,'keyData.mat'),'path','fileList','randObseV');
end

computeMetric=1;
NDiscardShots=3;
if computeMetric==1
    %waveSp=zeros(NCase,NModa,NReco);
    gradEn=zeros(NCase,NModa,NReco);
    T=zeros(NCase,NModa);
    gpu=1;
    for n=1:NCase
        pathVe=strcat(inpPath,filesep,inpFol,filesep,path{n},filesep,'An-Ve');    
        for m=1:NModa        
            x=readNII(strcat(pathVe,filesep,fileList{n}{m}{p}(1:end-7)),recoNam,gpu);
            Td=load(strcat(pathVe,filesep,fileList{n}{m}{p}(1:end-7),'_Tr.mat'));
            Td=squeeze(Td.T);            
            Td=Td(:,4:6);
            Taux=sqrt(normm(Td,[],2));
            [~,iaux]=sort(Taux,'descend');
            Td(iaux(1:NDiscardShots),:)=[];
            Td=std(Td,[],1);
            Td=sqrt(normm(Td));
            T(n,m)=Td;            
            for p=1:NReco
                %x{p}=abs(x{p});
                %x{p}=x{p}(floor(end/6)+1:2*floor(end/3),:,:);
                gradEn(n,m,p)=gradientEntropy(x{p});
                %waveSp(n,m,p)=waveletSparsity(x{p});%Note this does not
                %work with gpu arrays, seems highly sensitive to the used
                %area
            end
        end
    end
    %save(strcat(pathOu,filesep,'QAData.mat'),'path','gradEn','waveSp');
    save(strcat(pathOu,filesep,'QAData.mat'),'path','gradEn','T');
elseif computeMetric>=2
    load(strcat(pathOu,filesep,'QAData.mat'))
    if computeMetric==3
        for n=1:NCase
            pathVe=strcat(inpPath,filesep,inpFol,filesep,path{n},filesep,'An-Ve');    
            for m=1:NModa        
                Td=load(strcat(pathVe,filesep,fileList{n}{m}{p}(1:end-7),'_Tr.mat'));
                Td=squeeze(Td.T);            
                Td=Td(:,4:6);
                Taux=sqrt(normm(Td,[],2));
                [~,iaux]=sort(Taux,'descend');
                Td(iaux(1:NDiscardShots),:)=[];
                Td=std(Td,[],1);
                Td=sqrt(normm(Td));
                T(n,m)=Td;    
            end
        end
    end
    figure
    co={'r','g','b'};
    for m=1:NModa        
        subtightplot(1,3,m);
        for p=1:NReco
            plot(convertRotation(T(:,m),'rad','deg'),gradEn(:,m,p),strcat(co{p},'o'))        
            hold on
        end
        xlabel('Motion (deg)')
        ylabel('Gradient entropy')
        xlim([0 4])
        ylim([0.9 0.97])
        lsline
        grid on
        title(modaNam{m});
        legend(recoNam);
    end
    set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
end
