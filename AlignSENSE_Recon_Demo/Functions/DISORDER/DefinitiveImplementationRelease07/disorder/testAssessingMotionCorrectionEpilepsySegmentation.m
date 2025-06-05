
addpath(genpath(fileparts(mfilename('fullpath'))));
gpu=gpuDeviceCount;

% %IMAGE PATHS
% inpPath='/home/lcg13/Data/pnrawDe';
% inpFol='EpilepsyRelease07';
% pathor=studiesDetection(5);%Epilepsy cases
% path=pathor;path(8:11)=[];NCase=length(path);%Discard cases before 01/06/2019
% [user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
% modaNam={'mprage','t232','flair'};NModa=length(modaNam);
% recoNam={'Aq','Di','Re'};NReco=length(recoNam);
% %fileList=cell(NCase,NModa,NReco);
% fileList=[];
% for n=1:NCase
%     pathVe=strcat(inpPath,filesep,inpFol,filesep,path{n},filesep,'An-Ve');    
%     for m=1:NModa
%         for p=1:NReco
%             fileReco=strcat(pathVe,filesep,'*',modaNam{m},'*',recoNam{p},'.nii');
%             structReco=dir(fileReco);%Returns all the files and folders in the directory
%             nameReco={structReco.name};
%             if isempty(nameReco);error('File %s does not exist',fileReco);end
%             if strcmp(path{n},'2019_06_18/JA_106830') && m==1;sc=2;else sc=length(nameReco);end                  
%             nameReco=nameReco{sc};
%             fileList{n}{m}{p}=nameReco;
%         end
%     end
% end

%SEGMENTATION PATHS
%inpSegmPath='/home/lcg13/Data/pnrawDe';
%inpSegmFol='TestReconCode';

%            x=readNII(strcat(pathVe,filesep,fileList{n}{m}{p}(1:end-7)),recoNam,gpu);

fileList{1}{1}{1}='tu_01062019_0949441_6_2_mv3mpragegsttdisordersenseV4_Aq.nii';
fileList{1}{1}{2}='tu_01062019_0949441_6_2_mv3mpragegsttdisordersenseV4_Di.nii';
fileList{1}{1}{3}='tu_01062019_0949441_6_2_mv3mpragegsttdisordersenseV4_Re.nii';

inpSegmPath='/home/lcg13/Data';
inpPath='/home/lcg13/Data';
inpSegmFol='segmentation';
inpFol='segmentation';

pathSegm=strcat(inpSegmPath,filesep,inpSegmFol);
pathIma=strcat(inpPath,filesep,inpFol);
NCase=1;
segNam={'aq','di','re'};
recoNam={'Aq','Di','Re'};NReco=length(recoNam);
for n=1:NCase
    pathSegmCase=strcat(pathSegm,filesep,sprintf('seg%02d',n));
    pathCase=strcat(pathIma,filesep,sprintf('ima%02d',n));
    for p=1:NReco
        nii=load_untouch_nii(strcat(pathSegmCase,filesep,sprintf('aseg001%s.nii.gz',segNam{p})));
        MS=nii.hdr.dime.pixdim(2:4);
        MT=eye(4);MT(1,:)=nii.hdr.hist.srow_x;MT(2,:)=nii.hdr.hist.srow_y;MT(3,:)=nii.hdr.hist.srow_z;
        MS
        MT
        L{p}=nii.img;
        if gpu;L{p}=gpuArray(L{p});end
    end
    [I,MS,MT]=readNII(strcat(pathCase,filesep,fileList{n}{1}{p}(1:end-7)),recoNam,gpu);    
    for p=1:NReco
        MS{p}
        MT{p}
        %L{p}=shiftdim(L{p},1);
        L{p}=permute(L{p},[2 1 3]);
        %L{p}=flip(L{p},3);
        N=size(I{p});
        %L{p}=L{p}(size(L{p},1)-N(1)+1:end,:,:);
        %L{p}=L{p}(1:N(1),:,:);
        %L{p}=flip(L{p},1);
        L{p}=resampling(L{p},N,2);
    end
    visSegment(I{1},L{1}>0,0)
    return
    figure;imshow(I{1}(:,:,end/2),[])
    figure;imshow(L{1}(:,:,end/2)>0,[])
end

return
        


M=readtable(strcat(pathOu,filesep,'AE_01.xlsx'));
O1=M.(2);

load(strcat(pathOu,filesep,'QAData.mat'));
load(strcat(pathOu,filesep,'keyData.mat'));

modaNam={'mprage','t232','flair'};NModa=length(modaNam);
recoNam={'Aq','Di','Re'};NReco=length(recoNam);
NObse=1;
NCase=length(path);
scoreObs=zeros(NCase,NModa,NReco,NObse);

for o=1:NObse
    for l=1:length(O1)
        randObse=randObseV(o,l,:);
        scoreObs(randObse(1),randObse(2),randObse(3),o)=O1(l);
    end
end

for o=1:NObse
    figure
    co={'r','g','b'};
    for m=1:NModa        
        subtightplot(1,3,m);
        for p=1:NReco
            plot(convertRotation(T(:,m),'rad','deg'),scoreObs(:,m,p,o),strcat(co{p},'o'))
            hold on
        end
        xlabel('Motion (deg)')
        ylabel('Observer score')
        xlim([0 4])
        ylim([1 4])
        %lsline
        grid on
        title(modaNam{m});
        legend(recoNam);
        for p=1:NReco
            plot([min(convertRotation(T(:,m),'rad','deg')) max(convertRotation(T(:,m),'rad','deg'))],mean(scoreObs(:,m,p,o))*ones(1,2),strcat(co{p}))
        end
    end
    set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
end

