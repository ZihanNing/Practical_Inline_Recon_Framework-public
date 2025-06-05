
addpath(genpath(fileparts(mfilename('fullpath'))));

inpPath='/home/lcg13/Data/pnrawDe';
inpFol='EpilepsyRelease07';
outFol='ES';
pathOu=strcat(inpPath,filesep,outFol);

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

