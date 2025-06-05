addpath(genpath(fileparts(mfilename('fullpath'))));


[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

paths=studiesDetection(2);

suffIn{1}='Ic';
suffIn{2}='Re';
suffIn{3}='El';
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

peff=1;
snr=zeros(200,6,2);
rotCor=zeros(200,6,3);
for p=1:length(paths)%Till 112 ready
    pathRIn=strcat('/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SVR/',paths{p},'/An-S2/');
    structPath = dir(pathRIn);%Returns all the files and folders in the directory
    structPath(ismember( {structPath.name}, {'.', '..'}))=[];
    structFil=structPath(~[structPath.isdir]);

    c=0;
    for n=1:length(structFil)        
        if strcmp(structFil(n).name(end-6:end),'_Ic.nii')%Just in case there is no B1 as in SE_13_08_2018
            c=c+1;
            candF{c}=structFil(n).name;    
        end
    end
    if c==6            
        for l=1:c
            fileRECONIn=strcat(pathRIn,candF{l}(1:end-7));
            fprintf(sprintf('%s\n',fileRECONIn));
            tic
            [y,MS,MT]=readNII(fileRECONIn,suffIn,gpu);
            toc
            R=MT{1}(1:3,1:3);
            R=bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));          
            R=convertRotation(acos(abs(R(2,2))),'rad','deg');
  
            y{3}=y{3}>0.5;            
            s=y{2}(y{3});
            n=(y{2}(y{3})-y{1}(y{3}));
            snrrank=(abs(s).^2)./(abs(n).^2);
            snrrank=sort(snrrank);
            snr(peff,l,1)=gather(mean(abs(s).^2)./mean(abs(n).^2));
            snr(peff,l,2)=gather(snrrank(round(0.25*length(snrrank))));
            snr(peff,l,3)=R;          
            
            %y{3}=x>0.5;
            %s=y{2}(y{3});
            %n=(y{2}(y{3})-y{1}(y{3}));
            %snr(peff,l,2)=gather(mean(abs(s).^2)./mean(abs(n).^2));                        
        end                                    
        peff=peff+1;
    end    
end
cas={'Mean','Rank25'};
snr=snr(1:peff-1,:,:);
conta=1;contb=1;

figure(1)
for n=1:2
    subplot(2,2,conta)
    boxplot(sqrt(snr(:,:,n)))
    ylabel('Ratio')    
    title(cas{n});     
    grid on
    conta=conta+1;
    subplot(2,2,conta)
    boxplot(10*log10(snr(:,:,n)))
    ylabel('dB')    
    title(cas{n})
    conta=conta+1;
    grid on
end
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])

L=[2 3 5 6];
xlimits=[min(min(snr(:,L,3))) max(max(snr(:,L,3)))];
for n=1:2
    ylimitsa(n,:)=[min(min(sqrt(snr(:,L,n)))) max(max(sqrt(snr(:,L,n))))];
    ylimitsb(n,:)=[min(min(10*log10(snr(:,L,n)))) max(max(10*log10(snr(:,L,n))))];
end

figure(2)   
for n=1:2
    for l=L
        subplot(4,4,contb)
        plot(snr(:,l,3),sqrt(snr(:,l,n)),'*')
        xlim(xlimits)
        ylim(ylimitsa(n,:))
        xlabel('Rotation (deg)')
        ylabel('Ratio')
        title(sprintf('%s - Stack %d',cas{n},l));
        contb=contb+1;
        grid on
    end    
    for l=L
        subplot(4,4,contb)
        plot(snr(:,l,3),10*log10(snr(:,l,n)),'*')
        xlim(xlimits)
        ylim(ylimitsb(n,:))
        xlabel('Rotation (deg)')
        ylabel('dB')        
        title(sprintf('%s - Stack %d',cas{n},l))
        contb=contb+1;
        grid on
    end
end
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
