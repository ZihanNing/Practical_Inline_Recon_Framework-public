rawSource='/home/lcg13/Data/rawSource/pnraw/raw-nnu';
rawDestin='/home/lcg13/Data/rawDestin/ReconstructionsDebug06';

%path{1}='2018_08_16/MA_38430';%Prachi
%fil{1}='ma_16082018_1147282_6_2_fd1sensemtxV4';%Single echo
path{1}='2018_08_16/MA_38430';%Prachi
fil{1}='ma_16082018_1148360_7_2_fd1b1shimgclearp1o1V4';%Dual echo

%path{1}='2018_08_20/di_39531';%Santi
%fil{1}='di_20082018_1849291_4_2_wipethriveaxbhsenseV4';

%path{1}='2018_08_20/mu_39530';%Iñaki
%fil{1}='mu_20082018_1809237_4_2_wipdwi15sense1sh02senseV4';

%path{1}='2018_08_03/Bo_34530';%Jana
%fil{1}='bo_03082018_1653048_10_2_aslu4dtiaxsenseii1s9o76c1V4';




    
for n=1:length(fil)
    %matFil=fullfile(rawDestin,path,'ZZ-HE',sprintf('%s.mat',fil{n}));        
    %load(matFil);
    rawFill=fullfile(rawSource,path{n},sprintf('%s.raw',fil{n}));
    MR=MRecon(rawFill);
    Par=MR.Parameter;
    check=1;
    feedHeader(Par,check);
    %interpretLabels
    %feedHeader(Par,1);
    %feedHeader(Par,Names,1);
end
