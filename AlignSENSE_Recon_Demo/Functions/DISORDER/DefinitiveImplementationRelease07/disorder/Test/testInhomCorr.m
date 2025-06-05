%fil='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/2018_03_23/WI_4221/An-S2/wi_23032018_0946056_9_2_sf7zt2mbzax0sensemi8s1V4';
fil='/home/lcg13/Work/DataDefinitiveImplementationDebug06/wi_23032018_0946056_9_2_sf7zt2mbzax0sensemi8s1V4';


gpu=0;

suff={'SVR_After0Mot2.00'};
[y,MS,MT]=readNII(fil,suff,gpu);
x=y{1};
b=biasFieldEstimation(x,1);
y{1}=b;
suff={'SVR_BiasField'};
writeNII(fil,suff,y,MS,MT);
y{1}=real(x./b);%This does not work well
suff={'SVR_BiasFieldCorrected'};
writeNII(fil,suff,y,MS,MT);
