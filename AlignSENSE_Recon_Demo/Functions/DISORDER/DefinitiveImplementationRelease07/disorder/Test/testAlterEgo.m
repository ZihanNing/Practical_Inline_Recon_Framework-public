addpath(genpath('/home/lcg13/Work/DefinitiveImplemenationDebug07'))

tic
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/groupSliceYeShear.mat')
toc

[e,Tn,d]=groupwiseSliceRegistration(y,M,TT,parU.Lambda,MB,[],parU.fractionOrder,[],parU.interpolTime(2));


Tn=permute(Tn,[5 4 6 1 2 3]);

for n=1:6
Tv=Tn(:,:,n);
figure
plot(Tv(:))
hold on
end

return

tic
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/groupSliceYeShear.mat')
toc


[e,Ty,d]=groupwiseSliceRegistration(y,M,TT,parU.Lambda,MB,[],parU.fractionOrder,[],parU.interpolTime(2));

Ty=permute(Ty,[5 4 6 1 2 3]);
Tv=Ty(:,:,6);
plot(Tv(:))



legend('No','Ye')

return

load('/home/lcg13/Data/rawDestin/ReconstructionsDebug07/SlavaNewReconstruction/2018_03_09/MA_621/Dy-Fu/ma_09032018_0955126_19_2_sf7fmb3r350sensep1o1V4_Tr.mat')

TT=permute(T,[5 4 6 1 2 3]);
TTT=TT(:,:,6);
figure;plot(TTT(:))

return

tic
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/spin.mat')
toc

H=spinHistoryEstimation(y,T,T1,TR,MB,FA,3);        
