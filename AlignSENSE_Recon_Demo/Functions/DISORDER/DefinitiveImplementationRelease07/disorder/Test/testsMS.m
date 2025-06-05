addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'))

%return
%pathData='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SlavaNewReconstruction/2018_06_13/LA_19930/Dy-Fu/la_13062018_1221378_17_2_fd1fsbsensep1o1V4Init';
%pathData='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SlavaNewReconstruction/2018_06_13/LA_19930/Dy-Fu/la_13062018_1221378_17_2_fd1fsbsensep1o1V4';
pathData='/home/lcg13/Data/pnrawDe/ReconstructionsDebug07/2015_02_27/NG_12501/An-S2/ng_27022015_1406006_5_2_dhcp4t2tsesenseV4';
if ~exist('rec','var')
    tic
    load(strcat(pathData,'001.mat'));
    toc
    recS=rec;
end

rec=recS;
rec.Alg.parXT.threeD=1;
rec.Alg.parXT.alpha=[0.2 1];%Regularization
rec.Alg.parXT.correct=2;%To correct for motion
rec.Alg.parXT.threeD=1;%To deconvolve
rec.Alg.parXT.outlP=1.2;%To discard outliers

recT=solveXTMS(rec);%recT=[];

for n=0:2;visReconstruction(shiftdim(recT.x,n),0);end
for n=0:2;visReconstruction(shiftdim(recT.d,n),0);end

xx=cat(2,recT.x,recT.d);

for n=1:125
    visReconstruction(dynInd(xx,n,3))
    close all
end


%recS=solveXTMS(rec);

