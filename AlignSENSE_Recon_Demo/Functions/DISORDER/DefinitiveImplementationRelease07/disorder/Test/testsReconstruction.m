addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06'))


%pathData='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SlavaNewReconstruction/2018_06_13/LA_19930/Dy-Fu/la_13062018_1221378_17_2_fd1fsbsensep1o1V4Init';
pathData='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/SlavaNewReconstruction/2018_06_13/LA_19930/Dy-Fu/la_13062018_1221378_17_2_fd1fsbsensep1o1V4';

tic
if ~exist('rec','var');load(strcat(pathData,'001.mat'));end
toc
tic
rec.y=dynInd(rec.y,1,5);
toc
tic
recS=rec;
recS.Alg.parE.eigSc(1)=0.7;
recT=solveX(recS);
toc
visReconstruction(recT.x)
