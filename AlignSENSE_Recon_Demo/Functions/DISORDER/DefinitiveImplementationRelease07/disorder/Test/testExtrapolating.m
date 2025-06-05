load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/Saux.mat');

visReconstruction(Saux,0)
M=abs(Saux)~=0;
visReconstruction(single(M),0)
dev=gpuDevice;
wait(dev);tic
SauxExtr=extrapolating(Saux,M,'PG');
wait(dev);toc
visReconstruction(SauxExtr,0)
