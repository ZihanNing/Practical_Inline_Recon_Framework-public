load('/home/lcg13/Work/DataDefinitiveImplementationDebug07/rec.mat');

recN=solveX(rec);
visReconstruction(recN.x)

recN=solveS(rec);
recNN=solveX(recN);
visReconstruction(recNN.x)

