load('/home/lcg13/Matlab/Data/dat.mat','S','B','AcqDelta','maskThS');

addpath(genpath('/home/lcg13/Matlab/DISORDER7TR02'));
[S,M]=solveSFilt(B,S,AcqDelta,maskThS);

any(isnan(M(:)))
return
