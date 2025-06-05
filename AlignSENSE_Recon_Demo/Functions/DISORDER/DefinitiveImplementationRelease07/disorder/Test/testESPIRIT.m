%load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/recN.mat');
%clearvars
tic
if ~exist('rec','var')
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/recESPIRIT.mat');
end
toc
%rec=recN;

rec.Alg.parE.NCV=2;%Maximum number of eigenmaps
rec.Alg.parE.NC=[3 3 3];%[2 2 2];%Resolution of calibration area to compute compression
rec.Alg.parE.K=[100 100 100];%Resolution of targeted coil profiles
rec.Alg.parE.eigTh=0;%Threshold for picking singular vectors of the calibration matrix (relative to the largest singular value)
rec.Alg.parE.absPh=0;%1%Flag to compute the absolute phase
rec.Alg.parE.virCo=1;%Flag to use the virtual coil to normalize the maps
rec.Alg.parE.eigSc=[0.85 0.85];%Cut-off for the eigenmaps in soft SENSE
rec.Alg.parE.subSp=[1 1 1];%Subsampling in image space to accelerate%This could perhaps be taken as parE.NC, but the parameter is kept for testing
rec.Alg.parE.mirr=[2 2 2];%Whether to mirror along a given dimension
rec.Alg.parE.dimLoc=[];%Dimensions along which to localize to compute virtual coils
rec.Alg.parE.Kmin=5;%Minimum K-value before mirroring
rec.Alg.parE.Ksph=200;%Number of points for spherical calibration area, unless 0, it overrides other K's

Enc=rec.Enc;
overDec=rec.Alg.OverDec;
overDec(overDec<0)=1;
Enc.FOVSize=round(Enc.FOVSize./overDec);
x=margosianFilter(rec.x,Enc,1,1);
rec.x=dynInd(rec.x,1,4);
rec.y=dynInd(rec.y,1,5);
rec.X=[];
rec.x=[];



recN=rec;
recN.Alg.solverType='IRWLS';%'IRWLS';%Solver type. One of the following: 'CG','IRWLS'

recN=solveX(recN);

visReconstruction(recN.x,0)
return
%recN.y=recN.y-10*recN.X;

recN.X=[];
recN.x=[];
recN=solveX(recN);
visReconstruction(recN.x,0)
return

%return


N=size(recN.y);
%recN.y=fold(recN.y,2,N(2),ceil(N(2)/2));
%recN.y=ifold(recN.y,2,N(2),ceil(N(2)/2));
%recN.x=recN.x;
%recN.x=fold(recN.x,2,N(2),ceil(N(2)/2));
%recN.x=ifold(recN.x,2,N(2),ceil(N(2)/2));
recN.y=recN.X;
x=recN.x;
return

visReconstruction(recN.x,0)
recN.x=sum(recN.y,4);
visReconstruction(recN.x,0)


recN=solveESPIRIT(recN);
recN.x=[];
recN.X=[];
recNN=solveX(recN);
visReconstruction(recNN.x,0)


%recN=solveX(recN);


%visReconstruction(dynInd(recN.x,[29 5],3:4),0)
%return

%recN.S=permute(recN.S,[2 3 1 4 5 6]);
%recN.W=permute(recN.W,[2 3 1 4 5 6]);

N=size(recN.S);

visReconstruction(cat(2,cat(1,dynInd(recN.S,[1 1],[3 6]),dynInd(recN.S,[1 2],[3 6])),cat(1,dynInd(recN.S,[floor(N(3)/2)+1 1],[3 6]),dynInd(recN.S,[floor(N(3)/2)+1 2],[3 6])),cat(1,dynInd(recN.S,[N(3) 1],[3 6]),dynInd(recN.S,[N(3) 2],[3 6]))),0)
visReconstruction(cat(2,cat(1,dynInd(recN.W,[1 1],[3 6]),dynInd(recN.W,[1 2],[3 6])),cat(1,dynInd(recN.W,[floor(N(3)/2)+1 1],[3 6]),dynInd(recN.W,[floor(N(3)/2)+1 2],[3 6])),cat(1,dynInd(recN.W,[N(3) 1],[3 6]),dynInd(recN.W,[N(3) 2],[3 6]))),0)






