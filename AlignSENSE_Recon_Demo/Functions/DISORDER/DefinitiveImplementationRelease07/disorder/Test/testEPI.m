tic
if ~exist('rec','var')
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/rec.mat');
end
toc

rec.Alg.parE.eigSc=[0.85 0.3];%0.25
rec.Alg.UnfoldFat=2;

recN=solveX(rec);
y(:,:,1,1)=dynInd(recN.x,[20 7],3:4);
y(:,:,2,1)=dynInd(recN.x,[33 7],3:4);
y(:,:,1,2)=dynInd(recN.x,[23 8],3:4);
y(:,:,2,2)=dynInd(recN.x,[24 8],3:4);
y=y./multDimMax(y,1:2);

visReconstruction(cat(2,cat(1,y(:,:,1,1),y(:,:,2,1)),cat(1,y(:,:,1,2),y(:,:,2,2))),0)

return

gpuIn=1;

N=size(recN.x);
caux=log(dynInd(recN.X,7,4));
M=buildFilter(2*N(1:3),'tukeyIso',1,gpuIn,1,1);
caux=filtering(caux,M,1);
cauxa=2*(caux-median(caux))./(1.5*median(caux))-1;
cauxa=0.5+0.5*tanh(3*cauxa);%Standard sigmoid---1 indicates almost sure artifact, 0 indicates no artifact
cauxa=abs(filtering(cauxa,M,1));

visReconstruction(dynInd(cauxa,20,3),0);
visReconstruction(dynInd(cauxa,33,3),0);

cauxa=2*(caux-median(caux))./(3*median(caux))-1;
cauxa=0.5+0.5*tanh(3*cauxa);%Standard sigmoid---1 indicates almost sure artifact, 0 indicates no artifact
cauxa=abs(filtering(cauxa,M,1));

visReconstruction(dynInd(cauxa,20,3),0);
visReconstruction(dynInd(cauxa,33,3),0);




return

Wauxa=(1-cauxa)*rec.Alg.parE.eigSc(1)+cauxa*rec.Alg.parE.eigSc(2);                        
Waux=2*bsxfun(@times,bsxfun(@minus,E.WW,Wauxa),1./(1-Wauxa))-1;
Waux=0.5+0.5*tanh(3*Waux);%"Standard" sigmoid                      


