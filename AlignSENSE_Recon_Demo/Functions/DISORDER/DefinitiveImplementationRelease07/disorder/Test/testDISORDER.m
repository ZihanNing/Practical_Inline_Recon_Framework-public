
tic
if ~exist('CS','var')
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/CS.mat');
end
toc

%TO ACCELERATE
%NC=min(4,size(CS.E.Sf,4));
%CS.E.Sf=dynInd(CS.E.Sf,1:NC,4);CS.E.dS(2)=NC;
%CS.y=dynInd(CS.y,1:NC,4);

%CS.C.Ma=single(CS.R.Ti.la<50);
%CS.R=[];
%CS.x=[];

if isempty(CS.C);CS.A.Se=(normm(CS.E.Sf,[],4)+CS.R.Ti.la).^(-1);%Precondition
else CS.A.Se=(normm(CS.E.Sf,[],4)+1e-9).^(-1);
end

%CS.tol=1e-2;

[x,lambda]=CSsolver(CS.y,CS.E,CS.EH,CS.A,CS.C,CS.R,CS.x,CS.nIt,CS.tolType,CS.tol,2,CS.nOu);
return




tic
if ~exist('rec','var')
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/recFull.mat');
end
toc

%rec.Alg.parXT.exploreMemory=1;%To explore convergence without running the main methods
%rec.Alg.WriteSnapshots=0;

rec=solveXT(rec);

