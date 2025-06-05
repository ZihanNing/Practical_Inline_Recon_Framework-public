function [x,En,gF,cF,rF]=IRWLSsolver(y,E,EH,no,di,de,A,C,R,nIt,tolType,tol,estGF,estCF,sepnorm,deb)

%IRWLSSOLVER   Performs an iteratively reweighted least squares approach to
%solve fitting problems in an l-p norm along a given set of dimensions of
%the data
%   [X,EN,GF,CF]=IRWLSSOLVER(Y,E,EH,{NO},{DI},{A},{C},{R},{C},{X},{GF},{NIT},{TOLTYPE},{TOL},{ESTGF},{ESTCF},{SEPNORM},{DEB})
%   * Y is the array of measures
%   * E is an encoding structure
%   * EH is a decoding structure
%   * {NO} is the norm of fitting to be optimized
%   * {DI} are the dimensions over which the previous norm applies (for 
%   other dimensions the norm is supposed to be 2)
%   * {DE} is the regularization for weight computation (its role comes 
%   from a Huber loss)
%   * {A} is a preconditioner structure
%   * {C} is a constrain structure
%   * {R} is a regularizer structure
%   * {NIT} is the maximum number of iterations
%   * {TOLTYPE} is the type of tolerance used for convergence
%   * {TOL} is the tolerance set for convergence
%   * {ESTGF} indicates whether to estimate the g-factor
%   * {ESTCF} indicates whether to estimate the chi2-factor
%   * {SEPNORM} introduces a normalization separable along the dimensions
%   not included in the norm
%   * {DEB} indicates whether to print information about convergence
%   ** X is the reconstruction
%   ** EN are the energies of successive reconstructions
%   ** GF is a g-factor value
%   ** CHI2 is the chi2 factor
%   ** RF are the residuals in space
%
%The NormwiseBackwardError stopping condition is based on Petr Tichy. 
%ON ERROR ESTIMATION IN THE CONJUGATE GRADIENT METHOD: NORMWISE BACKWARD 
%ERROR, Proceedings of ALGORITMY, pp. 323-332 2016.

if nargin<4 || isempty(no);no=2;end
if nargin<5 || isempty(di);di=1:ndims(y);end
if nargin<6 || isempty(de);de=1e-5;end%1e-6 is unstable
if nargin<7;A=[];end
if nargin<8;C=[];end
if nargin<9;R=[];end
if nargin<10 || isempty(nIt);nIt=300;end
if nargin<11 || isempty(tolType);tolType='Energy';end%Only one implemented now
if nargin<12 || isempty(tol);tol=1e-2;end%1e-5 would be consistent with other pieces of the code, but it takes too long!
if nargin<13 || isempty(estGF);estGF=[0 0];end
if nargin<14 || isempty(estCF);estCF=0;end
if nargin<15 || isempty(sepnorm);sepnorm=0;end
if nargin<16 || isempty(deb);deb=1;end
cF=[];gF=[];rF=[];

gpu=isa(y,'gpuArray');
NitTest=2;
tolTyp=stopCondition(tolType,'IRWLS');

EH.Mb=single(ones(size(y)));
if gpu;EH.Mb=gpuArray(EH.Mb);end

non=2;%Initial homotopy value
redN=0.7;%Reduction factor to approach the target norm, it was 0.8 before!

err=inf;
En=inf(2,nIt);

nn=zeros(1,nIt+1);
[x,~,~,~,nn(1)]=CGsolver(y,E,EH,A,C,R,[],nIt,tolType,tol,0,0,deb);

if no<=-1;nIt=3;end

for n=1:nIt
    if no==2;break;end
    if no>-1;non=max(no,redN*non);else non=no;end
    EH.Mb(:)=1;
    [En(:,n),EH.Mb]=computeEnergy(y,x,E,R,EH,no,di,non,de,[],sepnorm);

    %if isinf(En(1,n))
    %    fprintf('THE ENERGY IS NOT FINITE. IRWLS IS ABORTED\n');break;
    %end
    if mod(n,NitTest)==1 && tolTyp==6
        if deb>=2;fprintf('IRWLS It %d - Ene Fid: %0.2g / Ene Reg: %0.2g / Ene Tot: %0.2g\n',n,En(1,n),En(2,n),sum(En(:,n)));end
    end    
    if tolTyp==6 && mod(n,NitTest)==1 && n>1;err=sum(En(:,n-NitTest)-En(:,n))/sum(En(:,n-NitTest));%err=sqrt(max(0,err));
    end

    if deb>=2 && (tolTyp<6 || mod(n,NitTest)==1);fprintf('It %d - Err %s: %0.2g\n',n,tolType,err);end
    if err<tol;break;end    
    if no<=-1;nItX=300;else nItX=nIt;end
    [x,~,~,~,nn(n+1)]=CGsolver(y,E,EH,A,C,R,x,nItX,tolType,tol,0,0,deb);
    if no>-1 && (non==no || non<0.1) && isinf(En(1,n));break;end
end

if n==nIt && no>-1;fprintf('IRWLS solver terminated without reaching convergence\n');end
En(:,(1:nIt)>n)=[];

if estGF(1)>0
    if deb>=2;fprintf('Number of iterations for g-factor calculations: %d\n',sum(nn));end
    gF=solveG(E,EH,A,C,R,sum(nn),tolType,0,deb,estGF,y,x);
end
if estCF;[cF,rF]=solveC(x,y,E,EH);end
