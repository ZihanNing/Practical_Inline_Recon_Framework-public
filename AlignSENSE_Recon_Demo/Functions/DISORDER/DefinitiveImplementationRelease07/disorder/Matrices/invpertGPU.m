function [Xp,fl]=invpertGPU(X,Xp,D,fl,nit,typ,th,tol,nmax)

%INVPERTGPU   Computes the inverse of a matrix for which the inverse of a
%another matrix obtained by a small perturbation of the matrix of interest
%is already known. It uses the formula (X+D)^{-1}=(I+X^{-1}D)X^{-1} which
%is a simplification of the formula in (3) of [1] FC Chang, "Inversion of 
%a perturbed matrix," Appl Math Let 19:169-173, 2006. The key idea is that 
%the new inverse may have better properties depending on the perturbation.
%In our case we promote diagonal dominance by factoring out X so that an
%inverse update is performed by using the Neumann series. Other usages
%could be for low rank or sparse perturbations
%   XP=INVPERTGPU(X,{D},{XP},{FL},{TYP},{TOL},{NMAX})
%   * X are the matrices plus the perturbations, to compute standard 
%   inverses
%   * {D} are the perturbations of a given set of matrices. If empty or not
%   existing the function computes the standard inverse on the input
%   * {XP} are the inverses of those matrices. If empty or not existing,
%   the function computes the standard inverse on the input
%   * {FL} are flags that indicate the number of iterations for which an
%   external iterative algorithm has kept on computing iterative inverses.
%   If set to -1 then standard inverse is computed in this step. Useful to
%   set them to minus one from time to time to prevent divergence of
%   iterative algorithms
%   * {NIT} are the number of iterations after which to refresh the inverse
%   by computing the full inverses
%   * {TYP} indicates the algorithm to compute the inverses. Currently 
%   only 'ddom' is supported
%   * {TH} is a certain threshold to activate accelerated inversion. It
%   defaults to 0.5, which was judged to be adequate for diagonal
%   dominances
%   * {TOL} is a certain tolerance for inversion convergence. Defaults to
%   0
%   * {NMAX} are the maximum number of iterations for inversion
%   convergence. Defaults to 1
%   * XP are the inverses of the perturbed matrices
%   * FL returns an update on the number of previous iterations for which
%   an iterative inverse was been computed
%

N=size(X);N(end+1:3)=1;
if nargin<2;D=[];end
if nargin<3;Xp=[];end
if nargin<4 || isempty(fl);fl=single(zeros([1 1 N(3)]));end
if nargin<5 || isempty(typ);typ='ddom';end
if nargin<6 || isempty(nit);nit=5;end
if nargin<7 || isempty(th);th=0.5;end
if nargin<8 || isempty(tol);tol=0;end
if nargin<9 || isempty(nmax);nmax=1;end

if isempty(Xp) || isempty(D)
    Xp=matfun(@inv,X);
    fl(:)=0;
else
    I=eye(N(1),'like',D);
    D=bsxfun(@plus,I,matfun(@mtimes,Xp,D));
    %RATIO OF DIAGONAL DOMINANCES    
    vR=gather(ddomGPU(D));
    fl(vR<th)=-1;
    stInv=fl==-1;
    if any(stInv(:));Xp(:,:,stInv)=matfun(@inv,X(:,:,stInv));end
    if any(~stInv(:))
        [Xpi,Di]=parUnaFun({Xp,D},@dynInd,~stInv,3);
        if strcmp(typ,'ddom')
            Di=invddomGPU(Di,tol,nmax);%Only a single iteration currently although more would be possible
        else
            error('Inversion method %d is not contemplated',typ);
        end       
        Xpi=matfun(@mtimes,Di,Xpi);
        Xp=dynInd(Xp,~stInv,3,Xpi);
    end    
    fl=fl+1;
    fl(fl==nit)=-1;
end
