% function [sigmahat,dhat,tauhat,speed,sigmahat2,dhat2,lambda, ...
%    lambdahat,exitflag,numiter,x0]=QuESTimate(Y,k,optspeed)
%
% Obtain nonlinear shrinkage estimator of the covariance matrix
% and an estimator of population eigenvalues
%
% Input: Y           : data set
%                          each row    represents an observation
%                          each column represents a  variable
%        k (optional): if omitted or NaN (default) then QuESTimate demeans Y
%                      if present then QuESTimate does not demean Y
%                      k=0: dataset Y has mean zero
%                      k=1: dataset Y has already been demeaned
%                      k>1: Y has k classes, and each one has been demeaned 
%        optspeed (optional): controls optimizer speed and precision
%                                 MajorOptimalityTolerance = optspeed(1)
%                                 MajorIterationsLimit     = optspeed(2)
%                             if NaN then default optimizer settings are used
%                             if omitted (default) then optspeed=[1e-5 50]
%
% Outputs: sigmahat : nonlinear shrinkage estimator of covariance matrix
%                     according to equation (3.10) of Reference (A)
%          dhat     : eigenvalues of sigmahat
%          tauhat   : column vector containing an estimate of population
%                     covariance matrix eigenvalues
%          speed    : 4*1 vector containing speed in seconds
%                        1: optimization time measured by Matlab
%                        2: CPU time measured by Tomlab
%                        3: real time measured by Tomlab
%                        4: time required to run the QuEST function
%          sigmahat2: nonlinear shrinkage estimator of covariance matrix
%                     according to Theorem (6.2) of Reference (B)
%          dhat2    : eigenvalues of sigmahat2
%          lambda   : observed sample eigenvalues
%          lambdahat: fitted sample eigenvalues
%          exitflag : exit flag returned by optimizer
%          numiter  : number of iterations inside the optimizer
%          x0       : linear shrinkage starting point
%
% References: 
% (A) "Spectrum Estimation: A Unified Framework for Covariance Matrix 
% Estimation and PCA in Large Dimensions" by Olivier Ledoit and Michael 
% Wolf, Journal of Multivariate Analysis (July 2015), volume 139, pages 360-384
% 
% (B) "Optimal Estimation of a Large-Dimensional Covariance Matrix under 
% Stein’s Loss"  by Olivier Ledoit and Michael Wolf, forthcoming in
% Bernoulli
%
% dependencies: functions QuEST and QuESTgrad
%               TOMLAB/SNOPT

function [sigmahat,dhat,tauhat,speed,sigmahat2,dhat2,lambda, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(Y,k,optspeed)

% initialize
[n,p]=size(Y);
if (nargin<2)||isnan(k)
   Y=Y-repmat(mean(Y),[n 1]);
   k=1;
end

% set precision of convergence criterion in the optimizer 
if nargin<3
   optspeed=[1e-5 50];
end
MajorOptimalityTolerance = optspeed(1); 
MajorIterationsLimit     = optspeed(2);

% set SNOPT scaling parameters (probably these should not be modified)
os0=100;     % adaptive scaling factor for SNOPT objective function
ts=100;      % scaling factor for optimization variables  

% extract sample covariance matrix eigenvalues
sample=(Y'*Y)./(n-k);
sample=(sample+sample')./2;
[u,lambdamat]=eig(sample);
[lambda,jsort]=sort(diag(lambdamat));
u=u(:,jsort);
lambda(lambda<0)=0;
lambda(1:p-n+k)=0;

% compute linear shrinkage coefficient
Z=Y.^2;
phi=sum(sum(Z'*Z/(n-k)-2*(Y'*Y).*sample/(n-k)+sample.^2));
gama=norm(sample-mean(diag(sample))*eye(p),'fro')^2;
shrinkage=max(0,min(1,phi/gama/(n-k)));

% select starting point by linearly shrinking sample eigenvalues
mlambda=mean(lambda);
tauscaling=ts/mlambda;
x0=mlambda+sqrt(1-shrinkage).*(lambda-mlambda);
x_0=x0.*tauscaling;

% specify linear constraints
A=ones(1,p)./p;                     % trace must be preserved
b_L=mlambda*tauscaling;
b_U=b_L;

% impose bounds on the variables
x_L=lambda(1).*ones(p,1).*tauscaling;
x_U=lambda(p).*ones(p,1).*tauscaling;

% set up optimization problem
Prob=conAssign('QuESTmse','QuESTdmse',[],[],x_L,x_U,[],x_0,[],0,A,b_L,b_U);
Prob.SOL.optPar(48)=max(500,p+1); % superbasics limit
if ~isnan(MajorOptimalityTolerance)
   Prob.SOL.optPar(10)=MajorOptimalityTolerance; 
end
if ~isnan(MajorIterationsLimit)
   Prob.SOL.optPar(35)=MajorIterationsLimit; 
end

% pass user-defined parameters as global variables
clear global USER_DEFINED_PARAMETERS
global USER_DEFINED_PARAMETERS
USER_DEFINED_PARAMETERS.lambdatarget=lambda;
USER_DEFINED_PARAMETERS.n=n-k;
USER_DEFINED_PARAMETERS.objscaling=os0/mean((lambda-QuEST(x_0./tauscaling,n-k)).^2);
USER_DEFINED_PARAMETERS.tauscaling=tauscaling;

% call optimizer
starttime1=now;
Result=tomRun('SNOPT',Prob);
endtime1=now;
speed(1,1)=(endtime1-starttime1)*24*60*60;
speed(2,1)=Result.CPUtime;
speed(3,1)=Result.REALtime;
exitflag=Result.ExitFlag;
numiter=Result.Iter;
tauhat=sort(Result.x_k)./tauscaling;

% find optimal nonlinear shrinkage formula
starttime2=now;
[lambdahat,dhat,dhat2]=QuEST(tauhat,n-k);
endtime2=now;
speed(4,1)=(endtime2-starttime2)*24*60*60;

% compute nonlinear shrinkage estimator
sigmahat=u*diag(dhat)*u';
sigmahat2=u*diag(dhat2)*u';
