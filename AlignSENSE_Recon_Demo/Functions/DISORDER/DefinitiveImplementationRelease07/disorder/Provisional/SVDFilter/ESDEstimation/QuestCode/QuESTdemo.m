% QuESTdemo.m - demonstrate usage of QuESTimate.m function
%               through a simulation
%
% Reference: "Spectrum Estimation: A Unified Framework
% for Covariance Matrix Estimation and PCA in Large Dimensions"
% by Olivier Ledoit and Michael Wolf (2013), Section 5.1.1, Table 2
%
% dependencies: functions QuEST, QuESTgrad, QuESTmse, QuESTdmse, QuESTimate,
%               and fmincon from the MATLAB Optimization Toolbox
%               (or SNOPT/TOMLAB third-party nonlinear optimizer)

clear

% set parameters
p=1000;
n=3000;
randn('state',0)

% specify population eigenvalues
z=linspace(0,1,p)';
z1=z(z<=1/2);
expo=3;
tau1=(1-(1-(2.*z1).^expo).^(1/expo))./2;
tau=1+9.*[tau1;1-flipud(tau1)];
sigma=diag(tau);

% simulate data set and estimate population eigenvalues
Y=randn(n,p)*sqrtm(sigma);
[sigmahat,dhat,tauhat,speed,sigmahat2,dhat2,lambda, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(Y,0);
% change here if you want to use TOMLAB/SNOPT instead of Matlab optimizer
%[sigmahat,dhat,tauhat,speed,sigmahat2,dhat2,lambda, ...
%   lambdahat,exitflag,numiter,x0]=QuESTimates(Y,0);

% compute finite-sample optimal rotation-equivariant estimator
Y=Y-repmat(mean(Y),[n 1]);
sample=(Y'*Y)./n;
[u,lambdamat]=eig(sample);
[lambda,jsort]=sort(diag(lambdamat));
u=u(:,jsort);
dstar=diag(u'*sigma*u);
dstar2=1./diag(u'*diag(1./tau)*u);

% plot eigenvalues results
figure(1)
if median(get(gcf,'color'))<0.5
   whitebg
end
plot(1:p,tau,'.b',1:p,tauhat,'.r')
set(gcf,'position',[100 150 700 500])
xl=xlabel('Index of Eigenvalues');
yl=ylabel('Population Eigenvalues');
ti=title('Accuracy of Eigenvalues Estimator');
leg=legend('True','Estimated','Location','SouthEast');
set(gca,'fontsize',14)
set(xl,'fontsize',14)
set(yl,'fontsize',14)
set(ti,'fontsize',14)
set(leg,'fontsize',10)

% plot nonlinear shrinkage results
figure(2)
plot(1:p,dstar,'.b',1:p,dhat,'-r')
set(gcf,'position',[200 100 700 500])
xl=xlabel('Index of Eigenvalues');
yl=ylabel('Shrunk Eigenvalues');
ti=title('Accuracy of Nonlinear Shrinkage Estimator');
leg=legend('True','Estimated','Location','SouthEast');
set(gca,'fontsize',14)
set(xl,'fontsize',14)
set(yl,'fontsize',14)
set(ti,'fontsize',14)
set(leg,'fontsize',10)
kids=get(gca,'children');
set(kids(1),'linewidth',2)

% plot nonlinear shrinkage results
figure(3)
plot(1:p,dstar2,'.b',1:p,dhat2,'-r')
set(gcf,'position',[300 50 700 500])
xl=xlabel('Index of Eigenvalues');
yl=ylabel('Shrunk Eigenvalues (Alternative Formula)');
ti=title('Accuracy of Alternative Nonlinear Shrinkage Estimator');
leg=legend('True','Estimated','Location','SouthEast');
set(gca,'fontsize',14)
set(xl,'fontsize',14)
set(yl,'fontsize',14)
set(ti,'fontsize',14)
set(leg,'fontsize',10)
kids=get(gca,'children');
set(kids(1),'linewidth',2)

