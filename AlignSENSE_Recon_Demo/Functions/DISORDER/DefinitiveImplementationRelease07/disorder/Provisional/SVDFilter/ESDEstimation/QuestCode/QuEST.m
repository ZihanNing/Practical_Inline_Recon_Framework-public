% function [lambda,d,d2]=QuEST(tau,n)
%
% Quantized Eigenvalues Sampling Transform
% 
% See article "Spectrum Estimation: A Unified Framework for Covariance 
% Matrix Estimation and PCA in Large Dimensions" by Olivier Ledoit and 
% Michael Wolf, Journal of Multivariate Analysis (July 2015), volume 139,
% pages 260-384, Section 2.2, equations (2.11)-(2.17)
%
% Inputs : tau   : column vector containing population eigenvalues
%          n     : positive scalar integer containing sample size
%
% Outputs: lambda: sample eigenvalues
%          d     : nonlinear shrinkage according to equation (3.10)
%          d2    : nonlinear shrinkage according to Theorem 6.2 on page 19
% of "Optimal Estimation of a Large-Dimensional Covariance Matrix under 
% Stein’s Loss"  by Olivier Ledoit and Michael Wolf, University of Zurich, 
% Department of Economics, Working Paper No. 122, revised version Nov 2014
%
% To see help message, type: QuEST('help')
%     
% Copyright Olivier Ledoit and Michael Wolf (2013-2016)

function [lambda,d,d2]=QuEST(tau,n)

if nargin==1
   if isequal(tau,'help')
      disp('  function [lambda,d,d2]=QuEST(tau,n)')
      disp(' ')
      disp('  Quantized Eigenvalues Sampling Transform')
      disp(' ')
      disp('  See article "Spectrum Estimation: A Unified Framework for Covariance')
      disp('  Matrix Estimation and PCA in Large Dimensions" by Olivier Ledoit and')
      disp('  Michael Wolf, Journal of Multivariate Analysis (July 2015), volume 139,')
      disp('  pages 260-384, Section 2.2, equations (2.11)-(2.17)')
      disp(' ')
      disp('  Inputs : tau   : column vector containing population eigenvalues')
      disp('           n     : positive scalar integer containing sample size')
      disp(' ')
      disp('  Outputs: lambda: sample eigenvalues')
      disp('           d     : nonlinear shrinkage')
      disp('           d2    : nonlinear shrinkage according to Theorem 6.2 on page 19')
      disp('  of "Optimal Estimation of a Large-Dimensional Covariance Matrix under')
      disp('  Stein''s Loss"  by Olivier Ledoit and Michael Wolf, University of Zurich,') 
      disp('  Department of Economics, Working Paper No. 122, revised version Nov 2014')
      disp(' ')
      disp('  To see help message, type: QuEST(''help'')')
      disp(' ')
      disp('  Copyright Olivier Ledoit and Michael Wolf (2013-2016)')
      disp(' ')
   end
else
   % define inputs
   def.tau=tau;
   def.n=n;
   def=checkinputs01(def);

   % call the subfunctions
   [int,dis,den,sol,net,sup,def]=intfun03(def);

   % post-process the output
   lambda=int.lambda;
   d=int.delta;
   d2=int.deltainv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function def=checkinputs01(def)
%
% check that the function inputs obey our restrictions
%
% required inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%                if not in weakly ascending order, will be sorted
%   def.n      : number of observations
%
% optional fields that will be inserted if missing:
%   def.p      : length of def.tau (optional, will be inserted if missing)
%   def.c      : p/n (optional, will be inserted if missing)
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% dependency : none
%
% Version 01: fixed equal weights, variable population eigenvalues

function def=checkinputs01(def)

% check structure of input argument
if ~isstruct(def)
   error('input argument must be a structure')
end
if ~all(isfield(def,{'tau','n'}))
   error('input argument must have fields tau and n at least')
end
if (~isnumeric(def.tau))|(~isnumeric(def.n))
   error('fields tau and n must be numerical')
end
if (~isreal(def.tau))|(~isreal(def.n))
   error('fields tau and n must be real')
end

% check dimensions
if ~isvector(def.tau)
   error('tau must be vector')
end
if size(def.tau,2)>1
   error('tau must be column vector')
end
if ~isscalar(def.n)
   error('n must be scalar')
end
if def.n<=0
   error('n must be strictly positive')
end

% insert fields p and c if missing
if ~isfield(def,'p')
   def.p=size(def.tau,1);
end
if ~isfield(def,'c')
   def.c=def.p/def.n;
end

% check p and c
if (~isnumeric(def.p))|(~isnumeric(def.c))
   error('fields p and c must be numerical')
end
if (~isscalar(def.p))|(~isscalar(def.c))
   error('fields p and c must be scalar')
end
if (~isreal(def.p))|(~isreal(def.c))
   error('fields p and c must be scalar')
end
if def.p~=round(def.p)
   error('p must be integer')
end
if (def.p<=0)|(def.c<=0)
   error('p and c must be strictly positive')
end
if def.p~=size(def.tau,1)
   error('p must be length of tau')
end
if def.c~=def.p/def.n
   error('c must be equal to p/n')
end

% check n is an integer if it is strictly less than p
if def.n<def.p
   if def.n~=round(def.n)
      error('n must be integer if it is strictly less than p')
   end
end

% check tau
lefttoler=1e-8;
toler=1e-12;
def.tau((def.tau<0)&(abs(def.tau)<lefttoler))=0;
def.tau((def.tau>=0)&(abs(def.tau)<toler))=0;
if any(def.tau<0)
   disp(['min(tau) = ' num2str(min(def.tau))])
   error('the elements of tau must be non-negative')
end
if def.p>1
   if any(diff(def.tau,1,1)<0)
      def.tau=sort(def.tau,1,'ascend');
   end
end
for i=2:def.p
   if def.tau(i)<def.tau(i-1)+toler;
      def.tau(i)=def.tau(i-1);
   end
end

% insert field weight if missing
if ~isfield(def,'weight')
   def.weight=ones(def.p,1)./def.p;
end
if ~isnumeric(def.weight)
   error('field weight must be numerical')
end
if ~isreal(def.weight)
   error('weight must be real')
end
if ~isvector(def.weight)
   error('weight must be vector')
end
if ~isequal(size(def.weight),size(def.tau))
   error('weight and tau must have the same sizes')
end
if ~isequal(def.weight,ones(def.p,1)./def.p)
   error('weight must be equal to ones(p,1)./p')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [int,dis,den,sol,net,sup,def]=intfun03(def)
%
% interpolate inverse cdf of sample eigenvalues F^{-1}(x) onto linear grid
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   int.lambda : p*1 eigenvalues representing sample spectral distribution
%   int.quant  : numint*1 cell array mapping support intervals to quantiles
%   int.nquant : numint*1 vector containing the length of int.quant{i}
%   int.bin    : numint*1 cell array mapping quantiles to F(x)
%   int.indic : numint*1 cell array of matrices indicating how to integrate
%   int.delta : optimal nonlinear shrinkage for cov matrix eigenvalues
%   int.deltainv : optimal nonlinear shrinkage for inverse cov matrix
%   int.limit0 : if p=n then the limit of delta(x) as x->0+
%   int.limitinv0 : if p=n then the limit of deltainv(x) as x->0+
%
%   see also function supfun04.m for outputs in the sup structure,
%            function netfun02.m for outputs in the net structure,
%            function solfun01.m for outputs in the sol structure,
%            function denfun02.m for outputs in the den structure,
%        and function disfun03.m for outputs in the dis structure
%
% direct dependency : disfun04
%
% Version 03 : use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

% supporting calculations: Olivier's notepad 12 May and 9 June 2012
%                          and 6 March 2013

function [int,dis,den,sol,net,sup,def]=intfun03(def)

global US_D

% get limiting spectral cdf of sample eigenvalues
[dis,den,sol,net,sup,def]=disfun04(def);

% initialize
int.lambda=NaN+zeros(def.p,1);
int.delta=NaN+zeros(def.p,1);
int.deltainv=NaN+zeros(def.p,1);
int.nquant=NaN+zeros(sup.numint,1);
int.bin=cell(sup.numint,1);
int.quant=cell(sup.numint,1);
int.indic=cell(sup.numint,1);

% take care of null sample eigenvalues in p>=n case
if max(sum(def.tau==0),def.p-def.n)>0
   int.lambda(1:max(sum(def.tau==0),def.p-def.n))=0;
   if def.p-def.n>sum(def.tau==0)
      lobound=min(def.tau)-def.c*sum(def.weight.*def.tau)-1;
      hibound=(1-sup.def2.weight(1))*sup.def2.tau(1);
      uFbar0=fzero(@(u) sum(sup.def2.weight.*sup.def2.tau ...
         ./(sup.def2.tau-u))-1./sup.def2.c,[lobound hibound]);
      nulladjust=uFbar0/(1-def.c);
      nulladjustinv=1/(1/(1-1./def.c) ...
         *sum(sup.def2.weight./sup.def2.tau)+1/uFbar0);
      int.delta(1:def.p-def.n)=nulladjust;
      int.deltainv(1:def.p-def.n)=nulladjustinv;
   else
      int.delta(1:sum(def.tau==0))=0;
      int.deltainv(1:sum(def.tau==0))=0;
   end      
end
if def.p==def.n
   if all(def.tau>0)
      int.limit0=1/sum(sup.def2.weight./sup.def2.tau);
      int.limitinv0=sum(sup.def2.weight./sup.def2.tau) ...
         /sum(sup.def2.weight./sup.def2.tau.^2);
   else
      int.limit0=0;
      int.limitinv0=0;
   end
end

% loop over support intervals
for i=1:sup.numint
   
   % select the points of the distribution that belong to this interval
   idx=find(dis.interval==i);
   F=dis.G(idx);
   [F,i1]=unique(F);
   idx=idx(i1);
   x=dis.zeta(idx).^den.a;
   nidx=length(idx);
   
   % define the quantiles that belong to this interval
   int.nquant(i)=sup.numeig(i)+1;
   int.quant{i}=linspace(F(1),F(nidx),int.nquant(i))';
   
   % integrate x over the interval [F(k),F(k+1)]
   intxdF=(F(2:nidx)-F(1:nidx-1)).*(x(1:nidx-1)+x(2:nidx))./2;

   % find out what bin each quantile belongs to
   edges=[F(1:nidx-1);F(nidx)+1];
   [discard,int.bin{i}]=histc(int.quant{i},edges);
   if (int.bin{i}(1)~=1)|(int.bin{i}(int.nquant(i))~=nidx-1)
      error('unexpected bin allocation')
   end

   % integrate x over the interval [F(k),quant(i)] 
   intxdF2=(int.quant{i}-F(int.bin{i})) ...
      .*(x(int.bin{i})+(int.quant{i}-F(int.bin{i})) ...
      .*(x(int.bin{i}+1)-x(int.bin{i})) ...
      ./(F(int.bin{i}+1)-F(int.bin{i}))./2);
   
   % compute integral of x over the interval [F(1),quant(i)]
   int.indic{i}=(repmat(2:nidx,[int.nquant(i) 1]) ...
      <=repmat(int.bin{i},[1 nidx-1]));
   intxdF3=sum(repmat(intxdF',[int.nquant(i) 1]).*int.indic{i},2)+intxdF2;
   
   % compute sample eigenvalues
   int.lambda(round(F(1)*def.p+1):round(F(nidx)*def.p))=diff(intxdF3,1,1).*def.p;
   
   % do same for nonlinear shrinkage of covariance matrix eigenvalues
   y=NaN+zeros(size(x));
   denominator=abs(1-def.c.*dis.m_LF(idx)).^2;
   if def.p==def.n
      j=((x==0)&(denominator==0));
      y(j)=int.limit0;
      y(~j)=x(~j)./denominator(~j);
   else
      y=x./denominator;
   end
   intydF=(F(2:nidx)-F(1:nidx-1)).*(y(1:nidx-1)+y(2:nidx))./2;
   intydF2=(int.quant{i}-F(int.bin{i})) ...
      .*(y(int.bin{i})+(int.quant{i}-F(int.bin{i})) ...
      .*(y(int.bin{i}+1)-y(int.bin{i})) ...
      ./(F(int.bin{i}+1)-F(int.bin{i}))./2);
   intydF3=sum(repmat(intydF',[int.nquant(i) 1]).*int.indic{i},2)+intydF2;
   int.delta(round(F(1)*def.p+1):round(F(nidx)*def.p))= ...
      diff(intydF3,1,1).*def.p;
   if def.p==def.n
      int.delta(int.delta<int.limit0)=int.limit0;
   end

   % do same for nonlinear shrinkage of inverse covariance matrix
   z=NaN+zeros(size(x));
   denominatorinv=1+def.c-2.*def.c.*real(dis.m_LF(idx));
   if def.p==def.n
      jinv=((x==0)&(denominatorinv==0));
      z(jinv)=int.limitinv0;
      z(~jinv)=x(~j)./denominatorinv(~j);
   else
      z=x./denominatorinv;
   end
   intzdF=(F(2:nidx)-F(1:nidx-1)).*(z(1:nidx-1)+z(2:nidx))./2;
   intzdF2=(int.quant{i}-F(int.bin{i})) ...
      .*(z(int.bin{i})+(int.quant{i}-F(int.bin{i})) ...
      .*(z(int.bin{i}+1)-z(int.bin{i})) ...
      ./(F(int.bin{i}+1)-F(int.bin{i}))./2);
   intzdF3=sum(repmat(intzdF',[int.nquant(i) 1]).*int.indic{i},2)+intzdF2;
   int.deltainv(round(F(1)*def.p+1):round(F(nidx)*def.p))=diff(intzdF3,1,1).*def.p;
   if def.p==def.n
      int.deltainv(int.deltainv<int.limitinv0)=int.limitinv0;
   end
end

% make results available for the gradient function
US_D.int=int;
US_D.dis=dis;
US_D.den=den;
US_D.sol=sol;
US_D.net=net;
US_D.sup=sup;
US_D.def=def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dis,den,sol,net,sup,def]=disfun04(def)
%
% compute F(x) the cdf of sample eigenvalues
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   dis.x        : abscissa of limiting distribution of sample eigenvalues
%   dis.zeta     : dis.x.^(1/a)
%   dis.G        : G(zeta)=F(zeta^a) 
%   dis.numx     : common length of column vectors dis.F and dis.x
%   dis.interval : index of support interval to which dis.x belongs
%   dis.weight0  : mass point of sample cdf at zero
%   dis.m_LF     : modified Stieltjes transform of sample eigenvalues
%   dis.fromweight : sup.numint*1 vector containing value of F at the
%                    beginning of support interval 
%   dis.toweight : same for the end of the support interval
%   dis.Graw     : numint*1 cell array containing raw integral of g 
%
%   see also function supfun04.m for outputs in the sup structure,
%            function netfun02.m for outputs in the net structure,
%            function solfun01.m for outputs in the net structure,
%        and function denfun02.m for outputs in the den structure
%
% direct dependency : denfun03
%
% Version 04 : use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

function [dis,den,sol,net,sup,def]=disfun04(def)

global US_D

% get limiting spectral density of sample eigenvalues
[den,sol,net,sup,def]=denfun03(def);

% initialize variables
dis.numx=net.nxi+2*sup.numint;
dis.x=NaN+zeros(dis.numx,1);
dis.zeta=NaN+zeros(dis.numx,1);
dis.G=NaN+zeros(dis.numx,1);
dis.m_LF=NaN+zeros(dis.numx,1);
dis.interval=NaN+zeros(dis.numx,1);
dis.weight0=max(sum(def.weight(def.tau==0)),1-1/def.c);
dis.fromweight=NaN+zeros(sup.numint,1);
dis.toweight=NaN+zeros(sup.numint,1);
dis.Graw=cell(sup.numint,1);

% loop over support intervals
for i=1:sup.numint
   from=sum(sup.numeig(1:i-1)*net.mult+2)+1;
   intlength=sup.numeig(i)*net.mult+2;
   to=from+intlength-1;
   dis.interval(from:to)=i;
   isin=(net.interval==i);
   dis.x(from:to)=[sup.endpoint(i,1);den.x(isin);sup.endpoint(i,2)];
   dis.zeta(from:to)=[sup.endpoint(i,1)^(1/den.a);den.zeta(isin); ...
      sup.endpoint(i,2)^(1/den.a)];
   dis.m_LF(from:to)=[sup.m_LF(i,1);den.m_LF(isin);sup.m_LF(i,2)];
   g=[0;den.g(isin);0];
   dis.fromweight(i)=dis.weight0+sum(sup.numeig(1:i-1))/def.p;
   dis.toweight(i)=dis.weight0+sum(sup.numeig(1:i))/def.p;
   dis.Graw{i}=[0;cumsum(diff(dis.zeta(from:to),1,1) ...
      .*(g(1:intlength-1,:)+g(2:intlength,:))./2,1)];
   dis.G(from:to)=dis.fromweight(i) ...
      +(dis.toweight(i)-dis.fromweight(i)) ...
      .*dis.Graw{i}./dis.Graw{i}(intlength);
end   

% make results available for the gradient function
US_D.dis=dis;
US_D.den=den;
US_D.sol=sol;
US_D.net=net;
US_D.sup=sup;
US_D.def=def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [den,sol,net,sup,def]=denfun03(def)
%
% compute x and f(x) the density of sample eigenvalues
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   den.x    : abscissa of the plot of the density of sample eigenvalues
%   den.f    : density of sample eigenvalues evaluated at den.x
%   den.a    : exponent of the power transform of the abscissa
%   den.zeta : den.x.^(1/den.a)
%   den.g    : derivative of G(zeta)=F(zeta^den.a) wrt zeta
%   den.m_LF : modified Stieltjes transform evaluated at den.x
%   see also function supfun04.m for outputs in the sup structure,
%            function netfun01.m for outputs in the net structure
%        and function solfun02.m for outputs in the net structure
%
% direct dependency : solfun02
%
% Version 03 : use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

function [den,sol,net,sup,def]=denfun03(def)

global US_D

% set exponent of the power transform of the abscissa
den.a=4;

% get value of -1/m_Fbar(lambda) for lambda in the support of F
[sol,net,sup,def]=solfun02(def);

% compute Stieltjes transform of F evaluated at x
TAU=repmat(def.tau',[net.nxi 1]);
den.m_LF=sum(repmat(def.weight',[net.nxi 1]) ...
   .*TAU./(TAU-repmat(sol.zxi,[1 def.p])),2);

% compute abscissa
den.x=real(sol.zxi-def.c.*sol.zxi.*den.m_LF);

% compute ordinate
den.f=1./pi./def.c.*imag(sol.zxi)./abs(sol.zxi).^2;

% compute abscissa
den.zeta=den.x.^(1/den.a);

% compute ordinate
den.g=den.a.*den.zeta.^(den.a-1).*den.f;

% make results available for the gradient function
US_D.den=den;
US_D.sol=sol;
US_D.net=net;
US_D.sup=sup;
US_D.def=def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [sol,net,sup,def]=solfun02(def)
%
% compute imaginary part of the solution to the Marcenko-Pastur equation
% in -1/m_Fbar(z) space as per Section 2.3 of Ledoit and Wolf (2012)
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   sol.zxi    : (complex) value of -1/m_Fbar(z) when z is on the real line
%                in the support of the sample spectral distribution
%   see also function supfun04.m for outputs in the sup structure
%        and function netfun01.m for outputs in the net structure
%
% direct dependency : netfun02
%
% Version 02: integrate auxiliary function solopt01.m into solfun

% supporting calculations: Olivier's notepad 8 May 2012

function [sol,net,sup,def]=solfun02(def)

global US_D

% get grid of real parts of -1/m_Fbar(z)
[net,sup,def]=netfun02(def);

% initialize
yxi=zeros(net.nxi,1);
def.kernel2=def.c.*def.weight.*def.tau.^2;
lobound=zeros(net.nxi,1);
hibound=zeros(net.nxi,1);

% loop over points on the xi grid in u-space
for iz=1:net.nxi
   dif2=(def.tau-net.xi(iz)).^2;
   mindif2=min(dif2);
   ismin=find(dif2==mindif2);
   lobound(iz)=sqrt(max(0,def.c*sum(def.weight(ismin) ...
      .*def.tau(ismin).^2,1)-mindif2))/2;
   hibound(iz)=sqrt(sum(def.kernel2,1)-mindif2)+1;
   if lobound(iz)==hibound(iz)
      yxi(iz)=lobound;
   elseif lobound(iz)<hibound(iz)
      yxi(iz)=fzero(@(y) 1-sum(def.kernel2./((def.tau-net.xi(iz)).^2+y^2)), ...
         [lobound(iz) hibound(iz)]);
   else
      error('unexpected ordering of the bounds')
   end
end

% compute function output
sol.zxi=net.xi+sqrt(-1).*yxi;

% make results available to the gradient function
US_D.sol=sol;
US_D.net=net;
US_D.sup=sup;
US_D.def=def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [net,sup,def]=netfun02(def)
%
% compute grid over [u_1,u_2] which is the real part of -1/m_Fbar(z)
% as per Ledoit and Wolf (Annals of Statistics, 2012) Section2 2.2-2.3
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   net.xi       : evenly spaced grid over the image of the support
%                     of F or Fbar through the function -1/m_Fbar
%                     dimension nxi*q
%   net.sup2xi   : linear mapping from sup.u_Fbar(:) to net.xi
%   net.nxi      : length of net.xi
%   net.mult     : integer equal to net.xi/def.p
%   net.interval : index of support interval to which xi gridpoint belongs
%
%   see also function supfun07.m for outputs in the sup structure
%
% direct dependency : supfun07
%
% Version 02: tailor xi grid to spectrum separation

function [net,sup,def]=netfun02(def)

global US_D

% parameters
minnxi=100;
bpar=[0.5 0.5]; % more points near edges because square root behavior

% determine nxi as an integer multiple of the number of nonzero eigenvalues
nxi0=max(minnxi,min(sum(def.tau~=0),def.n));
net.mult=ceil(nxi0/min(sum(def.tau~=0),def.n));
net.nxi=net.mult*min(sum(def.tau~=0),def.n); 

% compute the support
[sup,def]=supfun07(def);

% construct linear map from support edges to grid in u-space
net.sup2xi=zeros(net.nxi,sup.numint*2);
net.interval=NaN+zeros(net.nxi,1);
for i=1:sup.numint
   grid01=(sin((pi/2).*linspace(0,1,sup.numeig(i)*net.mult+2)')).^2;
   grid01([1 sup.numeig(i)*net.mult+2])=[];
   from=sum(sup.numeig(1:i-1)*net.mult)+1;
   to=from+sup.numeig(i)*net.mult-1;
   net.sup2xi(from:to,i)=1-grid01;
   net.sup2xi(from:to,i+sup.numint)=grid01;
   net.interval(from:to)=i;
end

% construct xi grid   
net.xi=net.sup2xi*reshape(sup.u_Fbar,[US_D.sup.numint*2 1]);

% make results available to the gradient function
US_D.net=net;
US_D.sup=sup;
US_D.def=def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [sup,def]=supfun07(def)
%
% compute support in u=-1/m_Fbar(z) space
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   sup.endpoint: numint*2 matrix containing endpoints of the numint 
%          intervals constituting the support of sample eigenvalues density
%   sup.u_Fbar  : value of -1/m_Fbar(x) at the endpoints of the support
%   sup.m_LF    : value of    m_LF(x)   at the endpoints of the support
%   sup.numeig  : numint*1 vector containing the number of eigenvalues in 
%                  each interval of the support
%   sup.numint  : number of distinct intervals in the support of F
%   sup.def2    : alternative problem definition where all the 
%                 tau's are unique and strictly positive
%
% direct dependency : none
%
% Version 07: integrate auxiliary function supdif01.m into supfun

% Supporting calculations in Olivier's notepad around 6 May 2012

function [sup,def]=supfun07(def)

global US_D

% find list of unique population eigenvalues
t=unique(def.tau);
zeroweight=0;
if length(t)==1
   w=sum(def.tau==t)./def.p;
else
   w=hist(def.tau,t)'./def.p;
   t(w==0)=[];
   w(w==0)=[];
   if any(t==0)
      zeroweight=w(t==0);
      w(t==0)=[];
      t(t==0)=[];
   end
end
if abs(sum(w)-mean(def.tau~=0))>1e-5./def.p
   error('weights must sum up to mean(tau~=0)')
end
K=length(t);
sup.def2.weight=w;
sup.def2.tau=t;
sup.def2.c=def.c;
sup.def2.ntau=length(t);
sup.def2.WT2=sup.def2.weight'.*sup.def2.tau'.^2;

% get numerical tolerance
opt=optimset('fzero');

% look for spectral separation if it exists
separation=NaN+zeros(K-1,4);
if K>1
   
   % preliminary test to see where spectral separation might occur
   xhat=(t(1:K-1).*t(2:K)).^(2/3) ...
      .*(w(1:K-1).^(1/3).*t(2:K  ).^(1/3)+w(2:K).^(1/3).*t(1:K-1).^(1/3)) ...
      ./(w(1:K-1).^(1/3).*t(1:K-1).^(2/3)+w(2:K).^(1/3).*t(2:K  ).^(2/3));
   warnstate=warning;
   warning('off','MATLAB:divideByZero')
   theta=w(1:K-1).*t(1:K-1).^2./(t(1:K-1)-xhat).^2 ...
      +  w(2:K  ).*t(2:K  ).^2./(t(2:K  )-xhat).^2;
   warning(warnstate)
   theta(~isfinite(theta))=max(theta(~isfinite(theta)));
   if K>2
      warnstate=warning;
      warning('off','MATLAB:divideByZero')
      minleft =repmat((w.*t.^2)',[K-1 1]) ...
         ./(repmat(t',[K-1 1])-repmat(t(2:K),[1 K])).^2;
      minright=repmat((w.*t.^2)',[K-1 1]) ...
         ./(repmat(t',[K-1 1])-repmat(t(1:K-1),[1 K])).^2;
      warning(warnstate)
      minleft( logical(triu(ones(K-1,K),0)))=0;
      minright(logical(tril(ones(K-1,K),1)))=0;
      theta=theta+sum(minleft,2);
      theta=theta+sum(minright,2);
   end
   kk=find(theta<1/def.c);
   nk=length(kk);
   if nk>0
      
      % find out where the minimum of phi lies relative to xhat
      dpsi=2.*sum(repmat(sup.def2.WT2,[nk 1]) ...
         ./(repmat(sup.def2.tau',[nk 1]) ...
         -repmat(xhat(kk),[1 sup.def2.ntau])).^3,2);
      
      % loop over points where separation might occur
      for i=1:nk
         k=kk(i);
         
         % three cases to find out the minimum of phi over the interval
         if abs(dpsi(i))<=eps(xhat(k))
            xstar=xhat(k);
         elseif dpsi(i)>0
            lobound=t(k  )+(2*w(k  )*t(k  )^2/( 2*w(k+1) ...
               *t(k+1)^2/(t(k+1)-xhat(k))^3+dpsi(i)))^(1/3);
            if abs(lobound-xhat(k))<=opt.TolX
               xstar=(lobound+xhat(k))/2;
            else
               xstar=fzero(@(x) 2.*sum(sup.def2.WT2 ...
                  ./(sup.def2.tau'-repmat(x,[1 sup.def2.ntau])).^3,2), ...
                  [lobound xhat(k)]);
            end
         elseif dpsi(i)<0
            hibound=t(k+1)-(2*w(k+1)*t(k+1)^2/(-2*w(k  ) ...
               *t(k  )^2/(t(k  )-xhat(k))^3-dpsi(i)))^(1/3);
            if abs(hibound-xhat(k))<=opt.TolX
               xstar=(hibound+xhat(k))/2;
            else
               xstar=fzero(@(x) 2.*sum(sup.def2.WT2 ...
                  ./(sup.def2.tau'-repmat(x,[1 sup.def2.ntau])).^3,2), ...
                  [xhat(k) hibound]);
            end
         else
            error('unexpected case for dpsi')
         end
         
         % if min(phi) on the interval is <1/c then spectrum separates
         phi=sum(sup.def2.WT2./(sup.def2.tau' ...
            -repmat(xstar,[1 sup.def2.ntau])).^2,2)-1./sup.def2.c;
         
         if phi<0
            
            % find zero of phi-1/def.c between t(i) and xstar
            blo=t(k)+sqrt(w(k)*t(k)^2/(w(k)*t(k)^2/(t(k)-xstar)^2 ...
               -phi+sum(w(k+1:K).*t(k+1:K).^2 ...
               .*(1./(t(k+1:K)-xstar).^2-1./(t(k+1:K)-t(k)).^2),1)));
            separation(k,1)=fzero(@(u) sum(sup.def2.WT2./(sup.def2.tau' ...
               -repmat(u,[1 sup.def2.ntau])).^2,2)-1./sup.def2.c, ...
               [blo xstar]);
            
            % find zero of phi-1/def.c between xstar and t(i+1)
            bhi=t(k+1)-sqrt(w(k+1)*t(k+1)^2/(w(k+1)*t(k+1)^2/(t(k+1)-xstar)^2 ...
               -phi+sum(w(1:k).*t(1:k).^2 ...
               .*(1./(t(1:k)-xstar).^2-1./(t(1:k)-t(k+1)).^2),1)));
            separation(k,2)=fzero(@(u) sum(sup.def2.WT2./(sup.def2.tau' ...
               -repmat(u,[1 sup.def2.ntau])).^2,2)-1./sup.def2.c, ...
               [xstar bhi]);
            
            % keep eigenvalues weights
            separation(k,3)=sum(w(1:k))+zeroweight;
            separation(k,4)=sum(w(k+1:K));
         end
      end
   end
end

% find leftmost interval bound
expect2=sum(w.*t.^2);
leftlo=t(1)-sqrt(def.c*expect2)-1;
lefthi=t(1)-sqrt(def.c*w(1)*t(1)^2)/2;
leftmost=fzero(@(u) sum(sup.def2.WT2./(sup.def2.tau' ...
   -repmat(u,[1 sup.def2.ntau])).^2,2)-1./sup.def2.c,[leftlo lefthi]);

% find rightmost interval bound
rightlo=t(K)+sqrt(def.c*w(K)*t(K)^2)/2;
righthi=t(K)+sqrt(def.c*expect2)+1;
rightmost=fzero(@(u) sum(sup.def2.WT2./(sup.def2.tau' ...
   -repmat(u,[1 sup.def2.ntau])).^2,2)-1./sup.def2.c,[rightlo righthi]);

% construct support
separation(isnan(separation(:,1)),:)=[];
sup.numint=1+size(separation,1);
sup.u_Fbar=NaN+zeros(sup.numint,2);
sup.u_Fbar(1,1)=leftmost;
sup.u_Fbar(2:sup.numint,1)=separation(:,2);
sup.u_Fbar(sup.numint,2)=rightmost;
sup.u_Fbar(1:sup.numint-1,2)=separation(:,1);

% get modified Stieltjes transform evaluated at support endpoints
sup.m_LF=NaN+zeros(sup.numint,2);
sup.m_LF(:)=sum(repmat((w.*t)',[sup.numint*2 1]) ...
   ./(repmat(t',[sup.numint*2 1])-repmat(sup.u_Fbar(:),[1 K])),2);

% get support of sample eigenvalues distribution
sup.endpoint=max(0,sup.u_Fbar-def.c.*sup.u_Fbar.*sup.m_LF);
if def.p==def.n
   sup.endpoint(1,1)=0;
end

% get number of eigenvalues per interval in the support
sup.numeig=NaN+zeros(sup.numint,1);
if sup.numint>1
   sup.numeig(1)=separation(1,3)-max(zeroweight,1-def.n/def.p);
   sup.numeig(2:sup.numint-1)=diff(separation(1:sup.numint-1,3));
   sup.numeig(sup.numint)=separation(sup.numint-1,4);
   sup.numeig=round(sup.numeig.*def.p);
else
   sup.numeig=def.p-max(sum(def.tau==0),def.p-def.n);
end

% make results available to the gradient function
US_D.sup=sup;
US_D.def=def;

