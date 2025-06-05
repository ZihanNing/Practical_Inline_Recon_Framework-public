% function dlambda=QuESTgrad(tau,n)
%
% Gradient of Quantized Eigenvalues Sampling Transform
%
% See article "Spectrum Estimation: A Unified Framework for Covariance 
% Matrix Estimation and PCA in Large Dimensions" by Olivier Ledoit and 
% Michael Wolf, Journal of Multivariate Analysis (July 2015), volume 139,
% pages 260-384, Section 2.2, equations (2.11)-(2.17)
%
% Inputs: tau   : column vector containing population eigenvalues
%         n     : positive scalar integer containing sample size
%
% Output: dlambda: matrix whose entry (i,j) is the partial derivative
%                   of the ith sample eigenvalue with respect to the 
%                   jth population eigenvalue
%
% To see help message, type: QuESTgrad('help')
%    
% Copyright Olivier Ledoit and Michael Wolf (2013-2016)

function dlambda=QuESTgrad(tau,n)

if nargin==1
   if isequal(tau,'help')
      disp('  function dlambda=QuESTgrad(tau,n)')
      disp(' ')
      disp('  Gradient of Quantized Eigenvalues Sampling Transform')
      disp(' ')
      disp('  See article "Spectrum Estimation: A Unified Framework for Covariance')
      disp('  Matrix Estimation and PCA in Large Dimensions" by Olivier Ledoit and')
      disp('  Michael Wolf, Journal of Multivariate Analysis (July 2015), volume 139,')
      disp('  pages 260-384, Section 2.2, equations (2.11)-(2.17)')
      disp(' ')
      disp('  Inputs : tau   : column vector containing population eigenvalues')
      disp('           n     : positive scalar integer containing sample size')
      disp(' ')
      disp('  Output: dlambda: matrix whose entry (i,j) is the partial derivative')
      disp('                    of the ith sample eigenvalue with respect to the')
      disp('                    jth population eigenvalue')
      disp(' ')
      disp('  To see help message, type: QuESTgrad(''help'')')
      disp(' ')
      disp('  Copyright Olivier Ledoit and Michael Wolf (2013-2016)')
      disp(' ')
   end
else
   % define inputs
   def.tau=tau;
   def.n=n;
   def=checkinputs01(def);

   % compute gradient by calling the subfunctions
   [dint,ddis,dden,dsol,dnet,dsup,def]=intgrd02(def);

   % post-process the output
   dlambda=dint.lambda;
   
   % clear global variable
   global US_D
   US_D=[];
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
   disp(min(def.tau))
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

% function [dint,ddis,dden,dsol,dnet,dsup,def]=intgrd02(def)
%
% compute gradient of sample eigenvalues wrt population eigenvalues 
%
% inputs:
%   def.tau     : column vector of dimension p with location of Diracs
%   def.n       : number of observations
%   def.p       : length of def.tau
%   def.c       : p/n
%   def.weight  : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   dint.lambda : gradient of sample eigenvalues wrt population eigenvalues 
%
% uses global variable US_D set by intfun01
%
% direct dependency : disgrd03
%
% Version 02 : use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

% supporting calculations: Olivier's notepad 25 April, 9 June & 24 Nov 2012

function [dint,ddis,dden,dsol,dnet,dsup,def]=intgrd02(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14) ...
      |(~isequal(def.weight,US_D.def.weight)) ...
      |(~isequal(def.p,US_D.def.p))|(~isequal(def.c,US_D.def.c))
   error('incompatible definitions of the problem')
end

% get gradient of zeta and G(zeta)
[ddis,dden,dsol,dnet,dsup,def]=disgrd03(def);

% initialize
dint.lambda=NaN+zeros(def.p,def.p);
p=def.p; % to shorten some formulas

% take care of null sample eigenvalues 
n0tau=sum(def.tau==0); % number of null population eigenvalues
if def.p-def.n>n0tau
   dint.lambda(1:def.p-def.n,:)=0;
else
   if n0tau>0
      % reference: Research Olivier\Regularization\Matlab\testgrad2.m
      dint.lambda(max(1,def.p-def.n+1):n0tau,1:n0tau)=1-def.c; 
      dint.lambda(max(1,def.p-def.n+1):n0tau,n0tau+1:def.p)=0;
   end
end

% loop over support intervals
for i=1:US_D.sup.numint
   
   % prepare variables relative to this interval
   idx=find(US_D.dis.interval==i);
   F=US_D.dis.G(idx);
   [F,i1]=unique(F);
   idx=idx(i1);
   x=US_D.dis.x(idx);
   nidx=length(idx);
   dF=ddis.G(idx,:);
   dx=US_D.den.a.*repmat(US_D.dis.zeta(idx,:).^(US_D.den.a-1),[1 def.p]) ...
      .*ddis.zeta(idx,:);
   k=US_D.int.bin{i};
   q=US_D.int.quant{i};
   nq=US_D.int.nquant(i);
   
   % compute gradient of integral of x from F(k) to F(k+1)
   dintxdF=(dF(2:nidx,:)-dF(1:nidx-1,:)) ...
      .*repmat((x(1:nidx-1)+x(2:nidx))./2,[1 def.p]) ...
      +repmat((F(2:nidx)-F(1:nidx-1)),[1 def.p]) ...
      .*(dx(1:nidx-1,:)+dx(2:nidx,:))./2;

   % compute gradient of integral of x from F(k) to quant(i)
   dintxdF2=repmat(q,[1 p]).*dx(k,:) ...
      -dx(k,:).*repmat(F(k),[1 p])-dF(k,:).*repmat(x(k),[1 p]) ...
      -2.*dF(k,:).*repmat((q-F(k)).*(x(k+1)-x(k))./(F(k+1)-F(k))./2,[1 p]) ...
      +(dx(k+1,:)-dx(k,:)).*repmat((q-F(k)).^2./(F(k+1)-F(k))./2,[1 p]) ...
      -(dF(k+1,:)-dF(k,:)).*repmat((q-F(k)).^2.*(x(k+1)-x(k))./(F(k+1)-F(k)).^2./2,[1 p]);
   
   % compute integral of x over the interval [F(1),quant(i)]
   dintxdF3=US_D.int.indic{i}*dintxdF+dintxdF2;
   
   % compute sample eigenvalues
   dint.lambda(round(F(1)*def.p+1):round(F(nidx)*def.p),:)= ...
      diff(dintxdF3,1,1).*def.p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [ddis,dden,dsol,dnet,dsup,def]=disgrd03(def)
%
% compute gradients of x and F(x) 
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   ddis.zeta : gradient of the power-transformed abscissa
%   ddis.G    : gradient of the power-transformed cdf 
%
% uses global variable US_D set by disfun04
%
% direct dependency : dengrd02
%
% Version 03 : use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

function [ddis,dden,dsol,dnet,dsup,def]=disgrd03(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14)|(~isequal(def.weight,US_D.def.weight)) ...
      |(~isequal(def.p,US_D.def.p))|(~isequal(def.c,US_D.def.c))
   error('incompatible definitions of the problem')
end

% get the gradient of -1/m_Fbar(z)
[dden,dsol,dnet,dsup,def]=dengrd02(def);

% initialize variables
ddis.zeta=NaN+zeros(US_D.dis.numx,def.p);
ddis.G=NaN+zeros(US_D.dis.numx,def.p);

% loop over support intervals
for i=1:US_D.sup.numint
   
   % compute gradient of dis.zeta
   indx=(US_D.dis.interval==i);
   intlength=sum(indx);
   isin=(US_D.net.interval==i);
   ddis.zeta(indx,:)=[reshape((1/US_D.den.a).*US_D.sup.endpoint(i,1,:) ...
      .^(1/US_D.den.a-1).*dsup.endpoint(i,1,:),[1 def.p]); ...
      dden.zeta(isin,:);reshape((1/US_D.den.a).*US_D.sup.endpoint(i,2,:) ...
      .^(1/US_D.den.a-1).*dsup.endpoint(i,2,:),[1 def.p])];
   
   % special treatment for left-hand side bound of the support when p=n
   if indx(1)
      if def.p==def.n
         ddis.zeta(1,:)=0;
      end
   end

   % compute gradient of dis.G
   g=[0;US_D.den.g(isin);0];
   dg=[zeros(1,def.p);dden.g(isin,:);zeros(1,def.p)];
   dG_=cumsum(diff(ddis.zeta(indx,:),1,1) ...
      .*repmat((g(1:intlength-1)+g(2:intlength,:))./2,[1 def.p]),1) ...
      +cumsum(repmat(diff(US_D.dis.zeta(indx),1,1),[1 def.p]) ...
      .*(dg(1:intlength-1,:)+dg(2:intlength,:))./2,1);
   dGraw=[zeros(1,def.p);dG_];
   
   % force value of F at interval endpoint to be the correct one
   wdif=US_D.dis.toweight(i)-US_D.dis.fromweight(i);
   ddis.G(indx,:)=wdif.*dGraw./US_D.dis.Graw{i}(intlength) ...
      -wdif.*repmat(US_D.dis.Graw{i},[1 def.p]) ...
      .*repmat(dGraw(intlength,:),[intlength 1]) ...
      ./US_D.dis.Graw{i}(intlength).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dden,dsol,dnet,dsup,def]=dengrd02(def)
%
% compute gradients of x and f(x) density of sample eigenvalues wrt weights
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   dden.x   : gradient of the abscissa 
%   dden.f   : gradient of the density of sample eigenvalues
%   dden.m_LF: gradient of modified Stieltjes transform evaluated at x
%   dden.zeta : gradient of the power-transformed abscissa 
%   dden.g   : gradient of the power-transformed density 
%
% uses global variable US_D set by denfun01
%
% direct dependency : solgrd01
%
% Version 02: use zeta=x^(1/a) and G(zeta)=F(zeta^a) instead of x and F

% supporting calculations: Olivier's notepad 24 April, 10 May & 10 June 2012 

function [dden,dsol,dnet,dsup,def]=dengrd02(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14)|(~isequal(def.weight,US_D.def.weight)) ...
      |(~isequal(def.p,US_D.def.p))|(~isequal(def.c,US_D.def.c))
   error('incompatible definitions of the problem')
end

% get the gradient of -1/m_Fbar(z)
[dsol,dnet,dsup,def]=solgrd01(def);

% make variables of the right dimensions
WEIGHT = repmat(def.weight',[US_D.net.nxi 1]);
TAU    = repmat(def.tau',[US_D.net.nxi 1]);
Z      = repmat(US_D.sol.zxi,[1 def.p]);

% compute gradient of f
dden.f=1./pi./def.c.*imag(dsol.zxi./Z.^2);

% compute gradient of m_LF
dden.m_LF=-Z.*WEIGHT./(TAU-Z).^2 ...
   +dsol.zxi.*repmat(sum(WEIGHT.*TAU./(TAU-Z).^2,2),[1 def.p]);

% compute the gradient of x
dden.x=real(dsol.zxi.*repmat(1-def.c.*US_D.den.m_LF,[1 def.p]) ...
   -def.c.*Z.*dden.m_LF);

% compute the gradient of zeta
dden.zeta=(1/US_D.den.a).*dden.x ...
   .*repmat(US_D.den.x.^(1/US_D.den.a-1),[1 def.p]);

% compute gradient of g
dden.g=US_D.den.a.*(US_D.den.a-1).*dden.zeta ...
   .*repmat(US_D.den.f.*US_D.den.zeta.^(US_D.den.a-2),[1 def.p]) ...
   +US_D.den.a.*repmat(US_D.den.zeta.^(US_D.den.a-1),[1 def.p]).*dden.f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dsol,dnet,dsup,def]=solgrd01(def)
%
% compute derivative of -1/m_Fbar(z) with respect to the tau's
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.n      : number of observations
%   def.p      : length of def.tau
%   def.c      : p/n
%   def.weight : p*1 vector where each entry is equal to 1/p
%
% outputs:
%   dsol.zxi   : gradient of -1/m_Fbar(z) with respect to the tau's
%                dimension net.nxi*def.p (complex)
%
% uses global variable US_D set by solfun02
%
% direct dependency : netgrd02

% supporting calculations: Olivier's notepad 21-24 April and 9 May 2012

function [dsol,dnet,dsup,def]=solgrd01(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14)|(~isequal(def.weight,US_D.def.weight)) ...
      |(~isequal(def.p,US_D.def.p))|(~isequal(def.c,US_D.def.c))
   error('incompatible definitions of the problem')
end

% get the gradient of the real part of -1/m_Fbar(z)
[dnet,dsup,def]=netgrd02(def);

% make variables of the right dimensions
WEIGHT = repmat(def.weight',[US_D.net.nxi 1]);
TAU    = repmat(def.tau',[US_D.net.nxi 1]);
XI     = repmat(US_D.net.xi,[1 def.p]);
YXI    = repmat(imag(US_D.sol.zxi),[1 def.p]);

% compute intermediary quantities
t_x         = TAU-XI;
cplxnorm    = t_x.^2+YXI.^2;
cplxnorm2   = cplxnorm.^2;
wt2         = WEIGHT.*TAU.^2;
numerator   = 2.*sum(wt2.*t_x./cplxnorm2,2);
denominator = 2.*sum(wt2.*YXI./cplxnorm2,2);

% compute d(imag(zxi))/d(tau)
dyxi_dtau=(2.*WEIGHT.*TAU./cplxnorm-2.*wt2.*t_x./cplxnorm2) ...
   ./repmat(denominator,[1 def.p]);

% compute d(imag(zxi))/d(real(zxi))
dyxi_dxi=numerator./denominator;

% put it all together
dyxi=dyxi_dtau+repmat(dyxi_dxi,[1 def.p]).*dnet.xi;

% reconstruct the derivative of the complex number -1/m_Fbar(z)
dsol.zxi=dnet.xi+sqrt(-1).*dyxi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dnet,dsup,def]=netgrd02(def)
%
% compute derivative of the grid over [u_1,u_2] which is the real part 
% of -1/m_Fbar(z) with respect to the weights
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
%   dnet.xi    : gradient of evenly spaced grid over the image of the 
%                support of F or Fbar through the function -1/m_Fbar
%                dimension nz*q*def.p
%
%   see also function supgrd01.m for outputs in the dsup structure
%
% uses global variable US_D set by netfun02
%
% direct dependency : supgrd01
%
% Version 02 : tailor xi grid to spectrum separation

function [dnet,dsup,def]=netgrd02(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14)|(~isequal(def.p,US_D.def.p))
   error('incompatible definitions of the problem')
end

% get the gradient of the support
[dsup,def]=supgrd01(def);

% compute the gradient of the grid
dnet.xi=US_D.net.sup2xi*reshape(dsup.u_Fbar,[US_D.sup.numint*2 def.p]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dsup,def]=supgrd01(def)
%
% compute gradient of support in u=-1/m_Fbar(z) space
%
% inputs:
%   def.tau    : column vector of dimension p with location of Diracs
%   def.p      : length of def.tau
%
% outputs:
%   dsup.u_Fbar  : gradient of -1/m_Fbar(x) at support endpoints
%
% direct dependency : none
%
% Version 01: differentiate with respect to tau instead of weights
%
% uses global variable US_D set by supfun06

% Supporting calculations in Olivier's notepad on 7 & 9 May 2012

function [dsup,def]=supgrd01(def)

global US_D

% check problem definition is compatible with global variables
if (max(abs(def.tau-US_D.def.tau))>1e-14)|(~isequal(def.p,US_D.def.p))
   error('incompatible definitions of the problem')
end

% initialize matrices
TAU=repmat(shiftdim(def.tau,-2),[US_D.sup.numint 2 1]);
U=repmat(US_D.sup.u_Fbar,[1 1 def.p]);
WEIGHT=repmat(shiftdim(def.weight,-2),[US_D.sup.numint 2 1]);

% compute derivative of support endpoints in u-space with respect to tau
dsup.u_Fbar=WEIGHT.*TAU.*U./(TAU-U).^3 ...
   ./(repmat(sum(WEIGHT.*TAU.^2./(TAU-U).^3,3),[1 1 def.p]));

% compute derivative of the endpoints of the support of F
dsup.endpoint=def.c.*U.^2.*WEIGHT./(TAU-U).^2;
if def.p==def.n
   dsup.endpoint(1,1,:)=0;
end
