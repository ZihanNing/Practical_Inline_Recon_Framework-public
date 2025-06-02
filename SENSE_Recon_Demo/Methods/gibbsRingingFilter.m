function x=gibbsRingingFilter(x,sub,fact,co,di,hard)

% GIBBSRINGINGFILTER applies a Gibbs ringing filter to smooth
%   X=GIBBSRINGINGFILTER(X,{SUB},{FACT},{CO},{DI})
%   * X is the data to be filtered
%   * {SUB} is the subsampling factor for filter resolution (see
%   buildFilter.m)
%   * {FACT} is the filter shape (see buildFilter.m)
%   * {CO} serves to apply the filter in the cosine domain
%   * {DI} is the type of filtering, 1 (default) for low pass, 0 for high
%   pass
%   * {HARD} indicates hard filtering
%   ** X is the filtered data
%

if nargin<2 || isempty(sub);sub=1;end
if nargin<3 || isempty(fact);fact=1;end
if nargin<4 || isempty(co);co=1;end
if nargin<5 || isempty(di);di=1;end
if nargin<6 || isempty(hard);hard=0;end

gpu=isa(x,'gpuArray');

NX=size(x);NX(end+1:3)=1;NX=NX(1:3);
NX=NX*(1+co);
if length(fact)==1;typ='tukeyIso';else typ='tukey';end
H=buildFilter(NX,typ,sub,gpu,fact,co);
if hard;H(H<1)=0;end
if ~di;H=1-H;end
x=filtering(x,H,co);
