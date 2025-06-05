function p2=phaseVariance(x2,n2,ND,ph)

%PHASEVARIANCE computes the phase variance
%   P2=PHASEVARIANCE(X2,N2,{ND},{WI},{PH}) computes the phase variance
%   X2 is the squared observed envelope
%   N2 is the noise variance
%   {ND} are the number of dimensions, it defaults to 3
%   {PH} is the type of variance to compute (0 for trigonometric, 1 for
%   angular), defaults to 0
%   P2 is the estimated phase variance
%

if nargin<3;ND=[];end
if nargin<4 || isempty(ph);ph=0;end

gpu=isa(x2,'gpuArray');

N=size(x2);N(end+1:3)=1;
spacing=1;strength=0;
H=buildFilter(N(1:ND),'tukeyIso',spacing,gpu,strength);

x2=real(filtering(x2,H));
%n2=real(filtering(n2,H));
x2(n2<1e-6)=0;
x2=max(bsxfun(@minus,x2,n2),0);
K=x2./(n2+eps);
spacing=1;strength=0;
H=buildFilter(N(1:ND),'tukeyIso',spacing,gpu,strength);
K=max(real(filtering(K,H)),1e-2);
[~,p2]=phaseDistribution(0,K,ph);
for n=1:length(p2);p2{n}=reshape(p2{n},size(x2));end
%spacing=1/par;strength=1;
%H=buildFilter(N(1:ND),'tukeyIso',spacing,gpu,strength);
%p2{1}=max(real(filtering(p2{1},H)),1e-2);




