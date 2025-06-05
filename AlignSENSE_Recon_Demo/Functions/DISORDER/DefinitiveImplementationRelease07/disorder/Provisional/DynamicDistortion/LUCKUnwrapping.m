function B=LUCKUnwrapping(x,M,modal)

% LUCKUNWRAPPING calls the Laplacian Unwrapping in Continuous K-space in 
%[1] MA Schofield, Y Zhu. "Fast phase unwrapping algorithm for 
%interferometric applications," Optics Letters, 28(14):1194-1196, July 
%2003 for dynamic field mapping
%   B=LUCKUNWRAPPING(X,M,MODAL)
%   * X is the complex data used to estimate the field
%   * W is a soft mask
%   * {MODAL} is the data modality (9->fMRI / 10->DWI), it defaults to 9
%   * B is the estimated field
%

if ~exist('modal','var') || isempty(modal);modal=9;end

gpu=isa(x,'gpuArray');

%UNWRAPPING PARAMETERS
N=size(x);
ND=numDims(x);
w=[5;10;2;1;0.5];w=w(1:ND);    
%w=[5;10;1;0.1;0.5];w=w(1:ND);
w=[1 1 1 1];
%fprintf('Used weights:%s\n',sprintf(' %.2f',w));

upSample=1;

%filterType='2ndFinite';
filterType='2ndFiniteDiscrete';

%PHASE UNWRAPPING ENCODING STRUCTURES
%Encoder
E.mi=1;%Mirror boundary conditions
E.ma=0;
Nup=round(N*upSample);
x=resampling(x,Nup);

E.Ff{1}=buildFilter(Nup*(1+E.mi),filterType,1./sqrt(w),gpu,0,E.mi);
%Decoder
EH=[];
%Constrain
C=[];
%Regularizer
R=[];
%R.Ti.la=1e-9;

%PHASE ENCODING PRECOMPUTATIONS
B=angle(x);
E.w=w;
if strcmp(filterType,'2ndFinite')
    B=cos(B).*filtering(sin(B),E.Ff{1},E.mi)-sin(B).*filtering(cos(B),E.Ff{1},E.mi);    
else
     Bn=B;Bn(:)=0;
     for n=1:ND;
         mirr=zeros(1,ND);mirr(n)=1;
         Baux=mirroring(B,mirr,1);
         siz=ones(1,ND);siz(n)=2*Nup(n);
         E.HF{n}=buildFilter(siz,'1stFiniteDiscreteForward',1,gpu);
         E.HB{n}=buildFilter(siz,'1stFiniteDiscreteBackward',1,gpu);        
         Baux=real(filtering(Baux,E.HF{n}));
         Baux=atan2(sin(Baux),cos(Baux));
         Baux=real(filtering(Baux,E.HB{n}));
         Baux=mirroring(Baux,mirr,0)/w(n);
         Bn=Bn+Baux;
     end
     B=Bn;
end

%Precondition
E.Ff{1}(1)=1e-9;%For inversion
A.Ps=E.Ff{1}.^(-1);
A.mi=E.mi;

%%PHASE UNWRAPPING SOLVER
B=CGsolver(B,E,EH,A,C,R);

%B=B+2*pi*round((Bn-B)/(2*pi));
%B=Bn;
B=resampling(B,N);



