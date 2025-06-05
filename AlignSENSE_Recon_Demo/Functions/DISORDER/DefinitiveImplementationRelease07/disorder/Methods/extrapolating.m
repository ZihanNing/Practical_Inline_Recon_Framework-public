function S=extrapolating(S,W,typ,mirr,gpu)

%EXTRAPOLATING   Extrapolates a coil map using the Papoulis-Grechberg 
%algorithm from [1] A Papoulis, "A new algorithm in spectral analysis and 
%band-limited extrapolation," IEEE Trans Circuits Syst, ASSP-22:735-742,
%Sep 1975 and [2] RW Gerchberg, "Super-resolution through error energy 
%reduction," Opt Acta, 21:709-720, Sep 1974. At the moment only 
%extrapolation in image domain is performed (based on [1]) but the plan is
%to extend this code to also perform superresolution (based on [2]).
%An alternative is to use Cadzow's algorithm [3] JA Cadzow, "An
%extrapolation procedure for band-limited signals," IEEE Trans Acoust, 
%Speech, Signal Processing, ASSP-27:4-12, Feb 1979 which avoids the
%truncation errors of [1] or to filter the data after performing the
%extrapolation. On top of this, the extrapolation matrix can be computed
%explicitly as suggested in [4] M Shaker, W Steenart, "An Approach to 
%Band-Limited Signal Extrapolation: The Extrapolation Matrix ," IEEE 
%Trans Circuits Syst, CAS-25(2):74-78, Feb 1978 and the required inversion
%can be accomplished by CG as in [5] T Strohmer, "On discrete band-limited 
%signal extrapolation," Contempor Math, Math Anal, Wavelets, Signal Process, 
%190:323-337, which is also implemented here in a multidimensional scenario 
%(and where the sampled area and band-pass of the signal may not be 
%separable). This provides exact solutions for the overdetermined problem 
%much quicker  than iterating and allows far reaching extrapolation (see 
%remark c at the end of Section 2 in [5]). Overdetermination is guaranteed 
%by proper construction of the bandpass filter. Neumann boundary conditions
%are promoted by using a DCT filtering strategy, as studied in [6] SA 
%Martucci, "Symmetric Convolution and the Discrete Sine and Cosine 
%Transforms," IEEE Trans Sign Process, 42(5):1038-1051, May 1994.
%   S=EXTRAPOLATING(S,W,{MIRR},{GPU})
%   * S are the sensitivities
%   * W is a mask that indicates the reliable coil estimates
%   * TYP is the algorithm used, 'PG' or 'C', 'C' by default
%   * {MIRR} indicates whether to use mirroring for extrapolation (defaults to 0)
%   * {GPU} determines whether to use gpu computations 
%   ** S are the extrapolated sensitivities
%

if ~exist('gpu','var') || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end;
if ~exist('mirr','var') || isempty(mirr);mirr=0;end
if ~exist('typ','var') || isempty(typ);typ='C';end

assert(ndims(W)<=3,'Number of dimensions for extrapolating (%d) should not be bigger than 3',ndims(W));

N=numel(W);
NP=gather(sum(W(:)));
assert(NP~=0,'No data is provided in input so no extrapolation is possible');
if N==NP;return;end

cutOffFact=[1 sqrt(2) 2];%Heuristic estimates of the filter bandwidth, but probably fine numerically
gibbsRingFact=1;

ND=numDims(W);vND=ones(1,ND);
harmV=[1 1/pi 3/(4*pi)];
cutOff=(NP*cutOffFact(ND)*harmV(ND))^(1/ND);%Cut off of the filter to preserve cutOffFact*NP samples when filtering

NF=(2*cutOff+1)*vND;
NS=size(W);
sp=NF./NS;

sp=sp/5;
%gibbsRingFact

H=buildFilter(NS*(1+mirr),'tukeyIso',sp,gpu,gibbsRingFact,mirr);
NS=size(S);NS(end+1:4)=1;
W=repmat(W,[1 1 1 NS(4:end)]);
S(~W)=0;

nIt=50;
tol=1e-3;
normConv=norm(S(W))*(N-NP)/NP;

if strcmp(typ,'PG')%Papoulis-Grechberg algorithm
    for n=1:nIt
        x=filtering(S,H,mirr);
        upS=norm(x(~W)-S(~W))/normConv;
        fprintf('Iteration %d. Update %.2g\n',n,upS);
        S(~W)=x(~W);
        if upS<tol;break;end
    end
    S=filtering(S,H);%To avoid truncation issues
elseif strcmp(typ,'C')%Cadzow's algorithm
    S0=filtering(S,H,mirr);S=S0;
    for n=1:nIt
       x=S;x(~W)=0;
       x=filtering(x,H,mirr);
       S=(S-x)+S0;
       upS=norm(S0(:)-x(:))/normConv;
       fprintf('Iteration %d. Update %.2g\n',n,upS);
       if upS<tol;break;end
    end
elseif strcmp(typ,'PG-CG')%Conjugate gradient Papoulis-Grechberg algorithm
    E.Xf=-(1-single(W));
    E.Gf=H;
    E.mi=mirr;E.ma=0;
    E.Ti.la=1;
    S=CGsolver(S,E,[],[],[],[],S);    
    S=filtering(S,H,mirr);%To avoid truncation issues
elseif strcmp(typ,'C-CG')%Conjugate gradient Cadzow's algorithm
    y=filtering(S,H,mirr);
    E.Xf=single(W);
    E.Ff{1}=H;
    E.mi=mirr;E.ma=0;
    S=CGsolver(y,E,[],[],[],[],y);
end
