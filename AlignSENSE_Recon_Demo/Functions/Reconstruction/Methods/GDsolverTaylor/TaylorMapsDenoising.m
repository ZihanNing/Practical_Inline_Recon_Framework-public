
function [D, x, En] = TaylorMapsDenoising(D, W_m, W_t, MT, nIt, filt, lambda, debug)

N = size(D);
nD = ndims(D);

if nargin<2 || isempty(W_m); W_m = ones(N(1:3));end
if nargin<3 || isempty(W_t); W_t = ones(ones([1 length(N)]) );end
if nargin<4 || isempty(MT); MT =eye(4); warning('TaylorMapsDenoising:: MT not provided: will results in suboptimal performance.');end
if nargin<5 || isempty(nIt); nIt=50;end
if nargin<6 || isempty(filt); filt.sp=0.7; filt.gibs=.5;end
if nargin<7 || isempty(lambda); lambda.Smooth=1;lambda.Ridge=0.1;end
if nargin<8 || isempty(debug); debug=0;end

%%% ALGORITHM PARAMETERS
gpu = isa(D,'gpuArray');
H = buildFilter(N(1:3),'tukeyIso',filt.sp,gpu,filt.gibs);
%W_t = ones([1 1 1 2]); W_t = dynInd(W_t,i,4,0);
%W_m = (abs(W_m)/max(abs(W_m(:))));
toler = 0.01;
MT(1:3,4)=0; MS = sqrt( multDimSum((MT(1:3,1:3)).^2,1));

%%% LINEAR TERMS
FOV = N(1:3).*MS;
type = class(gather(x));
N=cast(N, type);
FOV=cast(FOV, type);
[Kx, Ky, Kz] = linearTerms(N(1:3), N(1:3)/2./FOV);%*N/2 to compensate for scaling in linearTerms.m and ./FOV to account for anisotropy
[Kx, Ky, Kz] = rotateCoordinates(Kx, Ky, Kz, inv(MT));

%%% MAKE KERNELS
Kroll = Kx.*Kz./(Kx.^2+Ky.^2+Kz.^2);
Kpitch = -Ky.*Kz./(Kx.^2+Ky.^2+Kz.^2);
K = cat(nD, Kroll, Kpitch); Kx=[];Ky=[];Kz=[];
K(isnan(K))=0;
K = bsxfun(@times, fftshift(H), K);
for i=1:3; K = ifftshiftGPU(K,i);end

%%% INITIALISE CG
%make it least square CG
D = bsxfun(@times, D, sqrt(W_m).*sqrt(W_t));%because of weighted LS
Dnew  = decode(D);%because of LS

%initialise variables
x = zeros(N(1:3),'like',real(D));%susceptibility distribution
EHEx = decode(encode(x));
r = Dnew - EHEx;% r = b - A * x;
p = r;
rsold = multDimSum(r.^2); %r' * r;

%%%RUN CG ITERATIONS
En = [];
for i = 1:nIt % https://en.wikipedia.org/wiki/Conjugate_gradient_method
    [Ep,Rp] = encode(p);
    EHEp= decode(Ep); [~,RHRp] = decode(Rp);
    
    Ap = EHEp + RHRp + lambda.Ridge*p ; %A * p;

    alpha = rsold / multDimSum( p.* Ap) ; % (p' * Ap);
    x = x + real(alpha * p);
    r = r - alpha * Ap;
    rsnew = multDimSum( r.^2); %r' * r;

    if multDimMax(abs(r))<toler; break;end
    beta = (rsnew / rsold);
    p = r + beta * p;
    rsold = rsnew;
    
    %%%REPORT RESIDUALS
    En = cat(2,En,rsold);
    if debug; figure(1);plot(En);end
end
    

%%% RETRUN RENOISED D - not encode because W_t and W_m will affect this
 %Fourier transform + kernel
for i=1:3; x=fftGPU(x,i);end
x = bsxfun(@times, x, K);
for i=1:3; x=ifftGPU(x,i);end
D=real(x);

%% FUNCTIONS
    function [Ex, Rx] = encode(x)%%%ENCODE
        %Fourier transform + kernel
        for q=1:3; x=fftGPU(x,q);end
        x = bsxfun(@times, x, K);
        for q=1:3; x=ifftGPU(x,q);end
        x=real(x);
        Ex = bsxfun(@times, x, sqrt(W_t).*sqrt(W_m));
        if nargout>1
            % Weighting + regularisation
            Rx = lambda.Smooth * (FiniteDiff(x,1,1) + ...
                                 FiniteDiff(x,2,1) + ...
                                 FiniteDiff(x,3,1) ) ;        

            Rx = bsxfun(@times, Rx, sqrt(W_m));  
        end
        
    end

    function [EHx, RHx] = decode(x)%%%DECODE
        % Weighting + regularisation
        x = bsxfun(@times, x, sqrt(W_m));
        if nargout>1
            RHx = lambda.Smooth * (FiniteDiff(x,1,0) + ...
                                 FiniteDiff(x,2,0) + ...
                                 FiniteDiff(x,3,0) ) ;
            %Fourier transform + kernel
            for q=1:3; RHx=fftGPU(RHx,q);end
            RHx = bsxfun(@times, RHx, K);
            for q=1:3; RHx=ifftGPU(RHx,q);end
            %Sum
            RHx = real( multDimSum(RHx,nD));
            EHx=[];
        else
            x = bsxfun(@times, x, sqrt(W_t));
            %Fourier transform + kernel
            for q=1:3; x=fftGPU(x,q);end
            x = bsxfun(@times, x, K);
            for q=1:3; x=ifftGPU(x,q);end
            %Sum
            EHx = real(multDimSum(x,nD));
            RHx = [];
        end


    end
   

    
end