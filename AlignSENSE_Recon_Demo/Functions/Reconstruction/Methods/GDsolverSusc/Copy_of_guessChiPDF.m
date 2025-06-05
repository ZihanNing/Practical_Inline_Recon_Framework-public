

function [chi, B0, B0Est] = guessChiPDF(chiInit, x, B0GT, TE, MT, padDim, WDataConsist, nIt, filt, Reg, H0, gamma, debug)

N = size(x);
nD = numDims(x);

if nargin<4 || isempty(TE); TE=1; warning('guessChiPDF:: TE not provided: will results in suboptimal performance.'); end
if nargin<5 || isempty(MT); MT=eye(4); warning('guessChiPDF:: MT not provided: will results in suboptimal performance.'); end
if nargin<6 || isempty(padDim); padDim=[0 0 0]; end
if nargin<7 || isempty(WDataConsist); WDataConsist=abs(x); end
if nargin<8 || isempty(nIt); nIt=50; end
if nargin<9 || isempty(filt); filt.sp=0.7; filt.gibs=.4; end
if nargin<10 || isempty(Reg); Reg.Smooth=0;Reg.Ridge=0;Reg.Tissue=0; end
if nargin<11 || isempty(H0); H0=7;warning('guessChiPDF:: H0 not provided. Assumed to be 7T.'); end%In T
if nargin<12 || isempty(gamma); gamma=42.58;warning('guessChiPDF:: gamma not provided. Assumed to be 42.58MHz/T.'); end%MHz/T
if nargin<13 || isempty(debug); debug=2; end

gpu = isa(x,'gpuArray');
MS = sqrt(sum(MT(1:3,1:3).^2,1));

if ~isfield(Reg,'Smooth');Reg.Smooth=0;end
if ~isfield(Reg,'WSmooth') && Reg.Smooth>0;Reg.WSmooth=ones(size(x));end

if ~isfield(Reg,'Ridge');Reg.Ridge=0;end
if ~isfield(Reg,'WRidge') && Reg.Ridge>0;Reg.WRidge=ones(size(x));end

if ~isfield(Reg,'Tissue');Reg.Tissue=0;end
if ~isfield(Reg,'WTissue')&& Reg.Tissue>0;Reg.WTissue=ones(size(x));end

if ~isfield(Reg,'tissueSusc');Reg.tissueSusc=-9.2;end
WPTFI = ones(size(x));


%%% Convert image phase to B0 in Hz
B0 = CNCGUnwrapping(x,MS,'MagnitudeGradient4','LSIt')/2/pi/TE;
if debug>1
    phOrig = angle(x);
    phPred = angle(exp(+1i*2*pi*TE*B0));
    rr = [-200 200];
    plotND({abs(x),1},cat(4,B0,rescaleND(phOrig,rr),rescaleND(phPred,rr)),rr,[],0,{[],2},MT,{'Estimated B0';'Image phase';'B0 predicted phase'});
    plotND({abs(x),1},cat(4,B0GT,B0),rr,[],0,{[],2},MT,{'GT B0';'Estimated B0'});
end

%%% Delete lower order SH to remove the physiological B0
addpath(genpath('/home/ybr19/Projects/dB0_analysis'));
orderRemove = 1;
B = SH_basis(N,orderRemove,1);
%[~,~,B0] = interpolateBasis(B0, B);
%B0=B0GT;

usePTFI = 0;

% Padding
x  = padArrayND(x, padDim, [], 0);
WDataConsist  = padArrayND(WDataConsist, padDim, [], 0);
if isfield(Reg,'WSmooth'); Reg.WSmooth = padArrayND(Reg.WSmooth, padDim, [], 0);end
if isfield(Reg,'WRidge'); Reg.WRidge = padArrayND(Reg.WRidge, padDim, [], 0);end
if isfield(Reg,'WTissue'); Reg.WTissue = padArrayND(Reg.WTissue, padDim, [], 0);end

WPTFI = padArrayND(WPTFI, padDim, [], 0);
B0  = padArrayND(B0, padDim, []);%what you pad should not matter since weighting is set to 0
[MS,MT] = mapNIIGeom(MS,MT,'padArrayND',padDim, N); 
N = size(x);

if debug>1; plotND({abs(x),1},cat(4,WDataConsist),[],[],0,{[],2},MT,{'Image';'Mask'});end

% Create kernels and filters
FOV = N(1:3).*MS;
type = class(gather(x));
N=cast(N, type);
FOV=cast(FOV, type);
[kx, ky, kz] = linearTerms(N(1:3), N(1:3)/2./FOV);%*N/2 to compensate for scaling in linearTerms.m and ./FOV to account for anisotropy
[kx, ky, kz] = rotateCoordinates(kx, ky, kz, inv(MT));

K = 1/3 - kz.^2./(kx.^2+ky.^2+kz.^2);
K(isnan(K))=1/3; %If 0, mean induced field=0, if 1/3, mean induced field = 1/3*mean(susc)

H = buildFilter(N(1:3),'tukeyIso',filt.sp,gpu,filt.gibs);
H = ones(size(kz),'like', kz);

K = bsxfun(@times, fftshift(H), K);
for i=1:3; K = ifftshiftGPU(K,i);end

%%% INITIALISE CG
%make it least square CG
bNew = bsxfun(@times, B0, (WDataConsist.^2));%because of weighted LS - weighting in encode is applied in the encode function
bNew = decode(bNew);%because of LS
if Reg.Tissue>0; bNew = bNew + Reg.Tissue * Reg.WTissue.^2* Reg.tissueSusc;end

%initialise variables
if isempty(chiInit); chi = zeros(N(1:3),'like',real(x));else chi = chiInit;end %susceptibility distribution
EHEx = decode( WDataConsist.^2 .* encode(chi));
if Reg.Smooth>0 
    [~,Rx] = encode(chi);
    [~,RHRx] = decode(Reg.WSmooth.^2.*Rx);
else
    RHRx=0;
end
if Reg.Tissue>0; THTx = Reg.WTissue.^2.*chi;else THTx=0;end
r = bNew - (EHEx + Reg.Smooth * RHRx + Reg.Tissue*THTx);% r = b - A * x;
p = r;
rsold = multDimSum(r.^2); %r' * r;
toler = .0;

%%%RUN CG ITERATIONS
En = [];
for i = 1:nIt % https://en.wikipedia.org/wiki/Conjugate_gradient_method
    EHEp = decode(WDataConsist.^2.*encode(p));
    if Reg.Smooth>0
        [~,Rp] = encode(p);
        [~,RHRp] = decode(Reg.WSmooth.^2.*Rp);
    else
        RHRp=0;
    end    
    if Reg.Tissue>0; THTp = Reg.WTissue.^2.*p;else THTp=0;end

    Ap = EHEp + Reg.Smooth* RHRp + Reg.Tissue* THTp ; %A * p;

    alpha = rsold / multDimSum( p.* Ap) ; % (p' * Ap);
    chi = chi + real(alpha * p);
    r = r - alpha * Ap;
    rsnew = multDimSum( r.^2); %r' * r;

    if multDimMax(abs(r))<toler; break;end
    beta = (rsnew / rsold);
    p = r + beta * p;
    rsold = rsnew;
    
    %%%REPORT RESIDUALS
    En = cat(2,En,rsold);
    if debug; figure(189);plot(En);end
    if debug>2; plotND({abs(x),1},chi,[],[],0,{[],2},MT,'Chi',[],[],89);end

end
    
if debug>1; plotND({abs(x),1},chi,[],[],0,{[],2},MT,'Chi');end

%%% RETRUN B0 
B0Est = encode(chi);%Ok since Weighting not included in encode

%%% REMOVE PADDING
x  = padArrayND(x, padDim, 0);
B0  = padArrayND(B0, padDim, 0);
B0Est  = padArrayND(B0Est, padDim, 0);
chi  = padArrayND(chi, padDim, 0);

if debug>1
    rr = [];
    plotND({abs(x),1},cat(4,B0,B0Est, abs(B0-B0Est) ),rr,[],1,{[],2},MT,{'Estiamted B0 from image phase';'Susceptibility predicted B0';'Error'},abs(x)>.3*multDimMax(abs(x)),{3});
end

%% FUNCTIONS
    function [Ex, Rx] = encode(x)%%%ENCODE
        if usePTFI; x=bsxfun(@times, x, WPTFI);end
        if nargout>1
            % Weighting + regularisation
            Rx = (FiniteDiff(x,1,1) + ...
                                  FiniteDiff(x,2,1) + ...
                                  FiniteDiff(x,3,1) ) ;        

            %Rx = bsxfun(@times, Rx, sqrt(W));  
        end
        %Fourier transform + kernel
        for q=1:3; x=fftGPU(x,q);end
        x = bsxfun(@times, x, K);
        for q=1:3; x=ifftGPU(x,q);end
        x=real(x);
        Ex = H0*gamma*x;%bsxfun(@times, x, sqrt(W));
    end

    function [EHx, RHx] = decode(x)%%%DECODE
        % Weighting + regularisation
        %x = bsxfun(@times, x, sqrt(W));
        if nargout>1
            RHx =  (FiniteDiff(x,1,0) + ...
                    FiniteDiff(x,2,0) + ...
                    FiniteDiff(x,3,0) ) ;
        else
            RHx=[];
        end
        %Fourier transform + kernel
        for q=1:3; x=fftGPU(x,q);end
        x = bsxfun(@times, x, K);
        for q=1:3; x=ifftGPU(x,q);end
        if usePTFI; x=bsxfun(@times, x, WPTFI);end
        EHx = H0*gamma*real(x);
    end
   

    
end












