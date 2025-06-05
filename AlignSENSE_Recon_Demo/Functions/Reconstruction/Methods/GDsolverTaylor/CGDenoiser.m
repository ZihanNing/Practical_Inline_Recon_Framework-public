
function [Ddenoised]= CGDenoiser(D, W)

nIter = 50;
toler = 0.01;
lambda_smooth = 10;

N = size(D);
nD = ndims(D);

% MAKE KERNELS
[Kx, Ky, Kz]=linearTerms(N(1:3));
Kroll = Kx.*Kz./(Kx.^2+Ky.^2+Kz.^2);
Kpitch = -Ky.*Kz./(Kx.^2+Ky.^2+Kz.^2);
K = cat(nD, Kroll, Kpitch);
K(isnan(K))=0;

%%% INITIALISE
    %make it least square CG
    Dnew = decode(D,K,W);
    
    x = zeros(N(1:3),'like',real(D));
    r = Dnew - decode(encode(x, K,W),K,W);% r = b - A * x;
    p = r;
    rsold = multDimSum(r.^2); %r' * r;
    
%%%RUN

for i = 1:nIter
    
    Ap = decode( encode(p, K,W),K,W)  ; %A * p;
    Ap = Ap +  lambda_smooth * ( FiniteDiff( FiniteDiff(p,1),1 ,0) + ...
                                 FiniteDiff( FiniteDiff(p,2),2 ,0) + ...
                                 FiniteDiff( FiniteDiff(p,3),3 ,0) ) / numel(p);%take into account smoothness regularisation
    alpha = rsold / multDimSum( p.* Ap) ; % (p' * Ap);
    x = x + real(alpha * p);
    r = r - alpha * Ap;
    rsnew = multDimSum( r.^2); %r' * r;

    if multDimMax(abs(r))<toler; break;end

    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
end
    
    
Ddenoised = encode(x, K, W); %denoised version of D

% plot_(W,[],[],dynInd(D,1,6), dynInd(Ddenoised,1,6),[],[],[],[],[],33);
% plot_(x,[],[],dynInd(D,2,6), dynInd(Ddenoised,2,6),[],[],[],[],[],34);

end


function x = encode(x, K, W)%%%ENCODE
    for i=1:3; x=fftGPU(x,i);end
    
    x = bsxfun(@times, x, K);
    
    for i=1:3; x=ifftGPU(x,i);end
    x=real(x);
    x = bsxfun(@times, x, W);

end

function x = decode(x, K, W)%%%DECODE
    x = bsxfun(@times, x, W);
    for i=1:3; x=fftGPU(x,i);end
    
    x = bsxfun(@times, x, K);
    
    for i=1:3; x=ifftGPU(x,i);end
    x = real( multDimSum(x,6));
    
end