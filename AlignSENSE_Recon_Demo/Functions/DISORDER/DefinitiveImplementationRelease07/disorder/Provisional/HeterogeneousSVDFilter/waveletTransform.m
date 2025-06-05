function x=waveletTransform(x,wA,di,T,h)

%WAVELETTRANSFORM performs the wavelet transform
%   X=WAVELETTRANSFORM(X,SH,{DI},{T},{H}) performs the wavelet transform
%   * X is the array to be transformed
%   * WA is a wavelet transform
%   * {DI} is the direction for transforming (1-> analysis / 0-> 
%   synthesis). Defaults to 1
%   * {T} is a thresholding array
%   * {H} indicates whether to use hard thresholding. Defaults to 0
%   ** X is the transformed array
%

if nargin<3 || isempty(di);di=1;end
if nargin<4 || isempty(T);T=0;end
if nargin<5 || isempty(h);h=0;end

if di
    for m=1:length(wA)
        if size(x,m)==size(wA{m},1);x=fftGPU(x,m,wA{m});else x=aplGPU(wA{m},x,m);end
    end
else
    if h==0;x=sign(x).*max(bsxfun(@minus,abs(x),T),0);%Soft thresholding
    elseif h==1;x(bsxfun(@le,abs(x),T))=0;
    else x=bsxfun(@times,x,T);
    end
    for m=1:length(wA)
        if size(wA{m},2)==size(x,m);x=ifftGPU(x,m,wA{m}');else x=aplGPU(wA{m}',x,m);end
    end
end
