function x=shearletTransform(x,sH,di,T,h,re,or)

%SHEARLETTRANSFORM performs the shearlet transform
%   X=SHEARLETTRANSFORM(X,SH,{DI},{T}) performs the shearlet transform
%   X is the array to be transformed
%   SH is a shearlet structure
%   {DI} is the direction for transforming (1-> analysis / 0-> synthesis).
%   Defaults to 1
%   {T} is a thresholding array
%   {H} indicates whether to use hard thresholding. Defauts to 0
%   {RE} indicates whether to take the real component
%   {OR} indicates whether to normalize the forward transform
%   X is the transformed array
%

if nargin<3 || isempty(di);di=1;end
if nargin<4 || isempty(T);T=0;end
if nargin<5 || isempty(h);h=0;end
if nargin<6 || isempty(re);re=0;end
if nargin<7 || isempty(or);or=0;end

if di==1;assert(size(x,4)==1,'Array non-singleton for dimension 4 (size %d), which should be reserved for shearlets',size(x,4));
else assert(size(x,4)==size(sH.S,4),'Size of dimension 4 of the array (%d) not matching number of shearlets (%d)',size(x,4),size(sH.S,4));
end

ND=numDims(sH.dfw);
if di==1
    for n=1:ND
        x=ifftshift(x,n);
        x=fftGPU(x,n);
    end
    x=bsxfun(@times,x,conj(sH.S));
    for n=1:ND
        x=ifftGPU(x,n);
        x=fftshift(x,n);
    end
    if re;x=real(x);end
else
    if h==0;x=sign(x).*max(bsxfun(@minus,abs(x),T),0);%Soft thresholding
    elseif h==1;x(bsxfun(@le,abs(x),T))=0;
    else x=bsxfun(@times,x,T);
    end
    for n=1:ND
        x=ifftshift(x,n);
        x=fftGPU(x,n);
    end
    if di~=-1;x=bsxfun(@rdivide,sum(bsxfun(@times,x,sH.S),4),sH.dfw);else x=sum(bsxfun(@times,x,sH.S),4);end
    for n=1:ND
        x=ifftGPU(x,n);
        x=fftshift(x,n);%%%THERE WAS A BUG HERE---IT WILL HAVE AFFECTED THE EXPERIMENTS ON DENOISING
    end
    if re;x=real(x);end
end
