function xou=shearletTransformGPUMemory(x,sH,di,T,h,re,or)

%SHEARLETTRANSFORMGPUMEMORY performs the shearlet transform considering
%potential GPU memory limitations
%   X=SHEARLETTRANSFORMGPUMEMORY(X,SH,{DI},{T},{H},{RE},{OR}) performs the 
%shearlet transform
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

NSH=size(sH.S);NSH(end+1:4)=1;
if di==1;assert(size(x,4)==1,'Array non-singleton for dimension 4 (size %d), which should be reserved for shearlets',size(x,4));
else assert(size(x,4)==size(sH.S,4),'Size of dimension 4 of the array (%d) not matching number of shearlets (%d)',size(x,4),size(sH.S,4));
end

gpu=single(gpuDeviceCount && ~blockGPU);
ND=numDims(sH.dfw);
if di==1
    for n=1:ND
        x=ifftshift(x,n);
        x=fftGPU(x,n);
    end
    xou=gather(x);
    xou=repmat(xou,[1 1 1 NSH(4)]);
    for s=1:NSH(4)
        sh=dynInd(sH.S,s,4);
        if gpu;sh=gpuArray(sh);end
        sh=bsxfun(@times,conj(sh),x);
        for n=1:ND
            sh=ifftGPU(sh,n);
            sh=fftshift(sh,n);
        end
        xou=dynInd(xou,s,4,gather(sh));
    end
    if re;xou=real(xou);end
else
    xou=dynInd(x,1,4);
    if gpu;xou=gpuArray(xou);end
    xou(:)=0;
    dfw=sH.dfw;
    if gpu;dfw=gpuArray(dfw);end
    for s=1:NSH(4)
        xs=dynInd(x,s,4);
        sh=dynInd(sH.S,s,4);
        if gpu;xs=gpuArray(xs);sh=gpuArray(sh);end
        if size(T,4)==NSH(4)
            Ts=dynInd(T,s,4);
            if gpu;Ts=gpuArray(Ts);end
        else
            Ts=T;
        end
        if h==0;xs=sign(xs).*max(bsxfun(@minus,abs(xs),Ts),0);%Soft thresholding
        elseif h==1;xs(bsxfun(@le,abs(xs),Ts))=0;
        else xs=bsxfun(@times,xs,Ts);
        end
        for n=1:ND
            xs=ifftshift(xs,n);
            xs=fftGPU(xs,n);
        end
        if di~=-1;x=bsxfun(@rdivide,bsxfun(@times,xs,sh),dfw);else xs=bsxfun(@times,xs,sh);end
        xou=xou+xs;
    end
    for n=1:ND
        xou=ifftGPU(xou,n);
        xou=fftshift(xou,n);%%%THERE WAS A BUG HERE---IT WILL HAVE AFFECTED THE EXPERIMENTS ON DENOISING
    end
    if re;xou=real(xou);end
end
