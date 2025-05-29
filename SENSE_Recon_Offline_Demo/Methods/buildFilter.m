function x=buildFilter(N,typ,sp,gpu,gibbsRing,c)

% BUILDFILTER builds a filter in the frequency domain
%   X=BUILDFILTER(N,TYP,{SP},{GPU},{GIBBSRING},{C})
%   * N are the dimensions of the filter
%   * TYP is the filter type
%   * {SP} is the physical spacing
%   * {GPU} determines whether to use gpu computations 
%   * {GIBBSRING} is a filter parameter to prevent Gibbs ringing, defaults 
%   to 0
%   * {C} indicates if the filter is to be applied in the Cosine domain 
%   (for Neumann instead of periodic boundary conditions), defaults to 0, 
%   i.e., the filter is applied in the Fourier domain
%   ** X is the filter profile
%

if nargin<3 || isempty(sp);sp=1;end
if nargin<4 || isempty(gpu);gpu=useGPU;end
if nargin<5 || isempty(gibbsRing);gibbsRing=0;end
if nargin<6 || isempty(c);c=0;end

ND=length(N);
if length(c)==1;c=c*ones(1,ND);end
if length(sp)==1;sp=sp*ones(1,ND);end
x=0;
if gpu;x=gpuArray(x);end
%%NOT SURE THAT THE SPACING IN KGRID IS EQUIVALENT TO THE SPACING IN
%%RGRID...
if strcmp(typ,'tukeyC') || strcmp(typ,'tukeyIsoC') || strcmp(typ,'linearC') || strcmp(typ,'linearIsoC') || strcmp(typ,'quadrC') || strcmp(typ,'quadrIsoC');ck=(N+1)/2;else ck=ceil((N+1)/2);end
kGrid=generateGrid(N,gpu,2*1i*pi./sp,ck);
rGrid=generateGrid(N,gpu);for n=1:ND;rGrid{n}(:)=0;end

if strcmp(typ,'1stFinite')%First order continuous finite differences   
    for n=1:ND;x=bsxfun(@plus,x,kGrid{n});end
elseif strcmp(typ,'2ndFinite')%Second order continuous finite differences
    for n=1:ND;x=bsxfun(@plus,x,kGrid{n}.^2);end
elseif strcmp(typ,'1stFiniteDiscreteForward')%First order forward finite differences
    assert(all(~c),'1st Order Finite Discrete Difference filter not defined for DCT\n');
    for n=1:ND
        if N(n)>1;rGrid{n}(1)=1;rGrid{n}(end)=-1;rGrid{n}=-rGrid{n}/(sp(n));end
    end
elseif strcmp(typ,'1stFiniteDiscreteBackward')%First order forward finite differences
    assert(all(~c),'1st Order Finite Discrete Difference filter not defined for DCT\n');
    for n=1:ND
        if N(n)>1;rGrid{n}(1)=-1;rGrid{n}(2)=1;rGrid{n}=-rGrid{n}/(sp(n));end
    end
elseif strcmp(typ,'1stFiniteDiscreteCentered')%First order centered finite differences
    assert(all(~c),'1st Order Finite Discrete Difference filter not defined for DCT\n');
    for n=1:ND
        if N(n)>1;rGrid{n}(1)=0;rGrid{n}(2)=1;rGrid{n}(end)=-1;rGrid{n}=-rGrid{n}/(2*sp(n));end
    end
elseif strcmp(typ,'2ndFiniteDiscrete')%This corresponds, for instance, to 
    %the cosine profile in DC Ghiglia and LA Romero, "Robust 
    %two-dimensional weighted and unweighted phase unwrapping that uses
    %fast transforms and iterative methods", J Opt Soc Am A, 11(1):107-117, 
    %Jan 1994
    for n=1:ND
        if N(n)>2;rGrid{n}(1)=-2;rGrid{n}(2)=1;rGrid{n}(end)=1;rGrid{n}=rGrid{n}/(sp(n)^2);end
    end
elseif strcmp(typ,'FractionalFiniteDiscrete')
    for n=1:ND
        if N(n)>1;x=bsxfun(@plus,x,(2*abs(sin(imag(kGrid{n})/2))).^gibbsRing);end        
    end
elseif strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm') || strcmp(typ,'FractionalIso') || ... 
        strcmp(typ,'tukeyIso') || strcmp(typ,'tukeyIsoC') || strcmp(typ,'linearIso') || strcmp(typ,'linearIsoC') || ...
        strcmp(typ,'quadrIso') || strcmp(typ,'quadrIsoC') || strcmp(typ,'butter') || strcmp(typ,'gauss')
    if length(gibbsRing)==1
        for n=1:ND
            if N(n)>1;x=bsxfun(@plus,x,abs(kGrid{n}).^2);end
        end
    else%THIS ONLY WORKS FOR TUKEY
        for n=1:ND
            kGrid{n}=kGrid{n}*(max(1-gibbsRing)/(1-gibbsRing(n)));
            if N(n)>1;x=bsxfun(@plus,x,abs(kGrid{n}).^2);end
        end
        gibbsRing=min(gibbsRing);
    end
    x=sqrt(x);    
    if strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm') || strcmp(typ,'FractionalIso')    
        if strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm')
            x(x>pi)=pi;
            if strcmp(typ,'FractionalFiniteDiscreteIso');mult=2;else mult=1;end
            x=(mult*abs(sin(x/2)));
        end
        x=x.^gibbsRing;
    elseif strcmp(typ,'butter')
        x=1./(1+(x/pi).^(2*gibbsRing));       
    elseif strcmp(typ,'gauss')
        sigma=gibbsRing/(2*sqrt(2*log(2)));%For FWHM as a ratio of temporal window and temporal reconstruction resolution    
        x=exp(-(x.^2)*(2*sigma^2));        
    elseif strcmp(typ,'tukeyIso') || strcmp(typ,'tukeyIsoC') || strcmp(typ,'linearIso') || strcmp(typ,'linearIsoC') || ...
            strcmp(typ,'quadrIso') || strcmp(typ,'quadrIsoC') || strcmp(typ,'butter') || strcmp(typ,'gauss')
        kk=x;
        x(:)=1;        
        alpha=1-gibbsRing;
        if gibbsRing~=0
            if strcmp(typ,'linearIso') || strcmp(typ,'linearIsoC')
                fkk=1-((kk-pi*alpha)/((1-alpha)*pi));
            elseif strcmp(typ,'quadrIso') || strcmp(typ,'quadrIsoC')
                fkk=1-sqrt(((kk-pi*alpha)/((1-alpha)*pi)));
            else
                fkk=0.5*(1+cos(pi*((kk-pi*alpha)/((1-alpha)*pi))));
            end
            x(kk>=pi*alpha)=fkk(kk>=pi*alpha);
        end
        x(kk>=pi)=0;
    end
elseif strcmp(typ,'tukey') || strcmp(typ,'tukeyC') || strcmp(typ,'linear') || strcmp(typ,'linearC')
    N(end+1:2)=1;x=single(ones(N));
    if gpu>0;x=gpuArray(x);end
    gibbsRing(end+1:ND)=gibbsRing(end);
    for m=1:ND
        Naux=ones(1,ND);Naux(m)=N(m);spaux=ones(1,ND);spaux(m)=sp(m);        
        if strcmp(typ,'tukey');tukAux='tukeyIso';
        elseif strcmp(typ,'tukeyC');tukAux='tukeyIsoC';
        elseif strcmp(typ,'linear');tukAux='linearIso';
        elseif strcmp(typ,'linearC');tukAux='linearIsoC';
        elseif strcmp(typ,'quadr');tukAux='quadrIso';
        elseif strcmp(typ,'quadrC');tukAux='quadrIsoC';
        end
        x=bsxfun(@times,x,buildFilter(Naux,tukAux,spaux,gpu,gibbsRing(m)));
    end
    x=fftshift(x);
elseif strcmp(typ,'CubicBSpline')
    x=1;
    for n=1:ND
        x=bsxfun(@times,x,6*bsxfun(@rdivide,sinc(imag(kGrid{n}/(2*pi))).^4,4+2*cos(imag(kGrid{n}))));
    end
elseif strcmp(typ,'Spline')
    x=1;
    for n=1:ND
        %imag(kGrid{n}(:))'
        %pause
        x=bsxfun(@times,x,sinc(imag(kGrid{n}/(2*pi))).^(gibbsRing+1));
    end
else
    error('Unknown filter %s',typ);
end


if strcmp(typ,'1stFiniteDiscreteForward') || strcmp(typ,'1stFiniteDiscreteBackward') || strcmp(typ,'1stFiniteDiscreteCentered') || strcmp(typ,'2ndFiniteDiscrete')
    for n=1:ND
        if numel(rGrid{n})>1
            rGrid{n}=fftGPU(rGrid{n},n);
            if c(n);rGrid{n}=real(rGrid{n});end
            rGrid{n}=fftshift(rGrid{n},n);
            x=bsxfun(@plus,x,rGrid{n});
        end
    end
end
x=ifftshift(x);
if any(c)%We assume the filter is symmetric
    N=size(x);    
    assert(isreal(x) || any(~c),'The filter is not real, so it cannot be applied in the cosine domain');
    assert(all(mod(N(c==1),2)==0 | N(c==1)==1),'The filter does not have an even size in all dimensions as its size is%s',sprintf(' %d',N)); 
    v=cell(1,ND);
    for m=1:ND
        if c(m);v{m}=1:ceil(N(m)/2);else v{m}=1:N(m);end
    end
    x=dynInd(x,v,1:ND);
end

if (strcmp(typ,'FractionalFiniteDiscreteIso') || strcmp(typ,'FractionalFiniteDiscreteIsoNorm') || strcmp(typ,'FractionalIso')) && gibbsRing<0;x(1)=max(x(2:end));end%To prevent numerical instabilities
