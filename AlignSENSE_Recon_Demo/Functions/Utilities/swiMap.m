function xou=swiMap(xin)

%SWIMAP   Computes a SWI map
%   XOU=SWIMAP(XIN)
%   * XIN is the input data
%   * XOU is the output data
%

gpu=isa(xin,'gpuArray');
ND=max(numDims(xin),3);
NX=size(xin);NX(end+1:3)=1;

%HIGH PASS FILTER
HY=buildFilter(NX(1:3),'tukeyIso',0.8,gpu,1);
K1=imag(buildFilter(NX(1),'1stFinite',[],gpu));
K2=imag(buildFilter([1 NX(2)],'1stFinite',[],gpu));
K3=imag(buildFilter([1 1 NX(3)],'1stFinite',[],gpu));
KK=bsxfun(@plus,bsxfun(@plus,K1.^2,K2.^2),K3.^2);


KKZ=K1.^2;
HG=1/3-bsxfun(@rdivide,KKZ,KK+1e-9);
a=0.3;%0.25;
HGinv=HG;
HGinv(abs(HG)<a)=sign(HG(abs(HG)<a)).*(abs(HG(abs(HG)<a)/a).^2)/a;
HGinv(abs(HG)>=a)=1./HG(abs(HG)>=a);
HGinv(1)=0;

%RESHAPING TO ITERATE
[xin,NXS]=resPop(xin,4:ND,[],4);NXS(end+1:4)=1;
NI=3;
xou=xin;xou=repmat(xou,[ones(1,ND) NI]);

%ACTUAL COMPUTATION
for s=1:NXS(4)
    xs=dynInd(xin,s,4);
    
    %LOW RES IN Z DIRECTION
    xsp=filtering(xs,buildFilter(NX(1),'tukeyIso',0.5,gpu,1));
    
    %UNWRAPPING
    xsp=abs(xsp).*exp(1i*CNCGUnwrapping(xsp));
    
    %HIGH PASS FILTER THE PHASE
    xb=filtering(xsp,HY);    
    phi=angle(xsp.*conj(xb));    
    
    %DECONVOLVE SUSCEPTIBILITY
    phi=filtering(phi,HGinv);        
    phiw=max(min(-phi/pi+1,1),0).^4;
    xou=dynInd(xou,[s 2],[4 ND+1],abs(xs).*phiw.*exp(1i*phi));
    xou=dynInd(xou,[s 3],[4 ND+1],cdfFilt(abs(xs),'min',[35 1 1]));
end

%RESHAPING BACK
xou=reshape(xou,[NX NI]);
