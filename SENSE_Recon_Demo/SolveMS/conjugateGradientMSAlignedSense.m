function [x,y,n]=conjugateGradientMSAlignedSense(y,S,F,G,M,x,nX,toler)

%CONJUGATEGRADIENTMSALIGNEDSENSE   Conjugate gradient of a multi-slice 
%aligned SENSE reconstruction for Siemens data
%   X=CONJUGATEGRADIENTMSALIGNEDSENSE(Y,{S},{F},{G},{X},{NX},{TOLER}) 
%   * Y is the data (k-space)
%   * {S} is the coil array sensitivity map
%   * {F} contains the Fourier encoding
%   * {G} contains the SMSr encoding
%   * {M} is a mask
%   * {X} is a starting reconstruction
%   * {NX} is the number of iterations of the CG algorithm
%   * {TOLER} is the maximum update for convergence. It defaults to 0
%   ** X the reconstructed image
%   ** Y are the predictions
%   ** N is the number of effective iterations
%

if nargin<2;S=[];end
if nargin<3;F=[];end
if nargin<4;G=[];end
if nargin<5;M=[];end
if nargin<6;x=[];end
if nargin<7 || isempty(nX);nX=300;end
if nargin<8 || isempty(toler);toler=1e-5;end


n=0;%Effective iterations
ND=12;%Maximum dimensionality

isSoftMask=any(~ismember(M(:),[0 1]));
%PRECONDITIONER
if isSoftMask;P=preconditionMSAlignedSense(S,M);else P=preconditionMSAlignedSense(S);end

NP=length(y);
%DECODE
r=decodeMSAlignedSense(y,S,F,G);n=n+1;
if ~isSoftMask && ~isempty(M);r=r.*M;end
for m=1:NP;y{m}=gather(y{m});end

NX=size(r);
if isempty(x);x=zeros(NX,'like',r);
else r=r-systemMSAlignedSense(x);n=n+2;        
end
z=P.*r;

p=z;
rsold=sum(conj(z).*r,1:ND);
l=true;
if sqrt(min(abs(rsold(:))))<1e-9;l=false;end

%ITERATIONS
while l
    %SYSTEM MATRIX
    Ap=systemMSAlignedSense(p);n=n+2;
    
    %UPDATES
    al=conj(rsold)./sum(conj(p).*Ap,1:ND);
    xup=al.*p;
    x=x+xup; 
    if n>3
        xup=max(abs(xup(:)))/(max(abs(x(:)))+1e-6);
        fprintf('Inner iteration: %d / Maximum update: %0.2g\n',n,xup);
    else
        xup=inf;
    end
    if xup<toler || n>=nX;break;end
    r=r-al.*Ap;
    z=P.*r;

    rsnew=sum(conj(z).*r,1:ND);
    if sqrt(min(abs(rsnew(:))))<1e-9;break;end
    be=rsnew./rsold;
    p=z+be.*p;
    rsold=rsnew;
end

%LOSS
if nargout>=2;y=encodeMSAlignedSense(x,S,F);n=n+1;end

function x=systemMSAlignedSense(x)
    if ~isempty(M);x0=x;end
    x=encodeMSAlignedSense(x,S,F,G);
    x=decodeMSAlignedSense(x,S,F,G);
    if ~isempty(M) 
        if isSoftMask;x=x+x0./M.^2;else x=x.*M;end
    end
end

end