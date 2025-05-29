function xo=encodeMSAlignedSense(xi,S,F,G,y)

%ENCODEMSALIGNEDSENSE   Forward model of a multi-slice aligned SENSE 
%reconstruction for Siemens data
%   XO=ENCODEMSALIGNEDSENSE(XI,{S},{F},{G},{Y})  
%   * X is the image
%   * {S} is the coil array sensitivity map
%   * {F} contains the Fourier encoding
%   * {G} contains the SMS encoding
%   * {Y} contains the samples, for computing the residuals
%   ** XO is the encoded data
%

if nargin<2;S=[];end
if nargin<3;F=[];end
NP=length(F);
if nargin<4;G=cell(1,NP);end
if nargin<5;y=cell(1,NP);end

xo=cell(1,NP);
for n=1:NP
    if ~isempty(S);xo{n}=xi.*dynInd(S,n,5);else xo{n}=xi;end
    if ~isempty(F{n});xo{n}=aplGPU(F{n}',xo{n},2);end
    if ~isempty(G{n})
        NG=size(G{n},3:4);
        xo{n}=xo{n}.*resPop(G{n},3:4,[],3);
        xo{n}=fold(xo{n},3,prod(NG),NG(1),[],[],2);
    end    
    if ~isempty(y{n});xo{n}=xo{n}-y{n};end
end
