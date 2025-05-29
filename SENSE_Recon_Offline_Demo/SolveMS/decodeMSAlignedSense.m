function xo=decodeMSAlignedSense(xi,S,F,G)

%DECODEMSALIGNEDSENSE   Backward model of a multi-slice aligned SENSE 
%reconstruction for Siemens data
%   XO=DECODEMSALIGNEDSENSE(XI,{S},{F},{G})  
%   * XI is the Fourier data
%   * {S} is the coil array sensitivity map
%   * {F} contains the Fourier encoding
%   * {G} contains the SMS encoding
%   ** XO is the decoded data
%

if nargin<2;S=[];end
if nargin<3;F=[];end
NP=length(F);
if nargin<4;G=cell(1,NP);end

for n=1:NP
    if ~isempty(G{n})
        NG=size(G{n},3:4);
        xi{n}=ifold(xi{n},3,prod(NG),NG(1),[],[],2);
        xi{n}=xi{n}.*resPop(conj(G{n}),3:4,[],3);
    end   
    if ~isempty(F{n});xi{n}=aplGPU(F{n},xi{n},2);end
    if ~isempty(S);xi{n}=sum(xi{n}.*conj(dynInd(S,n,5)),4);end
    if n==1;xo=xi{n};else xo=xo+xi{n};end
end
