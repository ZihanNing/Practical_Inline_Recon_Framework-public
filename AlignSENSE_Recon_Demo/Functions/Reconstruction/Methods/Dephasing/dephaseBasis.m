
function [x, p]=dephaseBasis(B,c,NX, TE, returnField)

%DEPHASEBASIS   Computes the dephasing for a given set of basis functions
%and coefficients
%   X=DEPHASEBASIS(B, C, NX, TE)
%   * B = set of basis functions ordered as [prod(NX) Ncoef]
%   * c are the coefficients for corresponing motion states/basis functions
%   * NX is the spatial dimensions
%   * {TE} is the echo time (in seconds)
%   ** x is a complex image with phase structure
%   ** p is the phase not yet applied to the image and wrapped in the exp() operator

if nargin<4 || isempty(TE); TE=1;end
if nargin<5 || isempty(returnField); returnField=0;end

gpu = isa(B,'gpuArray');
if gpu && ~isa(c,'gpuArray'); c = gpuArray(c);end

nT = size(c,5);
nC = size(c,6);
c = permute(c, [6 5 1:4]);
NX = NX(1:3);%Only spatial dimensions needed

BlSzC = 5; %For coefficients

x = zeros(prod(NX) , nT, 'like', B); 
for ii= 1:BlSzC:nC
    vI=ii:min(ii+BlSzC-1,nC);
    x = x + dynInd(B, vI,2) * dynInd(c,vI,1 );
end

x=reshape(x , [NX 1 nT] ); %Motion states should be in 5th dimension / x in Hz
if returnField; return; end %In Hz

if nargout>1; p = 2*pi*TE * x;end %phase in radians
x=exp(+1i *2*pi*TE * x);% wrapped phase

