

function [chi] = guessChi(x, MS, chi)
%x image
%MS, resolution
%chi, initial chi estimate
if nargin<2 || isempty(MS);MS=ones([1 3]);end
if nargin<3 || isempty(chi); chi = zeros(size(x), 'like', real(x));else chi = chi - multDimMea(chi);end
assert(all(isreal(chi(:))),'guessChi:: chi must be real.');

%%% Create mask for weighted interpolation
parS=[];parS.Otsu = [0 1];
parS.conComp=0;
parS.nDilate=0;%In mm
M = refineMask(x, parS,MS);

chi(M==1) =  0;
chi(M==0) = -9 ;
