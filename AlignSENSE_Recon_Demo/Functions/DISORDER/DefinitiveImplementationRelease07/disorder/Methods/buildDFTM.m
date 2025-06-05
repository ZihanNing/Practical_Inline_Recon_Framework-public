function [DFTM,DFTMH]=buildDFTM(rGrid,kGrid,drGrid,dkGrid,factRegR,factRegK)

%BUILDDFTM   Builds a discrete Fourier transform matrix (either
%unidimensional, for vector inputs, or multidimensional, for cell inputs)
%   DFTM=BUILDDFTM(RGRID,KGRID,{DRGRID},{DKGRID},{FACTREGR},{FACTREGK})
%   * RGRID is a matrix with columns giving spatial grids and rows different
%   coordinates
%   * KGRID is a spectral grid with columns giving spectral grids and rows
%   different coordinates
%   * {DRGRID} is a matrix that indicates the density of the RGRID
%   * {DKGRID} is a matrix that indicates the density of the KGRID
%   * {FACTREGR} is a factor for regularizing the density compensation in
%   space
%   * {FACTREGK} is a factor for regularizing the density compensation in
%   the spectrum
%   * DFTM is the Fourier transform matrix for the corresponding grids
%   * DFTMH is the inverse Fourier transform matrix for the corresponding 
%   grids
%

assert(ismatrix(rGrid),'Only matricial arguments are accepted and rGrid is of %d dimensions',ndims(rGrid));
assert(ismatrix(kGrid),'Only matricial arguments are accepted and kGrid is of %d dimensions',ndims(kGrid));
assert(size(rGrid,2)==size(kGrid,2),'Number of columns not matched for rGrid (%d) and kGrid (%d)',size(rGrid,2),size(kGrid,2));

DFTM=kron(-2*pi*1i*kGrid(:,1),permute(rGrid(:,1),[2 1]));
for m=2:size(rGrid,2)
    DFTM=DFTM+kron(-2*pi*1i*kGrid(:,m),permute(rGrid(:,m),[2 1]));
end
DFTM=exp(DFTM)/sqrt(size(rGrid,1));%/sqrt(size(rGrid,1));
DFTMH=(DFTM');

if exist('dkGrid','var')
    assert(all(size(dkGrid)==size(kGrid)),'The sizes of spectral locations and density compensation grids are not the same');
    if ~exist('factRegK','var');factRegK=0.5;end
    DFTMH=densityCompensation(DFTMH,dkGrid,factRegK);
end

if exist('rkGrid','var')
    assert(all(size(rkGrid)==size(rGrid)),'The sizes of spatial locations and density compensation grids are not the same');
    if ~exist('factRegR','var');factRegR=0.5;end
    DFTM=densityCompensation(DFTM,drGrid,factRegR);
end

function x=densityCompensation(x,d,factReg)
    d=prod(abs(d),2);
    dTh=mean(d)*factReg;
    d=max(dTh,d);
    x=bsxfun(@times,x,1./d');
end

end
