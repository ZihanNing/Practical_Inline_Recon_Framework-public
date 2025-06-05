function x=buildSliceProfile(NSl,FWHM,typ,ups,gpu)

%BUILDSLICEPROFILE generates an excitation profile
%   X=BUILDSLICEPROFILE({NSL},{FWHM},{TYP},{UPS},{GPU})
%   * {NSL} is the number of slices
%   * {FWHM} allows to change the resolution, by decreasing it the 
%   resolution should be smaller, defaults to 1
%   * {TYP} is the type of profile, defaults to 'Gauss', other options are 
%   'Sinc' and 'Hard'
%   * {UPS} is the upsampling factor, defaults to 1
%   * {GPU} is a flag for gpu array generation
%   ** X is the excitation profile
%

if nargin<1 || isempty(NSl);NSl=1;end
if nargin<2 || isempty(FWHM);FWHM=1;end
if nargin<3 || isempty(typ);typ='Gauss';end
if nargin<4 || isempty(ups);ups=1;end
if nargin<5 || isempty(gpu);gpu=useGPU;end

N=NSl*ups;

rGrid=generateGrid(N,gpu,2*pi);
rGrid=rGrid{1};
if strcmp(typ,'Gauss') || strcmp(typ,'Sinc')
    sigma=2*pi*(FWHM/ups)/(2*sqrt(2*log(2)));        
    x=exp(-(rGrid.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));

    x=ifftshift(x);
    x=x/x(1);
    if strcmp(typ,'Sinc');x(x>=0.5)=1;x(x<0.5)=0;end
    x=ifftGPU(x,1);
    x=x/x(1);
    x=abs(x);
elseif strcmp(typ,'Hard')
    x=rGrid;x(:)=0;
    x(1:round(ups/FWHM))=1;
    x=shifting(x,(ups-1)/2);
end

