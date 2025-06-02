function T=upsampleTransform(T,NSh,NShPrev)

%UPSAMPLETRANSFORM upsamples a given set of motion states to a finer
%temporal resolution
%   T=UPSAMPLETRANSFORM(T)
%   * T is the transform at the coarser resolution
%   * NSH is the new number of shots
%   * NSHSPREV is the previous number of shots
%   * T is the transform at the finer resolution
%

ND=numDims(T);

%TRANSFORM IN THE UPSAMPLED GRID
upFactor=NSh/NShPrev;
NS=size(T,ND-1);
T=permute(T,[ND-1:ND 1:ND-2]);
T=mirroring(T,1,1);
T=resampling(T,NS*2*ceil(upFactor));

%UPSAMPLING FILTER
kernUp=single(zeros([NS*2*ceil(upFactor) 1]));
remFactor=upFactor-floor(upFactor);
if remFactor==0          
    kernUp(1:upFactor)=1/upFactor;
else
    kernUp(2:floor(upFactor)+1)=1/upFactor;
    kernUp([1 ceil(upFactor)+1])=(1-floor(upFactor)/upFactor)/2;
end
kernUp=fft(kernUp);

%UPSAMPLE
T=real(filtering(T,kernUp));
T=real(resampling(T,NSh*2));
T=mirroring(T,1,0);
T=permute(T,[3:ND 1:2]);