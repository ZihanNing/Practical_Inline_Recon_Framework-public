function [x,P]=ridgeDetectionExt(x,dims,up,P)

%RIDGEDETECTION   Fits a linear phase to the spectrum by detecting its peak
%(slope) and the phase of the peak (intercept)
%   X=RIDGEDETECTION(X,{DIMS})
%   * X is the array on which to detect the ridge
%   * {DIMS} are the dimensions across which to detect the ridge. It defaults
%   to all dimensions in the data
%   * UP is an upsampling factor of the spectrum for more accurate location of 
%   the peak. It defaults to 1
%   * X is the ridge information
%   * P are the ridge detection parameters (constant phase and linear
%   ramps, center of coordinates is given by the first element in the
%   array)
%

if nargin<3 || isempty(up);up=1;end
ND=numDims(x);N=size(x);
if nargin<2 || isempty(dims);dims=1:ND;end
if nargin >3 && ~isempty(P); di=1;else di=0;end

N(end+1:max(dims))=1;
ND=length(N);

gpu=isa(x,'gpuArray');

N(dims)=N(dims)*up;
nodims=1:ND;nodims(dims)=[];
Nnodims=N;Nnodims(dims)=1;

LD=length(dims);

%A PROBLEM ONLY PARTIALLY CONSIDERED IS RINGING INTRODUCED IN THE SPECTRUM
%BY NON COMPACT SIGNALS, A BIT OF WINDOWING MAY HELP, ALSO ONE COULD USE
%THE STATISTICS OF THE PHASE TO MODULATE BY A SPATIAL PROFILE DIFFERENT
%FROM THE SIGNAL MAGNITUDE, ALTHOUGH THAT COULD BE DONE OUTSIDE THIS
%FUNCTION
x=resampling(x,N,2);%Zero padded
for n=1:LD;x=fftGPU(x,dims(n));end%Not fftshifted, so origin at first element of array

x=resSub(x,dims); %dims must be contiguous dimensions
[x,ix]=max(x,[],dims(1));%YB: take largest component in k-space
x=reshape(x,Nnodims);
ix=ind2subV(N(dims),ix);

%PARAMETERS
if di==0%detection
    P=real(x);P(:)=0;
    repm=ones(1,max(numDims(P),dims(1)));repm(dims(1))=LD+1;%the offset and all the slopes
    repm(end+1:2)=1; %YB addition
    P=repmat(P,repm);
end

%%% Store all the slopes
if di==0%detection
    for n=1:LD;P=dynInd(P,1+n,dims(1),resPop(wrapToPi(2*pi*(dynInd(ix,n,2)-1)/N(dims(n))),1,Nnodims(nodims),nodims));end
else
    for n=1:LD
    ix = dynInd(ix,n,2, (round(...
        resPop((dynInd(P,1+n,dims(1))+2*pi*single(reshape(dynInd(P,1+n,dims(1))<0, size(dynInd(P,1+n,dims(1))) )))*N(dims(n))/2/pi+1,nodims,[],1 ) ...
        ) ));
    end%Fill it
    x=ones(size(x));
end

%%% Calculate phase of the harmonics corresponding to peak identification
%%% to extract the phase terms
F=buildStandardDFTM(N(dims),0,gpu);
for n=1:LD
    F{n}=F{n}(ix(:,n),:);
    F{n}=shiftdim(F{n},-ND);
    F{n}=resPop(F{n},ND+1,Nnodims(nodims),nodims);
    F{n}=resPop(F{n},ND+2,N(dims(n)),dims(n));
    
    NF=size(F{n});
    NF(dims(n))=NF(dims(n))/up;
    F{n}=resampling(F{n},NF,2);
end
x=sign(x);
for n=1:LD;x=bsxfun(@times,x,conj(F{n}));end 
if di==1%apply
    %Need to take out the phase of the first element of the Fourier harmonics - needed as upsampling might result in first element not zero phase and this
    %element influences the phase in dynInd(P,1,dims(1))
    x = bsxfun(@times, x, conj(dynInd(x, ones(1,LD), dims) ));
    x = bsxfun(@times, x, exp(+1i*(dynInd(P, 1, dims(1))) ) );
else
    P=dynInd(P,1,dims(1),dynInd(angle(x),ones(1,LD),dims));%Take first element for the offset
end