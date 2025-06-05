function H=buildHalfScanFilter(N,pFF,di,filt,gpu)

% BUILDHALFSCANFILTER builds a ramped partial Fourier filter for a
% given dimension and partial fourier factor
%   H=BUILDHALFSCANFILTER(X,PFF)
%   * N are the dimensions of the grid
%   * PF is the Partial Fourier factor
%   * {DI} is the direction of the filter, 1 to apply Partial Fourier, 0 to
%   reverse, it defaults to 1
%   * {FILT} is the type of filter for covariance calculation. Possibilities are
%   'ramp' and 'zefi'. It defaults to ramp
%   * {GPU} determines whether to use gpu computations 
%   ** H is the generated filter
%

if nargin<3 || isempty(di);di=1;end
if nargin<4 || isempty(filt);filt='ramp';end
if nargin<5 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
    
kGrid=generateGrid(N,gpu,N,ceil((N+1)/2));
kGrid=kGrid{1};

cutOffInd=round(N*(1-pFF)+1);        
cutOff=gather(abs(kGrid(cutOffInd)-0.5));        
indMA{1}=(kGrid(:)<=-cutOff);indMA{2}=(kGrid(:)>-cutOff & kGrid(:)<=cutOff);indMA{3}=(kGrid(:)>cutOff);
               
kyRM=cell(1,3);for n=1:3;kyRM{n}=gather(kGrid(indMA{n}));end
H=zeros([N 1],'single');
invert=1e-3;
invert2=1e-6;if invert>=invert2;invert2=0;end
if strcmp(filt,'ramp');H(indMA{1})=invert;H(indMA{2})=1+(1-invert)*kyRM{2}/cutOff;H(indMA{3})=2-invert;
else H(indMA{1})=invert;H(indMA{2})=1;H(indMA{3})=1;
end
if ~di
    H=1./(H+invert2);
    H(indMA{1})=H(find(indMA{2},1,'first'));
    X=flip(buildFilter(2*N,'tukey',1,gpu,0.05,1));
    X=circshift(X,cutOffInd-1);
    X(indMA{1})=0;
    H=H.*X;
end
if gpu;H=gpuArray(H);end
H=ifftshift(H);


