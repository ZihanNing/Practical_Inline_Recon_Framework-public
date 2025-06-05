function x=margosianRegrid(x,indM,di,no)

% MARGOSIANREGRID resamples a given array using the FFT
%   X=MARGOSIANREGRID(X,INDM,{DI})
%   * X is the array to be regridded
%   * INDM is the information for regridding as provided by the
%   margosianFilter function
%   * {DI} indicates the direction for regridding, 1 to extract the sampled
%   asymmetric spectrum (default) / 0 to symmetrize the spectrum again
%   * {NO} indicates whether to perform noise normalization
%   * X is the regridded array
%

if ~exist('di','var') || isempty(di);di=1;end
if ~exist('no','var') || isempty(no);no=0;end

N=size(x);
ND=size(indM,1);
for m=1:ND        
    if ~isempty(indM{m}) && sum(indM{m}{1})~=0
        x=fftGPU(x,m)/(N(m)^(no/2));       
        if di
            x=dynInd(x,~indM{m}{1},m);           
        else
            Nres=size(x);Nres(m)=length(indM{m}{1});
            xRes=single(zeros(Nres),'like',x);
            x=dynInd(xRes,~indM{m}{1},m,x);
            xRes=[];          
        end
        N=size(x);
        x=ifftGPU(x,m)*(N(m)^(no/2));
    end
end
