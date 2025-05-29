function x=filtering(x,H,c,f)

% FILTERING filters an array in Fourier domain
%   X=FILTERING(X,H,{C},{F},{F},{FH})
%   * X is the array to be filtered
%   * H is the filter to be applied
%   * {C} indicates if the filter is to be applied in the Cosine domain 
%   (for Neumann instead of periodic boundary conditions), defaults to 0, 
%   i.e., the filter is applied in the Fourier domain
%   * {F} indicates if the data is in the Fourier domain already, defaults 
%   to 0
%   ** X is the filtered image
%

if nargin<3 || isempty(c);c=0;end
if nargin<4 || isempty(f);f=0;end

gpu=isa(x,'gpuArray');

NH=size(H);
nDimsH=ndims(H);
comp=~isreal(x) || ~isreal(H);

if length(c)==1;c=c*ones(1,nDimsH);end

if sum(NH~=1)==1%Filtering on single direction
    if ~f
        d=find(NH~=1);
        if c(d);F=dctmtx(NH(d));else F=dftmtx(NH(d))/sqrt(NH(d));end
        if (gpu && ~isaUnderlying(x,'double')) || (~gpu && ~isa(x,'double'));F=single(F);end
        if gpu;F=gpuArray(F);end
        F=F'*bsxfun(@times,H(:),F);
        if ~comp;F=real(F);end
        x=fftGPU(x,d,F);%THIS IS ACTUALLY THE APPLICATION OF A GENERIC SQUARE MATRIX           
    else
        x=bsxfun(@times,x,H);
    end
else
    for m=1:nDimsH
        if NH(m)~=1 && ~f
            if ~c(m);x=fftGPU(x,m);
            else x=fctGPU(x,m);
            end
        end
    end
   
    x=bsxfun(@times,x,H);

    for m=1:nDimsH
        if NH(m)~=1 && ~f
            if ~c(m);x=ifftGPU(x,m);
            else x=ifctGPU(x,m);
            end
        end
    end
end
if ~comp;x=real(x);end%%%THIS IS PROBLEMATIC, IT DOES NOT HOLD IF H IS NOT EVEN!!
