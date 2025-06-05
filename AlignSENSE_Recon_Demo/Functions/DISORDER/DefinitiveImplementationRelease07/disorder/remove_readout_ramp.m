

function [x, shift] = remove_readout_ramp(x, dim, shift,  do_shift)

if nargin < 3 || isempty(shift); shift = [];end
if nargin < 4 || isempty(do_shift); do_shift = 1;end

if do_shift

    x = fftn(fftshift(x));
    %x = fftn(x);
    
    if isempty(shift)
         cent = floor(size(x)/2) + 1;
         cent_line = cent; cent_line(dim)= [];
         dims = 1:3; dims(dim) = [];
    
        readout_line = dynInd(x, cent_line, dims);
        [~,idx] = max(abs(readout_line));
        shift = idx-cent(dim);
    end
            
    x = circshift(x, -shift, dim);


    x = ifftshift(ifftn(x));
    %x = ifftn((x));
    
else
    shift = 0;
end







