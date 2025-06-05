function data = kspace2image_gpu(data,order)
%%% This function is to transfer kspace to image by fft with GPU
%%%
%%% INPUT:
%%% data: kspace [Col, Cha, Lin, Par...]
%%% order: 
%%%   'ori' - first fftshift then ifftshift
%%%   'revert' - firt ifftshit then fftshit
%%% 
%%% OUTPUT:
%%% data: image domain
%%%
%%% by Zihan @ king's
%%% March-2025

if nargin < 2 || isempty(order); order='ori';end 

nDims = ndims(data);
permDim = [1 3 4 2, 5:nDims];
data = permute(data,permDim); % [Col Lin Par Cha...]

switch order
    case 'ori' % first fftshift then ifftshift (default)
        for l=1:3
            data = fftshiftGPU(data,l); 
            data = fftGPU(data,l);%*NY(l);  %YB: was ifft originally
            data = ifftshiftGPU(data,l);
        end
    case 'revert'
         for l=1:3
            data = ifftshiftGPU(data,l); 
            data = fftGPU(data,l);%*NY(l);  %YB: was ifft originally
            data = fftshiftGPU(data,l);
         end
end

fprintf('Transferred data from kspace to image domain. \n');