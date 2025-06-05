function [ output ] = gradient_calc( m, data_mode, di )
% Calculates the gradient of a 3D data cube in the frequency domain, and outputs the gradient 
% magnitude or vector field
%
% Parameters:
% m: input image. Typically the magnitude or R2* data for regularization local weighting in QSM.
% data_mode: 0 to return the vector field. 1 for the L1 norm, and 2 for the L2 norm.
%
% Output:
% output: gradient image.
%
%  Last Modified by Carlos Milovic, 07.07.2020
    
    if nargin <2  || isempty(data_mode); data_mode = 0; end
    if nargin <3  || isempty(di); di = 1; end
    N = size(m);

    % Gradient operators in the frequency domain
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    if di == 0
        E1 = conj(E1);
        E2 = conj(E2);
        E3 = conj(E3);
    end
    % Convolution
    Gm1 = real(ifftn(E1.*fftn(m)));
    Gm2 = real(ifftn(E2.*fftn(m)));
    Gm3 = real(ifftn(E3.*fftn(m)));
    
    % Generate output
    if data_mode == 0
        output = (Gm1); 
        output(:,:,:,2) = (Gm2);
        output(:,:,:,3) = (Gm3);
    elseif data_mode == 1
        output = abs( Gm1 ) + abs( Gm2 ) + abs( Gm3 );
    elseif data_mode == 2
        output = sqrt( Gm1.^2 + Gm2.^2 +Gm3.^2 );
    end
    
    
end

