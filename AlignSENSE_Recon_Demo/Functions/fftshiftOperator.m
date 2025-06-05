

function [x] = fftshiftOperator(x, num, invert, dim)

if num==1
    if ~invert%%% Define this
        x = ifftshiftGPU(x,dim);
    else
        x = fftshiftGPU(x,dim);
    end
    
elseif num==2
    if ~invert%%% Define this
        x = fftshiftGPU(x,dim);
    else
        x = ifftshiftGPU(x,dim);
    end    
end

%The TWIX original one was fftshift for num==1
%and ifftshift for num==2

%I think it should be ifftshift (num=1) and fftshift (num=2)

end

