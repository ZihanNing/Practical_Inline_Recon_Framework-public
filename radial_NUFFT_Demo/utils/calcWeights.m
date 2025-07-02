function [w, gauss] = calcWeights(kx,ky,kz) 
%% calculate weighting function

    [N,steps] = size(kx);
    k_space = [kx(:),ky(:),kz(:)]';
    weights = k_space.*k_space;
    weights = sum(weights,1);
    weights = reshape(weights,N,steps);
    wmax = weights(end);
    
    %% Weight the kspace data
    
    hannf = hann(steps*2);
    hannf = repmat(hannf(steps+1:end)',[N,1]);
    
    w = weights.*hannf;