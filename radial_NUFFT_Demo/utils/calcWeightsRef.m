function w = calcWeightsRef(kx,ky,kz) 
%% calculate weighting function

    [N,steps] = size(kx);
    k_space = [kx(:),ky(:),kz(:)]';
    weights = k_space.*k_space;
    weights = sum(weights,1);
    weights = reshape(weights,N,steps);
    
    %% Weight the kspace data
    %sig = 0.3;
    %l =steps/2;
   % d = [1:steps/2];
    %gauss = 1/((l*sig)*(2*pi)^0.5) * exp(-0.5*((d)./(l*sig)).^2);
    %gauss = [gauss(end:-1:1),gauss];
    %gauss = repmat(gauss,[N,1]);
    
   % w = weights.*gauss;
    
    hannf = hann(steps);
    hannf = repmat(hannf',[N,1]);
    
    w = weights.*hannf;