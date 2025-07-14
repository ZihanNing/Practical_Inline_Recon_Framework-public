function [Gx,Gy,Gz] = calcGradWvfrm(N, steps, rampSteps)
    
    k=1:N;
    
    % calculate phy and theta
    hk = -1+2*(k-1)/(N-1);
    theta = acos(hk);

    phy = zeros(size(theta));
    for i = 2:N-1
        phy(i) = mod(phy(i-1) + (3.6/N^0.5) /((1-hk(i)^2)^0.5),2*pi);
    end

    % calculate x,y,z gradients
    Gx = sin(theta).*cos(phy); 
    Gy = sin(theta).*sin(phy);
    Gz = cos(theta);

    rampUpFactor = [rampSteps/2,ones(1,steps-2*size(rampSteps,2)),rampSteps(end:-1:1)/2]; 
    rampUpFactor = repmat(rampUpFactor,N,1);

    Gx = repmat(Gx',[1,steps]).*rampUpFactor;
    Gy = repmat(Gy',[1,steps]).*rampUpFactor;
    Gz = repmat(Gz',[1,steps]).*rampUpFactor;
end