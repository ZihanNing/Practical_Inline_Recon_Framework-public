function [Gx,Gy,Gz] = calcGradWvfrmRef(N, steps, rampSteps)
    
    k=1:N;
    
    % calculate phy and theta
    hk = -1+2*(k-1)/(N-1);
    theta = acos(hk);

    phy = zeros(size(theta));
    for i = 2:N-1
        phy(i) = mod(phy(i-1) + (3.6/N^0.5) /((1-hk(i)^2)^0.5),2*pi);
    end

    % calculate x,y,z gradients (TODO: what is the gradient strength?)
    Gx = sin(theta).*cos(phy); 
    Gy = sin(theta).*sin(phy);
    Gz = cos(theta);

    %YB: changed this... cause ramp-up kinda weird according to POET
%     rampUpFactor = [linspace(0,1,rampSteps), ones(1,steps-rampSteps)];
%     rampUpFactor = repmat(rampUpFactor,N,1);
    
    %YB: new calculation TODO: generalise this!!!
    rs = [linspace(0.5,rampSteps-0.5,rampSteps)]; rs = [rs(:)';rs(:)']; rs = [rs(:)']/rampSteps; % YB: complicated ramp-up....  
    rs = [0,rs(2:end)]; %YB: changed after ADC adjustment!!!
    rampUpFactor = [rs, ones(1,steps-rampSteps*4),rs(end:-1:1)]; 
    rampUpFactor = repmat(rampUpFactor,N,1);

    Gx = repmat(Gx',[1,steps]).*rampUpFactor*-1;
    Gy = repmat(Gy',[1,steps]).*rampUpFactor*-1;
    Gz = repmat(Gz',[1,steps]).*rampUpFactor*-1;

    
    
    
end