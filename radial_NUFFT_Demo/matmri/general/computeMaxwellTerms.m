function maxOut = computeMaxwellTerms(orders,coilParams,Gx,Gy,Gz,B0)
%[phi, phiConcTerms] = computeConcPhase(obj,linPhase,f0,dt,nonLinSphHarm)
% Compute phase from maxwell fields (i.e., concomitant gradients). 
% Everything should be in SI units. 

if isempty(coilParams)
    warning('No coil parameters provided - assuming a symmetric coil');
    coilParams.alph = 0.5;
    coilParams.z0x = 0;
    coilParams.z0y = 0;
    coilParams.x0 = 0;
    coilParams.y0 = 0;
end

% Initialize
didPrep = 0;

% Compute zeroth order terms
if any(orders==0)
    [Gx2, Gy2, Gz2, GzGx, GzGy] = doPrepCalcs(Gx,Gy,Gz);
    didPrep = 1;
    maxOut.zeroth = coilParams.z0x^2 * Gx2 / (2*B0) + ...
        coilParams.z0y^2 * Gy2 / (2*B0) + ...
        coilParams.alph^2     * coilParams.x0^2 * Gz2 / (2*B0) + ...
        (1-coilParams.alph)^2 * coilParams.y0^2 * Gz2 / (2*B0) + ...
        -coilParams.alph      * coilParams.x0*coilParams.z0x * GzGx / B0 + ...
        -(1-coilParams.alph)  * coilParams.y0*coilParams.z0y * GzGy / B0;
else
    maxOut.zeroth = [];
end

% Compute first order terms
if any(orders==1)
    if ~didPrep
        [Gx2, Gy2, Gz2, GzGx, GzGy] = doPrepCalcs(Gx,Gy,Gz);
        didPrep = 1;
    end
    maxOut.x = coilParams.alph^2 * coilParams.x0 * Gz2 / B0 + ...
        -coilParams.alph * coilParams.z0x * GzGx / B0;
    maxOut.y = (1-coilParams.alph)^2 * coilParams.y0 * Gz2 / B0 + ...
        -(1-coilParams.alph) * coilParams.z0y * GzGy / B0;
    maxOut.z = coilParams.z0x * Gx2 / B0 + ...
        coilParams.z0y * Gy2 / B0 + ...
        -coilParams.alph     * coilParams.x0 * GzGx / B0 + ...
        -(1-coilParams.alph) * coilParams.y0 * GzGy / B0;
else
    maxOut.x = [];
    maxOut.y = [];
    maxOut.z = [];
end

% Compute 2nd order terms
if any(orders==2)
    if ~didPrep
        [Gx2, Gy2, Gz2, GzGx, GzGy] = doPrepCalcs(Gx,Gy,Gz);
        didPrep = 1;
    end
    % z^2
    maxOut.z2 = (Gx2 + Gy2)/(2*B0);
    % x^2
    maxOut.x2 = coilParams.alph^2 * Gz2 / (2*B0);
    % y^2
    maxOut.y2 = (1-coilParams.alph)^2 * Gz2 / (2*B0);
    % zx
    maxOut.zx = -coilParams.alph * GzGx / B0;
    % zx
    maxOut.zy = -(1-coilParams.alph) * GzGy / B0;
else
    maxOut.z2 = [];
    maxOut.x2 = [];
    maxOut.y2 = [];
    maxOut.zx = [];
    maxOut.zy = [];
end

end

function [Gx2, Gy2, Gz2, GzGx, GzGy] = doPrepCalcs(Gx,Gy,Gz)
    Gx2 = Gx.^2;
    Gy2 = Gy.^2;
    Gz2 = Gz.^2;
    GzGx = Gz.*Gx;
    GzGy = Gz.*Gy;
end



