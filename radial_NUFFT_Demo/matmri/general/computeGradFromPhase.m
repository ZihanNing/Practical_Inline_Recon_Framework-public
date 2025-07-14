function grad = computeGradFromPhase(phs,dt,gamma)
% Computes MRI gradients based on phase acrual = 2*pi*k [rad/m], where
% everything is in SI units. Works along the first dimension.

% Find derivative of phase using Gaussian diff kernel
Nwd = 6; % should be even to end up with a centered kernel
gwindiff = diff(gausswin(Nwd)); 
norm = -floor(length(gwindiff)/2):floor(length(gwindiff)/2);
norm = abs(sum(norm(:).*gwindiff(:)));
gwindiff = gwindiff/norm;
dkdt = convn(phs,gwindiff(:),'same'); 

% Compute gradient
grad = dkdt/(2*pi*gamma*dt);

end

