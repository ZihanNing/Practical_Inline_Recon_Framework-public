function gdat = rollKernel(gdata,MTX)
effMtx = ceil(1.5*MTX);
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
DCF   = [1.0];
numThread = 1; % only have 1 data point
rokern = grid3_MAT(delta',k_not',DCF,effMtx,numThread);

gdata = squeeze(gdata(1,:,:,:) + 1j*gdata(2,:,:,:));
gdata = fftn(gdata);
gdata = fftshift(gdata,1);
gdata = fftshift(gdata,2);
gdata = fftshift(gdata,3);

% ROLLOFF
% change to complex, shift, then fft
rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
rokern = fftn(rokern);
rokern = fftshift(rokern,1);
rokern = fftshift(rokern,2);
rokern = fftshift(rokern,3);
rokern = abs(rokern);

%'   apply rolloff and crop'
gdata(rokern > 0) = gdata(rokern > 0) ./ rokern(rokern > 0);
xs = floor(effMtx/2 - effMtx/1.5/2+1);
if (mod(MTX,2)) % odd matrix size
    xe = floor(effMtx/2 + effMtx/1.5/2-1); 
else % even matrix size
    xe = floor(effMtx/2 + effMtx/1.5/2);
end
gdat = gdata(xs:xe,xs:xe,xs:xe);
%gdat = single(abs(gdata));
end