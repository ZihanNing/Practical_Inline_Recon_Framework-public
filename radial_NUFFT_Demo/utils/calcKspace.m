function [kx,ky,kz] = calcKspace(Gx,Gy,Gz,NUFFT)

    kx = cumsum(Gx,2);
    ky = cumsum(Gy,2);
    kz = cumsum(Gz,2);
    
    if ~NUFFT %scale from -0.5 to 0.5 
        kx = (kx-min(kx(:)))/(max(kx(:)) - min(kx(:)));
        kx = kx-0.5;
        ky = (ky-min(ky(:)))/(max(ky(:)) - min(ky(:)));
        ky = ky-0.5;
        kz = (kz-min(kz(:)))/(max(kz(:)) - min(kz(:)));
        kz = kz-0.5;
    else
        kx = (kx-min(kx(:)))/(max(kx(:)) - min(kx(:)));
        kx = kx*2*pi - pi;
        ky = (ky-min(ky(:)))/(max(ky(:)) - min(ky(:)));
        ky = ky*2*pi - pi;
        kz = (kz-min(kz(:)))/(max(kz(:)) - min(kz(:)));
        kz = kz*2*pi - pi;
    end
end

