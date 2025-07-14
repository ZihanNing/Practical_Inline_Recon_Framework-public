function [kx,ky,kz] = calcKspaceRef(Gx,Gy,Gz,NUFFT)

    kxStart = cumsum(Gx(:,end/2:end),2);
    kyStart = cumsum(Gy(:,end/2:end),2);
    kzStart = cumsum(Gz(:,end/2:end),2);

    
    kx = repmat(kxStart(:,end),[1,size(Gx,2)])-cumsum(Gx,2);
    ky = repmat(kyStart(:,end),[1,size(Gx,2)])-cumsum(Gy,2);
    kz = repmat(kzStart(:,end),[1,size(Gx,2)])-cumsum(Gz,2);

   
    if ~NUFFT 
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

