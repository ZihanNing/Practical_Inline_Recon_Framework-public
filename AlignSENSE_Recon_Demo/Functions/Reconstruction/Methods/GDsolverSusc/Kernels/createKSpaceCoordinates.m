    

function [kx, ky, kz] = createKSpaceCoordinates(N,MS,MT)


    %%% LINEAR TERMS
    FOV = N(1:3).*MS;
    [kx, ky, kz]=linearTerms(N(1:3), N(1:3)/2./FOV);%*Nsusc/2 to compensate for scaling in linearTerms.m and ./FOV to account for anisotropy

    %%% TAKE INTO ACCOUNT IMAGE ORIENTATION
    R_MT = MT(1:3,1:3)/diag(MS);%Need to apply forward since defined from RAS-->ijk and need to rotate ijk to RAS

    K = [kx(:)';ky(:)'; kz(:)'];
    K = inv( inv(R_MT) )  * K;
    kx = reshape(K(1,:).',N(1:3));
    ky = reshape(K(2,:).',N(1:3));
    kz = reshape(K(3,:).',N(1:3));
