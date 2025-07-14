function [result,EchoTimes,res] = regrid_multiEcho(data, FOV, MTX, proj, Tro, numEchoes, TE)

    
    %% global parameters
    NUFFT       = 0; 
    rampUpTime  = 0.1; % TODO: rut from WIP parameter!!!!
    steps = MTX;
    smplPts = linspace(Tro/steps,Tro,steps);
    rampSmplPts = smplPts(smplPts<rampUpTime)/rampUpTime;
    res = FOV/MTX;
    

    %% calc k-space information
    [Gx,Gy,Gz] = calcGradWvfrm(proj,steps,rampSmplPts);
    [kx,ky,kz] = calcKspace(Gx,Gy,Gz,NUFFT);
    [weights]  = calcWeights(kx,ky,kz); 

    %% load data
    daten = [data(:,1:steps)];

    
    if NUFFT 
        np = [MTX,MTX,MTX];
        tic; FT = nufft_init([kx(:),ky(:),kz(:)],np,[3,3,3],2*np,np/2,'minmax:kb'); toc;
        gdata = nufft_adj(daten(:).*weights(:),FT);
    else
        dat = zeros([2,size(daten)]);
        dat(1,:,:) = real(daten);
        dat(2,:,:) = imag(daten);
        crds = zeros([3,size(kx)]);
        crds(1,:,:) = kx; crds(2,:,:) = ky; crds(3,:,:) = kz;    

        effMtx = ceil(1.5*MTX);
        numThread = 10;
        gdata = grid3_MAT(dat,crds,weights,effMtx,numThread);   
        gdata = rollKernel(gdata,MTX);
        
        disp(['finished reconstructing echo: 1/',num2str(numEchoes+1)]);
        disp('-----------------------------------------');
    end
    
    if (numEchoes>0) 
        [Gx,Gy,Gz] = calcGradWvfrm(proj,2*steps,rampSmplPts);
        [kx,ky,kz] = calcKspaceRef(Gx,Gy,Gz,NUFFT);%dt,gambar,gradMax,NUFFT,FOV,GxOld,GyOld,GzOld);     
        weights    = calcWeightsRef(kx,ky,kz);

        %% load data
        gdataEcho = zeros([size(gdata),numEchoes]);
        for echo = 1:numEchoes
            daten = [data(:,steps+1+(echo-1)*2*steps:steps+echo*2*steps)];

            if NUFFT 
                np = [MTX,MTX,MTX];
                tic; FT = nufft_init([kx(:),ky(:),kz(:)],np,[3,3,3],2*np,np/2,'minmax:kb'); toc;
                gdata2 = nufft_adj(daten(:).*weights(:)/2,FT);
            else
                dat = zeros([2,size(daten)]);
                dat(1,:,:) = real(daten);
                dat(2,:,:) = imag(daten);
                crds = zeros([3,size(kx)]);
                if (mod(echo,2))
                    crds(1,:,:) = kx; crds(2,:,:) = ky; crds(3,:,:) = kz;    
                else
                    crds(1,:,:) = -kx; crds(2,:,:) = -ky; crds(3,:,:) = -kz;     
                end

                effMtx = ceil(1.5*MTX);
                numThread = 10;
                gdata2 = grid3_MAT(dat,crds,weights/2,effMtx,numThread);
                gdata2 = rollKernel(gdata2,MTX);

                     end
            gdataEcho(:,:,:,echo) = gdata2;
            disp(['finished reconstructing echo: ', num2str(echo+1),'/',num2str(numEchoes+1)]);
            disp(['-----------------------------------------']);
        end
    else
        gdataEcho = [];
    end
    
    EchoTimes = zeros(numEchoes+1,1);
    
    for p = 1:numEchoes 
            EchoTimes(p+1) = p*2*Tro;
    end
    EchoTimes = EchoTimes + TE;
    
    result = cat(4,gdata,gdataEcho);
end