function [result, res, TE, shifts_TE1, shifts] = methods_recon(filename, filename2, Tro, TE1only, shift, shift_echoes, shift_in, shift_in_TE1, shift_img)
    %===========================================================
    % Main k-space to image function for 23Na MERINA multi-echo
    % reconstruction
    %-----------------------------------------------------------
    % Sam Rot (UCL) and Yasmin Blunck (University of Melbourne)
    %-----------------------------------------------------------
    % Takes up to 11 arguments: 
    % 1) str: filename of first image
    % 2) str: filename of second image if k-space averaging (leave 
    %    empty as '' if no) 
    % 3) single: readout window of data (ms)
    % 4) str: only recon TE 1
    % 5) str: compute trajectory shift ('n','y')
    % 6) int: number of echoes for computing shift
    % 7) array: precomputed trajectory shifts 
    % 8) array: precomputed trajectory shifts for TE1
    % 9) str: which data to use for shift computation
    %-----------------------------------------------------------
    % Returns:
    % a) array: reconstructed images
    % b) single: resolution (mm ?)
    % c) single: echo time (s)
    % d) array: trajectory shifts for TE1 if computed 
    % e) array: trajectory shifts if computed
    %-----------------------------------------------------------
    % Notes: 
    % computation of radial trajectory shift not supported in demo
    %===========================================================
        
    %% loading data and fetching some parameters
    disp(['Loading Data from ', filename]);
    measDat = mapVBVD(char(filename));
    data = squeeze(measDat.image);

    if contains(measDat.hdr.Dicom.SoftwareVersions, 'XA60') % new XA60 header  
        numROs = measDat.hdr.Meas.iNoOfFourierLines;
        FOV = measDat.hdr.Dicom.dPhaseFOV*1e-3;
        MTX = measDat.hdr.Meas.lImagesPerSlab;
        TE = measDat.hdr.Meas.alTE(1)*1e-6;
        numEchoes = floor((data.dataSize(1)-MTX)/(2*MTX)); % nr refocused echoes, i.e. half-spoke is echo 0
    else
        numROs = measDat.hdr.Meas.NoOfFourierLines; 
        FOV = measDat.hdr.Meas.PeFOV*1e-3;
        MTX = measDat.hdr.Meas.NImageCols;
        TE = measDat.hdr.Meas.alTE(1)*1e-6;
        numEchoes = (data.dataSize(1)-MTX)/(2*MTX); % nr refocused echoes, i.e. half-spoke is echo 0
    end
    
    if TE1only == "y" %force recon only 1st echo
        numEchoes = 0;
    end
    
    raw = measDat.image(); % comes in 3D: [readpoints*TEs, channels, spokes]
    numChan = size(raw,2);
    img_out = zeros(numChan,MTX,MTX,MTX,numEchoes+1);
    
    %% averaging of k-spaces
    if not(isempty(filename2))
        measDat2 = mapVBVD(char(filename2));
        raw2 = measDat2.image();
        raw = (raw + raw2) ./ 2;
        clear raw2
    end    
    
    % calculate shift of first echo (commented out, not supported in this demo)
    raw_TE1 = raw(1:MTX,:,:); % pick first TE
    raw_TE1 = permute(raw_TE1,[1 3 2]);
    %if shift == "y" % if shift is to be computed
    %    [dr_TE1,dx_TE1,dy_TE1,dz_TE1] = methods_eccshift_TE1(raw_TE1, shift_img);
    %    disp('Calculating trajectory shifts for TE1')
    %elseif not(isempty(shift_in_TE1)) % if shift was provided
    %    dr_TE1 = shift_in_TE1(:,:,:,1);
    %    dx_TE1 = shift_in_TE1(:,:,:,2);
    %    dy_TE1 = shift_in_TE1(:,:,:,3);
    %    dz_TE1 = shift_in_TE1(:,:,:,4);
    %    disp('Using precalculated trajectory shifts for TE1')
    %else % no traectory shifts
    dx_TE1 = [];
    disp('Not using trajectory shifts for TE1')
    %end

    if numEchoes > 0 % if we have a multiecho dataset
        raw_TErest = raw(MTX+1:end,:,:);
        raw_TErest = reshape(raw_TErest,[2*MTX,numEchoes,numChan,numROs]);
        raw_TErest = permute(raw_TErest,[1,4,3,2]);
        % calculate shift of remaining echoes (commented out, not supported in this demo)
        %if shift == "y" % if shift is to be computed
        %    disp('Calculating trajectory shifts for TErest')
        %    [dr,dx,dy,dz] = methods_eccshift(raw_TErest,shift_echoes, shift_img);
        %elseif not(isempty(shift_in)) % if shift was provided
        %    disp('Using precalculated trajectory shifts for TErest')
        %    dr = shift_in(:,:,:,:,1);
        %    dx = shift_in(:,:,:,:,2);
        %    dy = shift_in(:,:,:,:,3);
        %    dz = shift_in(:,:,:,:,4);
        %else % no traectory shifts
        dx = [];
        disp('Not using trajectory shifts for TErest')
        %end
    end

    %% some parameters for calculating k-space coordinates
    NUFFT = 0; 
    rampUpTime = 0.1; % TODO: rut from WIP parameter!!!!
    steps = MTX;
    smplPts = linspace(Tro/steps,Tro,steps);
    rampSmplPts = smplPts(smplPts<rampUpTime)/rampUpTime;
    res = FOV/MTX;

    %% calc k-space information
    [Gx,Gy,Gz] = calcGradWvfrm(numROs,steps,rampSmplPts);
    [kx,ky,kz] = calcKspace(Gx,Gy,Gz,NUFFT);
    
    % rescale from max 0.5 to max pi for finufft
    kx = pi.*kx./0.5;
    ky = pi.*ky./0.5;
    kz = pi.*kz./0.5;
    w = calcWeights(kx,ky,kz)'; % density compensation weights
    ws = zeros(size(w)); 
    
    %% apply k-space shifts if using (commented out for this demo)
    %if not(isempty(dx_TE1))
    %    kxs = kx - 0.10 .* dx_TE1(:,2); %use 0.1 fudge factor here to scale
    %    kys = ky - 0.10 .* dy_TE1(:,2); %shift accounting for ramp sampling
    %    kzs = kz - 0.10 .* dz_TE1(:,2);
    %    for spoke=1:numROs % shift the compensation weights as well
    %        ds = - squeeze(0.10.*dr_TE1(spoke,2)) .* (MTX / (pi));
    %        ws(:,spoke) = fraccircshift(w(:,spoke),ds);
    %    end
    %else
    kxs = kx;
    kys = ky;
    kzs = kz;
    ws = w;
    %end
    
    %% reshape and prepare for NUFFT
    kxs = permute(kxs,[2,1]);
    kys = permute(kys,[2,1]);
    kzs = permute(kzs,[2,1]);
    kxs = squeeze(reshape(kxs,numROs*MTX,1));
    kys = squeeze(reshape(kys,numROs*MTX,1));
    kzs = squeeze(reshape(kzs,numROs*MTX,1));
    raw_TE1 = squeeze(reshape(raw_TE1,[numROs*MTX,numChan]));
    ws = reshape(ws,[numROs*MTX,1]);
    ws = squeeze(repmat(ws,[1,numChan]));
    
    %% recon first echo using finufft
    out = finufft3d1(kxs, kys, kzs, double(squeeze(raw_TE1.*ws)),-1,1e-8,MTX,MTX,MTX);      

    % add channel dimension for volume coil data
    if ndims(out) < 4
        img_out(:,:,:,:,1) = reshape(out,[1,size(out)]); 
    else
        img_out(:,:,:,:,1) = permute(out,[4,1,2,3]);
    end
    disp(['Reconstructed echo 1 out of ', num2str(numEchoes+1)])

    if numEchoes > 0
        %% calc k-space information
        [Gx,Gy,Gz] = calcGradWvfrm(numROs,2.*steps,rampSmplPts);
        [kx,ky,kz] = calcKspaceRef(Gx,Gy,Gz,NUFFT);
        kx = pi.*kx./0.5;
        ky = pi.*ky./0.5;
        kz = pi.*kz./0.5;        
        w = calcWeightsRef(kx,ky,kz)';
        % flip the odd echoes (counting half spoke as nr. 1)
        for echo=1:numEchoes % here 1 is the first full spoke
            if mod(echo,2) == 0
                raw_TErest(:,:,:,echo) = flip(raw_TErest(:,:,:,echo),1);
            end
        end
        for echo = 1:numEchoes

            %% apply k-space shifts if using (commented out for this demo)
            ws = zeros(size(w)); 
            %if not(isempty(dx))
            %    kxs = kx + 1.0 .* dx(:,1,echo,2);
            %    kys = ky + 1.0 .* dy(:,1,echo,2);
            %    kzs = kz + 1.0 .* dz(:,1,echo,2);
            %    for spoke=1:numROs
            %        ds = squeeze(1.0.*dr(spoke,chanshift,echo,2)) .* (2*MTX / (2*pi));
            %        ws(:,spoke) = fraccircshift(w(:,spoke),ds);
            %    end
            %else
            kxs = kx;
            kys = ky;
            kzs = kz;
            ws = w;
            %end
            kxs = permute(kxs,[2,1]);
            kys = permute(kys,[2,1]);
            kzs = permute(kzs,[2,1]);
            kxs = squeeze(reshape(kxs,numROs*2*MTX,1));
            kys = squeeze(reshape(kys,numROs*2*MTX,1));
            kzs = squeeze(reshape(kzs,numROs*2*MTX,1));
            ws = reshape(ws,[numROs*2*MTX,1]);
            ws = repmat(ws,[1,numChan]);

            raw = reshape(raw_TErest(:,:,:,echo),[numROs*2*MTX,numChan]);
            % recon using finufft
            out = finufft3d1(kxs, kys, kzs, double(squeeze(raw.*(ws./2))),-1,1e-8,MTX,MTX,MTX);
            
            % add channel dimension in dim 1 for volume coil data
            if ndims(out) < 4
                img_out(:,:,:,:,echo+1) = reshape(out,[1,size(out)]);
            else
                img_out(:,:,:,:,echo+1) = permute(out,[4,1,2,3]);
            end
            disp(['Reconstructed echo ', num2str(echo+1), ' out of ', num2str(numEchoes+1)])
        end
    end
    
    result = permute(img_out,[2,3,4,1,5]); %permute to x,y,z,ch,echo
    
    %output shifts for use if calculated 
    if shift == 'y'
        shifts = cat(5,dr,dx,dy,dz);
        shifts_TE1 = cat(5,dr_TE1,dx_TE1,dy_TE1,dz_TE1);
    else
        shifts = [];
        shifts_TE1 = [];
    end
end