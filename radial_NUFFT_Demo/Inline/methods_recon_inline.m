function [result,res,TE,shifts,shifts_TE1] = methods_recon_inline(twix, Tro, b1, selfcal, calpoints, shift, shift_echoes, shift_in,shift_in_TE1,shift_img,TE1only,meth_nufft,meth_script)
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
    % 4) str: for b1 mapping ('n','y')
    % 5) str: for sense self calibration ('n','y')
    % 6) int: number of points for self calibration
    % 7) str: compute trajectory shift ('n','y')
    % 8) int: number of echoes for computing shift
    % 9) array: precomputed trajectory shifts 
    % 10) array: precomputed trajectory shifts for TE1
    % 11) str: which data to use for shift computation
    % 12) str: only recon TE 1
    % 13) str: nufft method
    % 14) str: recon script version/method
    %-----------------------------------------------------------
    % Returns:
    % a) array: reconstructed images
    % b) single: resolution (mm ?)
    % c) single: echo time (s)
    % d) array: trajectory shifts if computed 
    % e) array: trajectory shifts for TE1 if computed
    %-----------------------------------------------------------
    % Notes: 
    %
    %-----------------------------------------------------------
    % This is an inline implemented version, which takes the 
    % converted twix-like raw within as the input
    % for the PIPE, please see /Inline/Inline_NUFFT_radial.m
    % 
    % Implemented by Zihan Ning (KCL)
    %===========================================================
        
    %% loading data and fetching some parameters
%     disp(['Loading Data from ', filename]);
%     measDat = mapVBVD(char(filename));
%     data = squeeze(measDat.image);
    % instead of reading by mapVBVD, use the converted twix-like raw
    measDat = twix; 
    data = measDat.data;

    if contains(measDat.hdr.Dicom.SoftwareVersions, 'XA60') % new XA60 header  
        numROs = measDat.hdr.Meas.iNoOfFourierLines;
        FOV = measDat.hdr.Dicom.dPhaseFOV*1e-3;
        MTX = measDat.hdr.Meas.lImagesPerSlab;
        TE = measDat.hdr.Meas.alTE(1)*1e-6;
        numEchoes = floor((size(data,1)-MTX)/(2*MTX));
    else
        numROs = measDat.hdr.Meas.NoOfFourierLines; 
        FOV = measDat.hdr.Meas.PeFOV*1e-3;
        MTX = measDat.hdr.Meas.NImageCols;
        TE = measDat.hdr.Meas.alTE(1)*1e-6;
        numEchoes = (data.dataSize(1)-MTX)/(2*MTX);
    end
    
    if b1 == "y" || TE1only == "y" %force use only 1st echo
        numEchoes = 0;
    end
    
    raw = data; % comes in 3D: [readpoints*TEs, channels, spokes]
    numChan = size(raw,2);
    img_out = zeros(numChan,MTX,MTX,MTX,numEchoes+1);
    
    
%     %%% averaging of k-spaces
%     if not(isempty(filename2))
%         measDat2 = mapVBVD(char(filename2));
%         raw2 = measDat2.image();
%         raw = (raw + raw2) ./ 2;
%         clear raw2
%     end

    % follow yasmin's recon routine
    if contains(meth_script, 'Yasmin')

        for ch = 1:numChan
            disp(['Utilising recon script: ', char(meth_script)]);
            disp(['Reconstructing images for channel ',num2str(ch)]);
            
            dat = squeeze(data(:,ch,:));

            tic; [result, TEs, res] = regrid_multiEcho(permute(dat,[2,1]), FOV, MTX, numROs, Tro, numEchoes, TE);
            img_out(ch,:,:,:,:) = result;
        end


    elseif meth_script == "Sam"

        disp(['Utilising recon script: ', char(meth_script)]);       
    
        %raw = ones(size(raw)); %uncomment for psf
        raw_TE1 = raw(1:MTX,:,:); % pick first TE
        raw_TE1 = permute(raw_TE1,[1 3 2]);
        if shift == "y"
            [dr_TE1,dx_TE1,dy_TE1,dz_TE1] = methods_eccshift_TE1(raw_TE1, shift_img);
            disp('Calculating trajectory shifts for TE1')
        elseif not(isempty(shift_in_TE1))
            dr_TE1 = shift_in_TE1(:,:,:,1);
            dx_TE1 = shift_in_TE1(:,:,:,2);
            dy_TE1 = shift_in_TE1(:,:,:,3);
            dz_TE1 = shift_in_TE1(:,:,:,4);
            disp('Using precalculated trajectory shifts for TE1')
        else
            dx_TE1 = [];
            disp('Not using trajectory shifts for TE1')
        end
        
        %%% if we are reconstructing images for SENSE self-calibration 
        if selfcal=='y'
            raw_TE1(calpoints:end,:,:) = 0 + 0j;
        end
        
        if numEchoes > 0
            raw_TErest = raw(MTX+1:end,:,:);
            raw_TErest = reshape(raw_TErest,[2*MTX,numEchoes,numChan,numROs]);
            raw_TErest = permute(raw_TErest,[1,4,3,2]);
            %calculate ecc shift
            if shift == "y"
                disp('Calculating trajectory shifts for TErest')
                [dr,dx,dy,dz] = methods_eccshift(raw_TErest,shift_echoes, shift_img);
            elseif not(isempty(shift_in))
                disp('Using precalculated trajectory shifts for TErest')
                dr = shift_in(:,:,:,:,1);
                dx = shift_in(:,:,:,:,2);
                dy = shift_in(:,:,:,:,3);
                dz = shift_in(:,:,:,:,4);
            else
                dx = [];
                disp('Not using trajectory shifts for TErest')
            end
        end
                %% global parameters
        NUFFT       = 0; 
        rampUpTime  = 0.1; % TODO: rut from WIP parameter!!!!
        steps = MTX;
        smplPts = linspace(Tro/steps,Tro,steps);
        rampSmplPts = smplPts(smplPts<rampUpTime)/rampUpTime;
        res = FOV/MTX;

        %% calc k-space information
        [Gx,Gy,Gz] = calcGradWvfrm(numROs,steps,rampSmplPts);
        [kx,ky,kz] = calcKspace(Gx,Gy,Gz,NUFFT);
        
        % scale to max pi for finufft
        kx = pi.*kx./0.5;
        ky = pi.*ky./0.5;
        kz = pi.*kz./0.5;
        w = calcWeights(kx,ky,kz)';
        
        % apply shifts
        if not(isempty(dx_TE1))
            kxs = kx - 0.10 .* dx_TE1(:,2); %use 0.1 fudge factor here to scale
            kys = ky - 0.10 .* dy_TE1(:,2); %shift accounting for ramp sampling
            kzs = kz - 0.10 .* dz_TE1(:,2);
            for spoke=1:numROs %shift the compensation weights as well
                ds = - squeeze(0.10.*dr_TE1(spoke,2)) .* (MTX / (pi));
                weights(:,spoke) = fraccircshift(w(:,spoke),ds);
            end
        else
            kxs = kx;
            kys = ky;
            kzs = kz;
            weights = w;
        end
        
        % do some reshaping
        kxs = permute(kxs,[2,1]);
        kys = permute(kys,[2,1]);
        kzs = permute(kzs,[2,1]);
        kxs = squeeze(reshape(kxs,numROs*MTX,1));
        kys = squeeze(reshape(kys,numROs*MTX,1));
        kzs = squeeze(reshape(kzs,numROs*MTX,1));
        coords = [kxs,kys,kzs];
        raw_TE1 = squeeze(reshape(raw_TE1,[numROs*MTX,numChan]));
        weights = reshape(weights,[numROs*MTX,1]);
        weights = squeeze(repmat(weights,[1,numChan]));
        
        if meth_nufft == "finufft"
            out = finufft3d1(kxs, kys, kzs, double(squeeze(raw_TE1.*weights)),-1,1e-8,MTX,MTX,MTX);      
        elseif meth_nufft == "matmri"
            Nop = nufftOp([MTX,MTX,MTX],coords * (0.5/pi),weights, false);
            out = Nop' * raw_TE1;
        end
    
        %add channel dimension for volume images
        if ndims(out) < 4
            img_out(:,:,:,:,1) = reshape(out,[1,size(out)]); 
        else
            img_out(:,:,:,:,1) = permute(out,[4,1,2,3]);
        end
        disp(['Reconstructed echo 1 out of ', num2str(numEchoes+1)])
        clear Nop
        if numEchoes > 0
        %% calc k-space information
            [Gx,Gy,Gz] = calcGradWvfrm(numROs,2.*steps,rampSmplPts);
            [kx,ky,kz] = calcKspaceRef(Gx,Gy,Gz,NUFFT);
            kx = pi.*kx./0.5;
            ky = pi.*ky./0.5;
            kz = pi.*kz./0.5;        
            w = calcWeightsRef(kx,ky,kz)';
            %flip the odd echoes
            for echo=1:numEchoes
                if mod(echo,2) == 0
                    raw_TErest(:,:,:,echo) = flip(raw_TErest(:,:,:,echo),1);
                end
            end
            for echo = 1:numEchoes
                for chan=1:numChan
                    clear weights
                    if not(isempty(dx))
                        if size(dx,2) == 1
                            chanshift = 1;
                        else
                            chanshift = chan;
                        end
                        kxs = kx + 1.0 .* dx(:,chanshift,echo,2);
                        kys = ky + 1.0 .* dy(:,chanshift,echo,2);
                        kzs = kz + 1.0 .* dz(:,chanshift,echo,2);
                        for spoke=1:numROs
                            ds = squeeze(1.0.*dr(spoke,chanshift,echo,2)) .* (2*MTX / (2*pi));
                            weights(:,spoke) = fraccircshift(w(:,spoke),ds);
                        end
                    else
                        kxs = kx;
                        kys = ky;
                        kzs = kz;
                        weights = w;
                    end
                    kxs = permute(kxs,[2,1]);
                    kys = permute(kys,[2,1]);
                    kzs = permute(kzs,[2,1]);
                    kxs = squeeze(reshape(kxs,numROs*2*MTX,1));
                    kys = squeeze(reshape(kys,numROs*2*MTX,1));
                    kzs = squeeze(reshape(kzs,numROs*2*MTX,1));
                    weights = reshape(weights,[numROs*2*MTX,1]);
                    %weights = repmat(weights,[1,numChan]);

                    raw = reshape(raw_TErest(:,:,chan,echo),[numROs*2*MTX,1]);
                    if meth_nufft == "finufft"
                        out = finufft3d1(kxs, kys, kzs, double(squeeze(raw.*(weights./2))),-1,1e-8,MTX,MTX,MTX);
                    elseif meth_nufft == "matmri"
                        Nop = nufftOp([MTX,MTX,MTX],[kxs,kys,kzs] * (0.5/pi),weights./2, false);
                        out = Nop' * raw;
                    end
    %                 if ndims(out) < 4
    %                     img_out(:,:,:,:,echo+1) = reshape(out,[1,size(out)]);
    %                 else
    %                     img_out(:,:,:,:,echo+1) = permute(out,[4,1,2,3]);
    %                 end
                    img_out(chan,:,:,:,echo+1) = out;
                end
                disp(['Reconstructed echo ', num2str(echo+1), ' out of ', num2str(numEchoes+1)])
            end
            
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