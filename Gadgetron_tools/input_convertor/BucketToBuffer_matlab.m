function twix = BucketToBuffer_matlab(bucket,connection_hdr,noiseData,saveJSON_flag,compulsory_traj)
% ZN: modified based on create_slice_from_bucket from gadgetron matlab
% toolbox demo
% ZN: to note, this currently only take care of the case with single
% contrast/echoes/invs... For the cases with multiple averages, this might
% not work
% ZN: currently this might only works for data with sprial trajectory,
% and other non-cartesian trejectories related code still needs to be
% tested
% current support: 
% - spiral data (single echo & single average)
%
% last modification: 12-Mar-2025
% by Zihan Ning @ King's College London

    if nargin < 3 || isempty(noiseData); noiseData=[];end 
    if nargin < 4 || isempty(saveJSON_flag); saveJSON_flag = 0; end 
    if nargin < 5 || isempty(compulsory_traj); compulsory_traj = []; end % ZN: in case the traj is encoded wrong in the twix header
    

    type_traj = connection_hdr.encoding.trajectory; 
    if ~isempty(compulsory_traj) && ~isequal(compulsory_traj,type_traj)
        fprintf('Although the trajectory is recognized as %s by the header, it is treated as %s in compulsory. \n',type_traj,compulsory_traj);
        type_traj = compulsory_traj;
    end
    fprintf('Converting data with %s trajectory from ismrmrd raw readin to a twix-like raw! \n',type_traj);
    
    %% Read related information

    disp("Assembling buffer from bucket containing " + num2str(bucket.data.count) + " acquisitions");
    % Extract information from large header (connection_hdr)
    matrix_size = connection_hdr.encoding.encodedSpace.matrixSize;
    kspace_encoding_Lin = connection_hdr.encoding.encodingLimits.kspace_encoding_step_1;
    kspace_encoding_Par = connection_hdr.encoding.encodingLimits.kspace_encoding_step_2;

    %% Read buffer data   
    switch type_traj
        case 'radial'
            %% ZN: initialization
            N_RO = size(bucket.data.data, 1);
            N_Cha = size(bucket.data.data, 2);
            N_Par = matrix_size.z;
            N_shot = size(bucket.data.data, 3)/matrix_size.z; % Number of spokes for radial
            N_Ave = max(bucket.data.header.average)+1; % averages
            N_Con = max(bucket.data.header.contrast)+1; % contrast/echo/inv/etc..
            if mod(N_shot,1)~=0 % ZN: sometimes extra lines will be sampled (normally at the beginning), we simply remove those lines here
                N_ExtraLine = mod(size(bucket.data.data, 3),matrix_size.z); % number of the extra lines 
                bucket.data.data = bucket.data.data(:,:,N_ExtraLine+1:end); % remove the extra lines (we think they are always at the beginning)
                N_shot = size(bucket.data.data, 3)/matrix_size.z; % update the N_interleaves
                fprintf('%s extra line has been removed in RO direction, and current number of interleaves is %s \n',num2str(N_ExtraLine),num2str(N_shot));
            end
            % ZN: here I tried to map the data matrix align with the twix
            % format read by mapVBVD
            %                 RO Cha Lin Par Sli Ave Phs Eco Rep Set Seg 
            twix.data = complex(zeros(... % ZN: currently this will not work for a multi-average or contrast (echos/invs...) scenario
                N_RO,... % RO
                N_Cha,... % CHA
                N_shot,... % number of spokes
                N_Par,... % Par
                ... % ZN: ignore the other dimensions for radial data for now...
                'single' ...
                ));
            
            %% ZN: map some dimension related headers 
            id = 1:length(bucket.data.header.kspace_encode_step_1); % ZN: currently this only works for single average & contrast data
                Lin = dynInd(bucket.data.header.kspace_encode_step_1,id,2);
                Par = dynInd(bucket.data.header.kspace_encode_step_2,id,2); 
                Eco = dynInd(bucket.data.header.contrast,id,2);  
                Ave = dynInd(bucket.data.header.average,id,2); 
            id=[];
            
            twix.image.NAve = N_Ave;
            twix.image.Lin = double(Lin + 1); % ZN: +1 since MRD format Lin start from 0, instead of 1
            twix.image.Par = double(Par + 1); % ZN: +1 since MRD format Lin start from 0, instead of 1
            if isempty(Eco);twix.image.Eco = ones(1,length(Lin));else;twix.image.Eco = Eco + 1; end
            if isempty(Ave);twix.image.Ave = ones(1,length(Lin));else;twix.image.Ave = Ave + 1; end
            
            % ZN: I think for spiral trajectory, these parameter might be
            % not that important, just hard coded here to keep consistency
            twix.image.centerLin = ones(1,length(Lin));
            twix.image.centerPar = ones(1,length(Lin));
            twix.image.centerCol = ones(1,length(Lin));
            
            % ZN: if you'd like to have encoded & recon MatrixSize & FOV_mm
            % I recommend you use the following parameters directly
            % instead of trying to read them from twix-like data since the
            % calculation is already done by the ParameterMap
            % considering RO OS removal or slice OS, etc
            % double check before you use
            twix.hdr.encodedSpace.matrixSize = [connection_hdr.encoding.encodedSpace.matrixSize.x,connection_hdr.encoding.encodedSpace.matrixSize.y,connection_hdr.encoding.encodedSpace.matrixSize.z];
            twix.hdr.encodedSpace.FOVmm = [connection_hdr.encoding.encodedSpace.fieldOfView_mm.x,connection_hdr.encoding.encodedSpace.fieldOfView_mm.y,connection_hdr.encoding.encodedSpace.fieldOfView_mm.z];
            twix.hdr.reconSpace.matrixSize = [connection_hdr.encoding.reconSpace.matrixSize.x,connection_hdr.encoding.reconSpace.matrixSize.y,connection_hdr.encoding.reconSpace.matrixSize.z]; 
            % ZN: to note, I've modified the parametermap for spiral to
            % make the reconSpace.matrixSize.y to be correct, double check!
            twix.hdr.reconSpace.FOVmm = [connection_hdr.encoding.reconSpace.fieldOfView_mm.x,connection_hdr.encoding.reconSpace.fieldOfView_mm.y,connection_hdr.encoding.reconSpace.fieldOfView_mm.z];
                              
            % ZN: kspace line header related to geom
            twix.image.slicePos = bucket.data.header.position;
            % ZN: for slicePos, somehow I cannot get
            % twix.image.slicePos(4:7,:) and it seems that we can use
            % read_dir... to get something similar for geom computation
            % instead
            twix.hdr.Dicom.tPatientPosition = connection_hdr.measurementInformation.patientPosition; 
            % some extra headers related to geom will be user-defined
            % parameters, see below
            
            % ZN: following are why we need slicePos(4:7) and how we can
            % use another way to do the calculation...
            % ZN: currently I took hdr.phase_dir/read_dir/slice_dir from
            % MRD header - the geom computed based on this might have small
            % shift compared with scanner reconstruction, but they should
            % be consist inherently
            twix.image.read_dir = bucket.data.header.read_dir; % ZN: they should be the same for all the points, so just take the first one for geom computation should be enough
            twix.image.phase_dir = bucket.data.header.phase_dir;
            twix.image.slice_dir = bucket.data.header.slice_dir;
            
            % ZN: then for geom computation, replaced the above code (if
            % you used to working with twix raw) with the code below
            % (start to work with twix-like structure converted from ISMRMRD)
            
%             %     %ROTATION (when you compute geom with twix)
%         %     quaternionRaw = slicePos(4:7);
%         %     rec.Par.Mine.Arot = eye(4);
%         %     rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS
%         % 
%         %     %TRANSLATION
%         %     translationRaw = slicePos(1:3); %in mm for center FOV (I think)
%         %     rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';
% 
%             %ROTATION (when you compute geom with twix_like structure converted from MRD)
%             rec.Par.Mine.Arot = double( eye(4));
%             rec.Par.Mine.Arot(1:3,1:3)=double( cat(2, twix.image.read_dir(:,1), twix.image.phase_dir(:,1), twix.image.slice_dir(:,1))); %RO-PE-SL to PCS
%             Arot_tmp = rec.Par.Mine.Arot;
%             rec.Par.Mine.Arot(1:3,1) = Arot_tmp(1:3,2);rec.Par.Mine.Arot(1:3,2) = Arot_tmp(1:3,1); 
%             clear Arot_tmp
            
            % ZN: sequence related parameters
            twix.hdr.Meas.alTE = connection_hdr.sequenceParameters.TE*1e3; % TE
            twix.hdr.Meas.alTR = connection_hdr.sequenceParameters.TR*1e3; % TR
            twix.hdr.Dicom.adFlipAngleDegree = connection_hdr.sequenceParameters.flipAngle_deg; % FA
            %% ZN: Bucket2Buffer
            % handling Lin & Par & Contrast & Average 
            % ZN: currently just checked that Lin & Par works; better to check with multiple average or contrast/echo data later
            if exist('N_ExtraLine','var') && ~isempty(N_ExtraLine) && N_ExtraLine~=0 % if we've removed lines from bucket, the Lin & Par need to be aligned as well
                Lin = Lin(N_ExtraLine+1:end); % Lin
                Par = Par(N_ExtraLine+1:end); % Par
            end
            
            for i = 1:length(Lin)
                encode_step_1 = Lin(i);
                encode_step_2 = Par(i);
                twix.data(:, :, encode_step_1 + 1, encode_step_2 + 1) = permute(transpose(squeeze(bucket.data.data(:, :, i))),[2 1]);            
            end
                   
            %% ZN: map some added parameters (Userdefined parameters in ParameterMap)
            % ZN: user defined parameters will be saved in the field 'userParameters'
            % by its type (Double, Long, or String)
            
            % Double
            if isfield(connection_hdr.userParameters,'userParameterDouble')
                hdr_UP.Double = connection_hdr.userParameters.userParameterDouble;
                [twix.hdr.Protocol.dOverallImageScaleFactor,status] = searchUPfield(hdr_UP,'Double','UP_dOverallImageScaleFactor');
                [twix.hdr.Meas.flReadoutOSFactor,status] = searchUPfield(hdr_UP,'Double','UP_flReadoutOSFactor');
                [twix.hdr.Dicom.dPhaseFOV,status] = searchUPfield(hdr_UP,'Double','dPhaseFOV');
            end
            
            % Long
            if isfield(connection_hdr.userParameters,'userParameterLong')
                hdr_UP.Long = connection_hdr.userParameters.userParameterLong;
                [twix.hdr.MeasYaps.sWipMemBlock.alFree{32},status] = searchUPfield(hdr_UP,'Long','UP_sWipMemBlock_alFree_32');
                [twix.hdr.Meas.alDwellTime,status] = searchUPfield(hdr_UP,'Long','UP_dwellTime'); % ZN: I only took out the alDwellTime{1} from twix header here
%                 [twix.hdr.Meas.alTE,status] = searchUPfield(hdr_UP,'Long','UP_alTE'); % ZN: I only took out the alTE{1} from the twix header; current version might not used for multiecho cases
                [twix.hdr.Dicom.lGlobalTablePosCor,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosCor'); % geom
                [twix.hdr.Dicom.lGlobalTablePosTra,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosTra'); % geom
                [twix.hdr.Dicom.lGlobalTablePosSag,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosSag'); % geom
                [twix.hdr.Meas.iNoOfFourierLines,status] = searchUPfield(hdr_UP,'Long','iNoOfFourierLines'); 
                [twix.hdr.Meas.lImagesPerSlab,status] = searchUPfield(hdr_UP,'Long','lImagesPerSlab'); 
            end
            
            % String
            if isfield(connection_hdr.userParameters,'userParameterString')
                hdr_UP.String = connection_hdr.userParameters.userParameterString;
                [twix.hdr.Dicom.SoftwareVersions,status] = searchUPfield(hdr_UP,'String','SoftwareVersions'); % software version
            end
            
            %% ZN: handling the header for reconstructed image
            % usually since the scanner cannot provide reconstructions for
            % such non-cartesian imaging, we will use the metadata of the
            % retrievall dummy scan (which will be a dummy scan with
            % Cartesian traj with the matched matrix size of the
            % reconstructed image)

            % ZN: we maintain the small header for each kspace line just in case
            twix.bucket_header = bucket.data.header;
        case 'spiral'
            %% ZN: initializaton 
            N_RO = size(bucket.data.data, 1);
            N_Cha = size(bucket.data.data, 2);
            N_Par = matrix_size.z;
            N_shot = size(bucket.data.data, 3)/matrix_size.z; % N_interleaves
            N_Ave = max(bucket.data.header.average)+1; % averages
            N_Con = max(bucket.data.header.contrast)+1; % contrast/echo/inv/etc..
            if mod(N_shot,1)~=0 % ZN: sometimes extra lines will be sampled (normally at the beginning), we simply remove those lines here
                N_ExtraLine = mod(size(bucket.data.data, 3),matrix_size.z); % number of the extra lines 
                bucket.data.data = bucket.data.data(:,:,N_ExtraLine+1:end); % remove the extra lines (we think they are always at the beginning)
                N_shot = size(bucket.data.data, 3)/matrix_size.z; % update the N_interleaves
                fprintf('%s extra line has been removed in RO direction, and current number of interleaves is %s \n',num2str(N_ExtraLine),num2str(N_shot));
            end
            % ZN: here I tried to map the data matrix align with the twix
            % format read by mapVBVD
            %                 RO Cha Lin Par Sli Ave Phs Eco Rep Set Seg 
            twix.data = complex(zeros(... % ZN: currently this will not work for a multi-average or contrast (echos/invs...) scenario
                N_RO,... % RO
                N_Cha,... % CHA
                N_shot,... % interleaves
                N_Par,... % Par
                ... % ZN: ignore the other dimensions for spiral data for now...
                'single' ...
                ));
            
            %% ZN: map some dimension related headers 
            id = 1:length(bucket.data.header.kspace_encode_step_1); % ZN: currently this only works for single average & contrast data
                Lin = dynInd(bucket.data.header.kspace_encode_step_1,id,2);
                Par = dynInd(bucket.data.header.kspace_encode_step_2,id,2); 
                Eco = dynInd(bucket.data.header.contrast,id,2);  
                Ave = dynInd(bucket.data.header.average,id,2); 
            id=[];
            
            twix.image.NAve = N_Ave;
            twix.image.Lin = double(Lin + 1); % ZN: +1 since MRD format Lin start from 0, instead of 1
            twix.image.Par = double(Par + 1); % ZN: +1 since MRD format Lin start from 0, instead of 1
            if isempty(Eco);twix.image.Eco = ones(1,length(Lin));else;twix.image.Eco = Eco + 1; end
            if isempty(Ave);twix.image.Ave = ones(1,length(Lin));else;twix.image.Ave = Ave + 1; end
            
            % ZN: I think for spiral trajectory, these parameter might be
            % not that important, just hard coded here to keep consistency
            twix.image.centerLin = ones(1,length(Lin));
            twix.image.centerPar = ones(1,length(Lin));
            twix.image.centerCol = ones(1,length(Lin));
            
            % ZN: if you'd like to have encoded & recon MatrixSize & FOV_mm
            % I recommend you use the following parameters directly
            % instead of trying to read them from twix-like data since the
            % calculation is already done by the ParameterMap
            % considering RO OS removal or slice OS, etc
            % double check before you use
            twix.hdr.encodedSpace.matrixSize = [connection_hdr.encoding.encodedSpace.matrixSize.x,connection_hdr.encoding.encodedSpace.matrixSize.y,connection_hdr.encoding.encodedSpace.matrixSize.z];
            twix.hdr.encodedSpace.FOVmm = [connection_hdr.encoding.encodedSpace.fieldOfView_mm.x,connection_hdr.encoding.encodedSpace.fieldOfView_mm.y,connection_hdr.encoding.encodedSpace.fieldOfView_mm.z];
            twix.hdr.reconSpace.matrixSize = [connection_hdr.encoding.reconSpace.matrixSize.x,connection_hdr.encoding.reconSpace.matrixSize.y,connection_hdr.encoding.reconSpace.matrixSize.z]; 
            % ZN: to note, I've modified the parametermap for spiral to
            % make the reconSpace.matrixSize.y to be correct, double check!
            twix.hdr.reconSpace.FOVmm = [connection_hdr.encoding.reconSpace.fieldOfView_mm.x,connection_hdr.encoding.reconSpace.fieldOfView_mm.y,connection_hdr.encoding.reconSpace.fieldOfView_mm.z];
                              
            % ZN: kspace line header related to geom
            twix.image.slicePos = bucket.data.header.position;
            % ZN: for slicePos, somehow I cannot get
            % twix.image.slicePos(4:7,:) and it seems that we can use
            % read_dir... to get something similar for geom computation
            % instead
            twix.hdr.Dicom.tPatientPosition = connection_hdr.measurementInformation.patientPosition; 
            % some extra headers related to geom will be user-defined
            % parameters, see below
            
            % ZN: following are why we need slicePos(4:7) and how we can
            % use another way to do the calculation...
            % ZN: currently I took hdr.phase_dir/read_dir/slice_dir from
            % MRD header - the geom computed based on this might have small
            % shift compared with scanner reconstruction, but they should
            % be consist inherently
            twix.image.read_dir = bucket.data.header.read_dir; % ZN: they should be the same for all the points, so just take the first one for geom computation should be enough
            twix.image.phase_dir = bucket.data.header.phase_dir;
            twix.image.slice_dir = bucket.data.header.slice_dir;
            
            % ZN: then for geom computation, replaced the above code (if
            % you used to working with twix raw) with the code below
            % (start to work with twix-like structure converted from ISMRMRD)
            
%             %     %ROTATION (when you compute geom with twix)
%         %     quaternionRaw = slicePos(4:7);
%         %     rec.Par.Mine.Arot = eye(4);
%         %     rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS
%         % 
%         %     %TRANSLATION
%         %     translationRaw = slicePos(1:3); %in mm for center FOV (I think)
%         %     rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';
% 
%             %ROTATION (when you compute geom with twix_like structure converted from MRD)
%             rec.Par.Mine.Arot = double( eye(4));
%             rec.Par.Mine.Arot(1:3,1:3)=double( cat(2, twix.image.read_dir(:,1), twix.image.phase_dir(:,1), twix.image.slice_dir(:,1))); %RO-PE-SL to PCS
%             Arot_tmp = rec.Par.Mine.Arot;
%             rec.Par.Mine.Arot(1:3,1) = Arot_tmp(1:3,2);rec.Par.Mine.Arot(1:3,2) = Arot_tmp(1:3,1); 
%             clear Arot_tmp
            
            
            %% ZN: Bucket2Buffer
            % handling Lin & Par & Contrast & Average 
            % ZN: currently just checked that Lin & Par works; better to check with multiple average or contrast/echo data later
            if exist('N_ExtraLine','var') && ~isempty(N_ExtraLine) && N_ExtraLine~=0 % if we've removed lines from bucket, the Lin & Par need to be aligned as well
                Lin = Lin(N_ExtraLine+1:end); % Lin
                Par = Par(N_ExtraLine+1:end); % Par
            end
            
            for i = 1:length(Lin)
                encode_step_1 = Lin(i);
                encode_step_2 = Par(i);
                twix.data(:, :, encode_step_1 + 1, encode_step_2 + 1) = permute(transpose(squeeze(bucket.data.data(:, :, i))),[2 1]);            
            end
                   
            %% ZN: map some added parameters (Userdefined parameters in ParameterMap)
            % ZN: user defined parameters will be saved in the field 'userParameters'
            % by its type (Double, Long, or String)
            % Double
            if isfield(connection_hdr.userParameters,'userParameterDouble')
                hdr_UP.Double = connection_hdr.userParameters.userParameterDouble;
                [twix.hdr.Protocol.dOverallImageScaleFactor,status] = searchUPfield(hdr_UP,'Double','UP_dOverallImageScaleFactor');
                [twix.hdr.Meas.flReadoutOSFactor,status] = searchUPfield(hdr_UP,'Double','UP_flReadoutOSFactor');
                [twix.hdr.Dicom.dPhaseFOV,status] = searchUPfield(hdr_UP,'Double','dPhaseFOV');
            end
            
            % Long
            if isfield(connection_hdr.userParameters,'userParameterLong')
                hdr_UP.Long = connection_hdr.userParameters.userParameterLong;
                [twix.hdr.MeasYaps.sWipMemBlock.alFree{32},status] = searchUPfield(hdr_UP,'Long','UP_sWipMemBlock_alFree_32');
                [twix.hdr.Meas.alDwellTime,status] = searchUPfield(hdr_UP,'Long','UP_dwellTime'); % ZN: I only took out the alDwellTime{1} from twix header here
                [twix.hdr.Meas.alTE,status] = searchUPfield(hdr_UP,'Long','UP_alTE'); % ZN: I only took out the alTE{1} from the twix header; current version might not used for multiecho cases
                [twix.hdr.Dicom.lGlobalTablePosCor,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosCor'); % geom
                [twix.hdr.Dicom.lGlobalTablePosTra,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosTra'); % geom
                [twix.hdr.Dicom.lGlobalTablePosSag,status] = searchUPfield(hdr_UP,'Long','UP_lGlobalTablePosSag'); % geom
                [twix.hdr.Meas.iNoOfFourierLines,status] = searchUPfield(hdr_UP,'Long','iNoOfFourierLines'); 
                [twix.hdr.Meas.lImagesPerSlab,status] = searchUPfield(hdr_UP,'Long','lImagesPerSlab'); 
            end
            
            % String
            if isfield(connection_hdr.userParameters,'userParameterString')
                hdr_UP.String = connection_hdr.userParameters.userParameterString;
                [twix.hdr.Dicom.SoftwareVersions,status] = searchUPfield(hdr_UP,'String','SoftwareVersions'); % software version
            end
            
            % ZN: we maintain the small header for each kspace line just in case
            twix.bucket_header = bucket.data.header;
            
        otherwise 
            %% ZN: this is for Cartesian trajectory (still need to be tidied up)
            % currently work for MEGE
            %% ZN: initializaton 
            N_RO = size(bucket.data.data, 1);
            N_Cha = size(bucket.data.data, 2);
%             N_Lin = matrix_size.y;
%             N_Par = matrix_size.z;
            N_Lin = kspace_encoding_Lin.maximum + 1; % ZN: handling the zero-padding afterwards, in the twix-like structure, we provide unpadded data
            N_Par = kspace_encoding_Par.maximum + 1;
            N_Ave = max(bucket.data.header.average)+1; % averages
            N_Con = max(bucket.data.header.contrast)+1; % contrast/echo/inv/etc..
            % ZN: here I tried to map the data matrix align with the twix
            % format read by mapVBVD
            %                 RO Cha Lin Par Sli Ave Phs Eco Rep Set Seg 
            twix.data = complex(zeros(... % ZN: currently this will not work for a multi-average
                N_RO,... % RO
                N_Cha,... % CHA
                N_Lin,... % Lin
                N_Par,... % Par
                N_Ave,... % Ave
                N_Con,... % contrast/echo/inv/etc...
                ... % ZN: ignore the other dimensions for spiral data for now...
                'single' ...
                ));
            
            %% ZN: map some dimension related headers    
            % ZN: extract sample info from kspace line header
            id = 1:length(bucket.data.header.kspace_encode_step_1); 
                Lin = double(dynInd(bucket.data.header.kspace_encode_step_1,id,2));
                Par = double(dynInd(bucket.data.header.kspace_encode_step_2,id,2)); 
                Con = double(dynInd(bucket.data.header.contrast,id,2));  
                Ave = double(dynInd(bucket.data.header.average,id,2)); 
            id=[];
            
            % ZN: rearrange the sample info according to the matrix of data
            % encoded_lines AVE CON
            n_rep = length(Lin) / (N_Con*N_Ave);
            sort_Lin = zeros(n_rep,N_Con,N_Ave);
            sort_Par = sort_Lin;
            sort_Con = sort_Lin;
            sort_Ave = sort_Lin;
            for i = 0:N_Con-1
                for j = 0:N_Ave-1
                    idx = find(Con == i & Ave == j);
                    if numel(idx) ~= n_rep; fprintf('Mismatch in expected number of repeatation for contrast %s, average %s.\n',num2str(N_Con),num2str(N_Ave)); break; end
                    sort_Lin(:,i+1,j+1) = Lin(idx);
                    sort_Par(:,i+1,j+1) = Par(idx);
                    sort_Con(:,i+1,j+1) = Con(idx);
                    sort_Ave(:,i+1,j+1) = Ave(idx);
                end
            end
            
            % ZN: copy the twix.image field info here
            twix.image.NAve = N_Ave;
            % ZN: copy both original and sorted sample info to the
            % twix-like structure (original twix has unsorted info)
            twix.image.Lin = Lin + 1; twix.image.sort_Lin = sort_Lin + 1; % ZN: +1 since MRD header start from 0 instead of 1
            twix.image.Par = Par + 1; twix.image.sort_Par = sort_Par + 1; 
            twix.image.Ave = Ave + 1; twix.image.sort_Ave = sort_Ave + 1; 
            twix.image.Eco = Con + 1; twix.image.sort_Eco = sort_Con + 1; 
            % ZN: the following part are calculated by parameterMap
            % if any errors, modify the parameter maps instead
            % if any accelerations in Lin or Par, there could be 1 pixel mismatch
            twix.image.centerLin = (connection_hdr.encoding.encodingLimits.kspace_encoding_step_1.center + 1)*ones(1,length(Lin));
            twix.image.centerPar = (connection_hdr.encoding.encodingLimits.kspace_encoding_step_2.center + 1)*ones(1,length(Par));
            twix.image.centerCol = double(max(bucket.data.header.center_sample) + 1)*ones(1,length(Lin)); % ZN: to make centerCol identical to the twix raw, you should not remove RO OS with the default gadget RemoveROOversampling
            
            % ZN: if you'd like to have encoded & recon MatrixSize & FOV_mm
            % I recommend you use the following parameters directly
            % instead of trying to read them from twix-like data since the
            % calculation is already done by the ParameterMap
            % considering RO OS removal or slice OS
            % double check before you use
            % acquisition matrix size (without zero-padding, with OS)
            dataSize = size(twix.data); % ZN: [Col Cha Lin Par Ave Con]; 
            twix.hdr.acqSpace.matrixSize = [dataSize(1),dataSize(3),dataSize(4)];
            twix.hdr.acqSpace.NCha = dataSize(2);
            % encoded grid (with OS, including zero-padded points due to partial fourier or asymmetric echo)
            twix.hdr.encodedSpace.matrixSize = [connection_hdr.encoding.encodedSpace.matrixSize.x,connection_hdr.encoding.encodedSpace.matrixSize.y,connection_hdr.encoding.encodedSpace.matrixSize.z];
            twix.hdr.encodedSpace.FOVmm = [connection_hdr.encoding.encodedSpace.fieldOfView_mm.x,connection_hdr.encoding.encodedSpace.fieldOfView_mm.y,connection_hdr.encoding.encodedSpace.fieldOfView_mm.z];
            % recon grid (without OS, including zero-padded points due to
            % PF and asymmetric echo, considering reconstruct the image to
            % a higher resolution)
            % for reconstruct to a higher resolution, it means that the PE resolution in sequence parameter map is not 100%
            twix.hdr.reconSpace.matrixSize = [connection_hdr.encoding.reconSpace.matrixSize.x,connection_hdr.encoding.reconSpace.matrixSize.y,connection_hdr.encoding.reconSpace.matrixSize.z]; 
            twix.hdr.reconSpace.FOVmm = [connection_hdr.encoding.reconSpace.fieldOfView_mm.x,connection_hdr.encoding.reconSpace.fieldOfView_mm.y,connection_hdr.encoding.reconSpace.fieldOfView_mm.z];
    
            
      
            % ZN: kspace line header related to geom
            twix.image.slicePos = bucket.data.header.position;
            % ZN: for slicePos, somehow I cannot get
            % twix.image.slicePos(4:7,:) and it seems that we can use
            % read_dir... to get something similar for geom computation
            % instead
            twix.hdr.Dicom.tPatientPosition = connection_hdr.measurementInformation.patientPosition; 
            % some extra headers related to geom will be user-defined
            % parameters, see below
            
            % ZN: following are why we need slicePos(4:7) and how we can
            % use another way to do the calculation...
            % ZN: currently I took hdr.phase_dir/read_dir/slice_dir from
            % MRD header - the geom computed based on this might have small
            % shift compared with scanner reconstruction, but they should
            % be consist inherently
            twix.image.read_dir = bucket.data.header.read_dir; % ZN: they should be the same for all the points, so just take the first one for geom computation should be enough
            twix.image.phase_dir = bucket.data.header.phase_dir;
            twix.image.slice_dir = bucket.data.header.slice_dir;
            
            % ZN: then for geom computation, replaced the above code (if
            % you used to working with twix raw) with the code below
            % (start to work with twix-like structure converted from ISMRMRD)
            
%             %     %ROTATION (when you compute geom with twix)
%         %     quaternionRaw = slicePos(4:7);
%         %     rec.Par.Mine.Arot = eye(4);
%         %     rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS
%         % 
%         %     %TRANSLATION
%         %     translationRaw = slicePos(1:3); %in mm for center FOV (I think)
%         %     rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';
% 
%             %ROTATION (when you compute geom with twix_like structure converted from MRD)
%             rec.Par.Mine.Arot = double( eye(4));
%             rec.Par.Mine.Arot(1:3,1:3)=double( cat(2, twix.image.read_dir(:,1), twix.image.phase_dir(:,1), twix.image.slice_dir(:,1))); %RO-PE-SL to PCS
%             Arot_tmp = rec.Par.Mine.Arot;
%             rec.Par.Mine.Arot(1:3,1) = Arot_tmp(1:3,2);rec.Par.Mine.Arot(1:3,2) = Arot_tmp(1:3,1); 
%             clear Arot_tmp
            
            % ZN: add the hdr from connection_hdr as well
            twix.hdr.sequenceParameters = connection_hdr.sequenceParameters;
            twix.hdr.acquisitionSystemInformation = connection_hdr.acquisitionSystemInformation;
            twix.hdr.subjectInformation = connection_hdr.subjectInformation;
            twix.hdr.studyInformation = connection_hdr.studyInformation;
            twix.hdr.measurementInformation = connection_hdr.measurementInformation;
            twix.hdr.encoding = connection_hdr.encoding;
            
            %% ZN: Bucket2Buffer
            % ZN: divide the Bucket into dimensions according to contrasts and averages
            n_rep = length(Lin) / (N_Con*N_Ave);
            sort_bucket = zeros(N_RO,N_Cha,n_rep,N_Ave,N_Con);
            for contrast = 0:N_Con-1
                for average = 0:N_Ave-1
                    idx = find(Con == contrast & Ave == average);
                    if numel(idx) ~= n_rep; fprintf('Mismatch in expected number of repeatation for contrast %s, average %s.\n',num2str(N_Con),num2str(N_Ave)); break; end
                    sort_bucket(:,:,:,average+1,contrast+1) = bucket.data.data(:,:,idx);
                end
            end
            
            % ZN: to Cartesian grid
            % this data has already be zero-padded
            % if zero-padded in a wrong way, please modify the parameterMap
            % to correct this
            % twix.data: Col Cha Lin Par Ave Con
            % 
            % if you find a shift in Lin & Par direction, perhaps due to
            % acceleration, check the maximum of Lin or Par to see whether
            % it starts from 1 instead of 0
            for contrast = 1:N_Con
                for average = 1:N_Ave
                    for i = 1:size(sort_Lin,1)
                        encode_step_1 = sort_Lin(i,1);
                        encode_step_2 = sort_Par(i,1);
                        twix.data(:, :, encode_step_1 + 1, encode_step_2 + 1, average, contrast) = squeeze(sort_bucket(:, :, i,average, contrast));
                    end
                end
            end
            
            % ZN: for zero-padding due to Partial Fourier in Lin or Par
            % directions, we tend to handle it later
            % here calculate the pre and post padding points for Lin and Par directions
            % see function: zeroPadKspace.m for details
            % here is just for sampling_description (used for reconstructed
            % image header generation)
            kspace_PrePad = [0 0]; 
            kspace_PrePad(1) = centerIdx(matrix_size.y) - kspace_encoding_Lin.center - 1; % Pre-padding in Lin % -1 since it is start from 0 instead of 1
            kspace_PrePad(2) = centerIdx(matrix_size.z) - kspace_encoding_Par.center - 1; % Pre-padding in Par

            
            
            %% ZN: map some added parameters (Userdefined parameters in ParameterMap)
            % ZN: user defined parameters will be saved in the field 'userParameters'
            % by its type (Double, Long, or String)
            % If you cannot find the parameters you need, please modify the
            % ParameterMap to add them first
            hdr_UP.Double = connection_hdr.userParameters.userParameterDouble;
            hdr_UP.Long = connection_hdr.userParameters.userParameterLong;
            hdr_UP.String = connection_hdr.userParameters.userParameterString;
            
            % ZN: the way to add new parameters
%             [twix.hdr.Protocol.dOverallImageScaleFactor,status] = searchUPfield(hdr_UP,'Double','UP_dOverallImageScaleFactor');
            
            % Double
            [twix.hdr.Dicom.flReadoutOSFactor,status] = searchUPfield(hdr_UP,'Double','flReadoutOSFactor');
            [twix.hdr.Meas.flPhaseOS,status] = searchUPfield(hdr_UP,'Double','flPhaseOS_zihan');
            [twix.hdr.Meas.flSliceOS,status] = searchUPfield(hdr_UP,'Double','flSliceOS_zihan');
            [twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1},status] = searchUPfield(hdr_UP,'Double','alDwellTime_zihan');
            [twix.hdr.MeasYaps.sKSpace.dPhaseResolution,status] = searchUPfield(hdr_UP,'Double','dPhaseResolution');
            [twix.hdr.MeasYaps.sKSpace.dSliceResolution,status] = searchUPfield(hdr_UP,'Double','dSliceResolution');
             
            % Long
            [twix.hdr.Dicom.lGlobalTablePosTra,status] = searchUPfield(hdr_UP,'Long','lGlobalTablePosTra_zihan'); % geom
            [twix.hdr.MeasYaps.lScanTimeSec,status] = searchUPfield(hdr_UP,'Long','lScanTimeSec_zihan'); 
            [twix.hdr.MeasYaps.sKSpace.unReordering,status] = searchUPfield(hdr_UP,'Long','unReordering_zihan'); 
            [twix.hdr.Meas.ucEnableEllipticalScanning,status] = searchUPfield(hdr_UP,'Long','ucEnableEllipticalScanning');  % ACS mode: 2 - integrated; 4- seperate
            
            [twix.hdr.MeasYaps.sPat.ucPATMode,status] = searchUPfield(hdr_UP,'Long','ucPATMode');  % acceleration mode: 2- GRAPPA; 16 - CAIPI
            [twix.hdr.MeasYaps.sPat.ucRefScanMode,status] = searchUPfield(hdr_UP,'Long','ucRefScanMode');  % ACS mode: 2 - integrated; 4- seperate
            
            [twix.hdr.MeasYaps.sFastImaging.lTurboFactor,status] = searchUPfield(hdr_UP,'Long','lTurboFactor');  % tse releted
            [twix.hdr.MeasYaps.sFastImaging.lEchoTrainDuration,status] = searchUPfield(hdr_UP,'Long','lEchoTrainDuration');  % tse releted
            [twix.hdr.Meas.RFEchoTrainLength,status] = searchUPfield(hdr_UP,'Long','RFEchoTrainLength');  % tse releted
            
            [twix.hdr.MeasYaps.sWipMemBlock.alFree{2},status] = searchUPfield(hdr_UP,'Long','sWipMemBlock_alFree_2');  % WIP (user defined parameters in sequence)
            [twix.hdr.MeasYaps.sWipMemBlock.alFree{3},status] = searchUPfield(hdr_UP,'Long','sWipMemBlock_alFree_3');  % WIP (user defined parameters in sequence)
            
            [twix.hdr.Meas.alRegridMode(1),status] = searchUPfield(hdr_UP,'Long','alRegridMode');  
            
            [twix.hdr.Meas.lRefLinesPE,status] = searchUPfield(hdr_UP,'Long','lRefLinesPE');   % ACS line info
 
            
            % Stringsort_bucket = zeros(N_RO,N_Cha,n_rep,N_Ave,N_Con);
            [twix.hdr.Meas.tScanningSequence,status] = searchUPfield(hdr_UP,'String','tScanningSequence_zihan'); % sequency type
            
            %%% if more, add here...
            
            %% ZN: some other computed parameters that you'd like to add to the header
            % usually you just added the most-general headers here
            
            % correct the recon fov_mm
            % ZN: this should be fixed in the parameterMap, but now just
            % modify here temporarily
            % ZN: for the NX
            
            % If multiple contrasts, to distinguish whether multiple invs or echoes
            twix.image.NInv = 1; twix.image.NEco = 1; % single echo/inv GRE
            if ndims(twix.data)>=6
                if dataSize(6) > 1 && ~isempty(twix.hdr.Meas.tScanningSequence) % multi-contrast
                    isShotBased = strcmp(twix.hdr.Meas.tScanningSequence,'GR\IR');  
                    if isShotBased
                        twix.image.NInv = dataSize(6); twix.image.NEco = 1; % multiple invs
                    else
                        twix.image.NInv = 1; twix.image.NEco = dataSize(6); % multiple echoes
                    end
                else % single contrast
                    twix.image.NInv = 1; twix.image.NEco = 1; % single echo/inv GRE
                end
            end
            
            % acceleration related parameters
            % Try to put back the acceleration related parameters into
            % original field of twix, and at the same time, generate a new
            % field for the summary (easy to use)
            if ~isempty(twix.hdr.MeasYaps.sPat.ucPATMode) % acceleration flag & type
                switch twix.hdr.MeasYaps.sPat.ucPATMode
                    case 1
                        twix.hdr.acceleration.isAccel = 0; % No accel (1)
                    case 2
                        twix.hdr.acceleration.isAccel = 1; twix.hdr.acceleration.type = 'GRAPPA';
                    case 16
                        twix.hdr.acceleration.isAccel = 1; twix.hdr.acceleration.type = 'CAIPI';
                    otherwise
                        twix.hdr.acceleration.isAccel = 1; twix.hdr.acceleration.type = 'Unrecog'; twix.hdr.acceleration.ucPATMode = twix.hdr.MeasYaps.sPat.ucPATMode;
                end
            end
            if twix.hdr.acceleration.isAccel % accelerated
                twix.hdr.MeasYaps.sPat.lAccelFactPE = connection_hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
                twix.hdr.MeasYaps.sPat.lAccelFact3D = connection_hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2;
                twix.hdr.acceleration.R = [twix.hdr.MeasYaps.sPat.lAccelFactPE, twix.hdr.MeasYaps.sPat.lAccelFact3D]; % [Lin Par]
                twix.hdr.acceleration.calibrationMode = connection_hdr.encoding.parallelImaging.calibrationMode; % seperate; embeded (integrated); interleaved
                
            end
                     
            
            %% ZN: handling the header for reconstructed image
            % handling buffer.sampling_description (this usually used to align the
            % format for reconstructed image to be sent back to gadget chain)
            twix.sampling_description.encoded_fov = ... % encoded FOV (mm)
                [connection_hdr.encoding.encodedSpace.fieldOfView_mm.x/2,... % ZN: remove RO OS
                connection_hdr.encoding.encodedSpace.fieldOfView_mm.y,...
                connection_hdr.encoding.encodedSpace.fieldOfView_mm.z];
            twix.sampling_description.recon_fov = ... % recon FOV (mm)
                [connection_hdr.encoding.reconSpace.fieldOfView_mm.x,... 
                connection_hdr.encoding.reconSpace.fieldOfView_mm.y,...
                connection_hdr.encoding.reconSpace.fieldOfView_mm.z];
            twix.sampling_description.encoded_matrix = ... % encoded matrix Size
                [connection_hdr.encoding.encodedSpace.matrixSize.x/2,...  % ZN: remove RO OS
                connection_hdr.encoding.encodedSpace.matrixSize.y,...
                connection_hdr.encoding.encodedSpace.matrixSize.z];
            twix.sampling_description.recon_matrix = ... % recon matrix Size
                [connection_hdr.encoding.reconSpace.matrixSize.x,...  
                connection_hdr.encoding.reconSpace.matrixSize.y,...
                connection_hdr.encoding.reconSpace.matrixSize.z];
            twix.sampling_description.sampling_limits.min = ... % kspace encoding limit (min)
                [0,...  % min of encoding idx of RO [warning: this might not work for the asymmetric echo case]
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_1.minimum + kspace_PrePad(1,1),... % Lin direction; -1 for start from 0 instead of 1
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_2.minimum + kspace_PrePad(1,2)]; % similar above
            twix.sampling_description.sampling_limits.center = ... % kspace encoding limit (centre)
                [centerIdx(connection_hdr.encoding.reconSpace.matrixSize.x) - 1,...  % centre of encoding idx of RO [warning: this might not work for the asymmetric echo case]; -1 since start from 0
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_1.center + kspace_PrePad(1,1),... % Lin direction; -1 for start from 0 instead of 1
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_2.center + kspace_PrePad(1,2)]; % similar above
            twix.sampling_description.sampling_limits.max = ... % kspace encoding limit (max)
                [connection_hdr.encoding.reconSpace.matrixSize.x - 1,...  % max of encoding idx of RO [warning: this might not work for the asymmetric echo case]; -1 since start from 0
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_1.maximum + kspace_PrePad(1,1),... % Lin direction; -1 for start from 0 instead of 1
                connection_hdr.encoding.encodingLimits.kspace_encoding_step_2.maximum + kspace_PrePad(1,2)]; % similar above

            % ZN: we maintain the small header for each kspace line just in case
            twix.bucket_header = bucket.data.header;
        
    end
       
    %% Read ref data
    % ZN: currently this is not align with the twix like structure
    % and just an initial version works for Cartesian
    
    readref_flag = 0;
    if bucket.reference.count>0; readref_flag = 1; end % ZN: there's reference/calibration lines to be read
    if readref_flag
        Num_refLin = max(bucket.reference.header.kspace_encode_step_1) + 1;
        Num_refPar = max(bucket.reference.header.kspace_encode_step_2) + 1;
        Num_refCon = max(bucket.reference.header.contrast) + 1;
        Num_refAve = max(bucket.reference.header.average) + 1;
        twix.refscan = complex(zeros( ...
            size(bucket.reference.data, 1), ... % RO
            size(bucket.reference.data, 2), ... % CHA
            Num_refLin, ... % PE1, Lin
            Num_refPar, ... % PE2, Par
            Num_refAve,... % Ave
            Num_refCon,... % contrast/echo/inv/etc...
            'single' ...
        ));
    
        id = 1:length(bucket.reference.header.kspace_encode_step_1); 
            refLin = double(dynInd(bucket.reference.header.kspace_encode_step_1,id,2));
            refPar = double(dynInd(bucket.reference.header.kspace_encode_step_2,id,2)); 
            refCon = double(dynInd(bucket.reference.header.contrast,id,2));  
            refAve = double(dynInd(bucket.reference.header.average,id,2)); 
        id=[];
        
        n_rep = length(refLin) / (Num_refCon*Num_refAve);
        sort_refLin = zeros(n_rep,Num_refCon,Num_refAve);
        sort_refPar = sort_refLin;
        sort_refCon = sort_refLin;
        sort_refAve = sort_refLin;
        for i = 0:Num_refCon-1
            for j = 0:Num_refAve-1
                idx = find(refCon == i & refAve == j);
                if numel(idx) ~= n_rep; fprintf('Mismatch in expected number of repeatation for contrast %s, average %s.\n',num2str(N_Con),num2str(N_Ave)); break; end
                sort_refLin(:,i+1,j+1) = refLin(idx);
                sort_refPar(:,i+1,j+1) = refPar(idx);
                sort_refCon(:,i+1,j+1) = refCon(idx);
                sort_refAve(:,i+1,j+1) = refAve(idx);
            end
        end
            
        % ref: Bucket2Buffer
        n_rep = length(refLin) / (Num_refCon*Num_refAve);
        sort_refbucket = zeros(size(bucket.reference.data, 1),size(bucket.reference.data, 2),n_rep,Num_refAve,Num_refCon);
        for contrast = 0:Num_refCon-1
            for average = 0:Num_refAve-1
                refidx = find(refCon == contrast & refAve == average);
                if numel(refidx) ~= n_rep; fprintf('Mismatch in expected number of repeatation for contrast %s, average %s.\n',num2str(Num_refCon),num2str(Num_refCon)); break; end
                sort_refbucket(:,:,:,average+1,contrast+1) = bucket.reference.data(:,:,refidx);
            end
        end
        for contrast = 1:Num_refCon
            for average = 1:Num_refAve
                for i = 1:size(sort_refLin,1)
                    encode_step_1 = sort_refLin(i,1);
                    encode_step_2 = sort_refPar(i,1);
                    twix.refscan(:, :, encode_step_1 + 1, encode_step_2 + 1, average, contrast) = squeeze(sort_refbucket(:, :, i,average, contrast));
                end
            end
        end
        mask2 = squeeze( any( twix.refscan~=0, [1,3,4,5,6] ) );   % might be some zero rows in the PE directions, crop them
        mask3 = squeeze( any( twix.refscan~=0, [1,2,4,5,6] ) );   
        twix.refscan = twix.refscan( :, mask2, mask3, :, :, : ); 
        twix.refscan = twix.refscan(:,:,:,:,1,1); % usually just take out the first contrast/average as the reference
    end
     
    % ZN: dimension of read-out ref [Col, Cha, Lin, Par]
    
    %%%% usually, we dont need the header from reference, but if any, added
    %%%% here
    
    %% Read noise lines (if any)
    if ~isempty(noiseData); twix.noise = noiseData; end
    
    %% Save JSON for debugging
    if saveJSON_flag
        twixJSON = twix;
        twixJSON.data = [];
        if isfield(twixJSON,'noise');twixJSON.noise = []; end
        if isfield(twixJSON,'image');twixJSON.image = []; end
        if isfield(twixJSON,'bucket_header');twixJSON.bucket_header = [];end
        if isfield(twixJSON,'reference');twixJSON.reference = [];end
    %     twixJSON.image.Lin = [];
    %     twixJSON.image.Par = [];
    %     twixJSON.image.Eco = [];
    %     twixJSON.image.Ave = [];
    %     twixJSON.image.centerLin = twixJSON.image.centerLin(:,1);
    %     twixJSON.image.centerPar = twixJSON.image.centerPar(:,1);
    %     twixJSON.image.centerCol = twixJSON.image.centerCol(:,1);
    %     twixJSON.image.slicePos = twixJSON.image.slicePos(:,1);
    %     twixJSON.image.read_dir = twixJSON.image.read_dir(:,1);
    %     twixJSON.image.phase_dir = twixJSON.image.phase_dir(:,1);
    %     twixJSON.image.slice_dir = twixJSON.image.slice_dir(:,1);

        savejson('',twixJSON,'JSON_MRD_converted_twix');%writeJSON has specific fields 
    end
end