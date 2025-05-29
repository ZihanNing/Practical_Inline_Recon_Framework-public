
function handle_connection_disorder_MEGE_SENSE_subcall_bash
    % ZN: this is the SENSE (oridinary) recon pipeline
    
    opengl('save', 'none'); % stop showing X11 related warning
    seq_type = 'SWI';
    recon_type = 'SENSE'; % ZN: also has DISORDER recon version
    
    % ZN: write the log
    currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
    timestamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    diary_filename = sprintf('%s/log/subcall_log_SWI_%s_%s.txt', pwd, recon_type, timestamp);
    diary(diary_filename); 
    fprintf([currentTime,'\n']);
    fprintf('===================== Getting into the Subcall Function (MEGE) ====================\n');
    fprintf('===================== Performing SENSE reconstruction (MEGE) ======================\n');
    fprintf('Current folder:\n');
    disp(pwd);

    %%% Fix for CUBLAS issues
    warning off parallel:gpu:device:DeviceLibsNeedsRecompiling
    try
        gpuArray.eye(2)^2;
    catch; ME
    end
    % reset the gpu memory and reset gpu to solve the memory issue % by ZN
    gpuArray([]);
%     gpuDevice(1);

    %%% ADD SUPPORT FUNCTIONS PATH
    if ~exist('recon_type', 'var'); recon_type = 'SENSE';end
    if isequal(recon_type,'SENSE')
        addpath(genpath('Functions_SENSE')); % ZN: the support function box for SENSE recon
    else
        addpath(genpath('Functions')); % ZN: the support function box for DISORDER recon
    end 
    addpath(genpath('Control'));
    matlab_version_check  % ZN: in case of crashes due to matlab version conflicts

    %%% MULTI-GPU HANDLING 
    enable_multiGPU = 1; % ZN: enabling multi-GPU (only work on the server)
    gpuInfo = checkGPUAvailability(); % ZN: check GPU resources -> sometimes GPU driver update can cause disconnection
    if gpuInfo > 1
        if enable_multiGPU % ZN: enabling parallel compuation on multi-GPUs
            gpuIndex = ConfigurationGPU() % ZN: this will return the useable GPU in an ascending order
            gpuIndex_activate = gpuIndex(end);
            gpuDevice(gpuIndex_activate); % only activate one GPU
            fprintf('GPU %s has been activate!\n',sprintf('%d',gpuIndex_activate));
            updateGPUStatus(gpuIndex_activate, 'occupied'); % ZN: update the log file
        else % disable multiGPU compuation
            gpuDevice(1); % ZN: only run the code on GPU 1
        end
    elseif gpuInfo == 1
        enable_multiGPU = 0; % ZN: if there's only one GPU, disable multi-GPU calculation function
        warning('There is only one GPU accessable, the mutli-GPU calculation cannot be used!')
    else
        error('No GPU accessable - DISORDER recon aborted!')
    end


    %%% ADD PATHS AND REPOS
    try
        %%% LOAD THE ESSENTIAL INFO AND RAW DATA SAVED IN MAINCALL_BASH FUNC
        fprintf('=========== Loading the raw data.\n');
%         load('main_to_sub_SWI.mat') 
        [fileName, nameRef, patientID, save_path, save_time, sorting_info] = read_sortinginfo(seq_type); % ZN: the earliest registered data will be loaded for processing
        pause(5); load(nameRef); % 'bufferData','fileName','hdrConnection' are saved in nameRef
        flag_erase = erase_sortinginfo(sorting_info, save_path,patientID,fileName, save_time, seq_type); % ZN: the loaded data's sorting info will be wiped from the structure to prevent re-recon; corresponding main_to_sub file is updated
        if flag_erase
            fprintf('=========== Raw data successfully loaded.\n');
        else
            fprintf('!!!WARNING: ERROR IN RAW DATA READING AND SORTING INFO UPDATING!!!')
        end
    
%         % ZN: need to disable this timeline logging for the multi-GPU
%         % version, to prevent it abort the computation of other sequences
%         %%% PRINT THE TIME LINE (FOR ABORTION SECURITY)
%         logFilePath = fullfile(pwd, 'log', 'record_Disorder_timeline.txt');
%         check_timeline_for_abortion(logFilePath,fileName,patientID) % ZN: first check if there's an abortion
%         fileID = fopen(logFilePath, 'a');
%         if fileID ~= -1
%             currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
%             fprintf(fileID, [currentTime,'  ']);
%             fprintf(fileID, 'Start the MATLAB handle: disorder_subcall_bash. Sequence&Patient:%s&%s\n',fileName,patientID);
%             fclose(fileID);
%         end


        %%% ASSIGN NAMES & PARAMETERS 
        fprintf('=========== Assign names & parameters.\n');
        rec.Nam=struct('caseIn',save_path,'bodyIn',[],'surfIn',fileName,'dataIn',fileName,'bodyInNoExt',[],'surfInNoExt',fileName,'dataInNoExt',fileName);
        rec=reconAlgorithm(rec); % default parameter
        
        %%% READ HEADER AND BUFFER
        rec = ReadIsmrmrdBufferHdr(rec,bufferData, hdrConnection);
        bufferData_header.bits.buffer.headers = bufferData.bits.buffer.headers; % used for sending out the data
        bufferData_header.bits.buffer.sampling_description = bufferData.bits.buffer.sampling_description;
        bufferData = []; hdrConnection = [];
        
        %%% INVERT ISMRMRD TO REC (REFSCAN)
        fprintf('=========== Invert ismrmrd to rec for refscan.\n');
        useACS_flag = 1; % ZN: currently this version only support inherent calibration lines (=1; both seperate or integrated), not for external REF yet
        rec=reconInvert_gadg(rec,2); % ZN: typ = 1, body; =2, refscan; =3, high-res
        if rec.Alg.writeAfterReconInvert; rec=reconWrite_gadg(rec,2,'reconInvert');end % write the refscan in nii
        
        %%% COIL SENSITIVTY MAP ESTIMATION
        fprintf('=========== Estimate coil sensitivity map.\n');
        rec=reconCoilSens(rec);
        delete(gcp('nocreate')); % ZN: free the cpu cores for other calculation
        if rec.Alg.writeCoilSens; rec=reconWrite_gadg(rec,2,'reconCoilSens');end % write the refscan coil sensitivity map in nii
        
        %%% INVERT ISMRMRD TO REC (HIGHRES)
        rec=reconInvert_gadg(rec,3); % ZN: typ = 1, body; =2, refscan; =3, high-res
        if rec.Alg.writeAfterReconInvert; rec=reconWrite_gadg(rec,3,'reconInvert');end % write the refscan in nii
        
        %%% SENSE RECONSTRUCTION
        rec=reconSENSE_gadg(rec);
        rec=reconWrite_gadg(rec,3,'final'); 
        x = squeeze(gather(rec.x));% rec.x will be the reconstructed image

        %%% DEAL WITH ORIENTATION & SPLITE ECHOES 
        % ZN: currently only work for dual-echo SWI
        x_1 = x(:,:,:,1);
        x_2 = x(:,:,:,2);
        x_1 = flipPermute(x_1, [1 1 1]); % need to be double check
        x_1 = circshift(x_1, [1 1 1]);
        x_2 = flipPermute(x_2, [1 1 1]);
        x_2 = circshift(x_2, [1 1 1]);

        %%% CREATE FILE TO SEND
        % ZN: just for dual-echo MEGE
        %Phase
        % Magnitude
        fileToSend_1 = ...
            gadgetron.types.Image.from_data(...
                                            permute(abs(x_1), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                            referenceFromReconData(bufferData_header) ...
                                            );
        %fileToSend.header.image_type = gadgetron.types.Image.COMPLEX;
        fileToSend_1.header.image_type = gadgetron.types.Image.MAGNITUDE;
        fileToSend_1.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend_1.header.image_series_index = 1;
        fileToSend_1.header.image_index = 1;

        fileToSend_2 = ...
            gadgetron.types.Image.from_data(...
                                            permute(abs(x_2), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                            referenceFromReconData(bufferData_header) ...
                                            );
        %fileToSend.header.image_type = gadgetron.types.Image.COMPLEX;
        fileToSend_2.header.image_type = gadgetron.types.Image.MAGNITUDE;
        fileToSend_2.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend_2.header.image_series_index = 1;
        fileToSend_2.header.image_index = 2;


        %Phase
        fileToSend_3 = ...
            gadgetron.types.Image.from_data(...
                                            permute(angle(x_1), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                            referenceFromReconData(bufferData_header) ...
                                            );
        fileToSend_3.header.image_type = gadgetron.types.Image.PHASE;
        fileToSend_3.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend_3.header.image_series_index = 2;
        fileToSend_3.header.image_index = 1;

        fileToSend_4 = ...
            gadgetron.types.Image.from_data(...
                                            permute(angle(x_2), [4 1:3]), ...%Coil profiles first by convention of this (gadgetron) function
                                            referenceFromReconData(bufferData_header) ...
                                            );
        fileToSend_4.header.image_type = gadgetron.types.Image.PHASE;
        fileToSend_4.header.data_type = gadgetron.types.Image.FLOAT;
        fileToSend_4.header.image_series_index = 2;
        fileToSend_4.header.image_index = 2;

        % save for delayed pick
        save_path = rec.Nam.caseIn;
        image_name = ['image_to_be_send_highres_',fileName,'.mat'];
        nameRef = fullfile(save_path, image_name); % ZN: this version is saved for backup, but not read by auto dummy scan (can be used for callback dummy by manually selecting this file)
        save(nameRef, 'fileToSend_1','fileToSend_2','fileToSend_3','fileToSend_4','-v7.3');
        nameRef = fullfile(save_path, 'image_to_be_send_highres_SWI.mat'); % ZN: this version is saved for auto dummy scan, but will be overwritten and save the last one with the same seq type
        save(nameRef, 'fileToSend_1','fileToSend_2','fileToSend_3','fileToSend_4','-v7.3');


        fprintf('\n')
        currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
        fprintf([currentTime,'\n']);
        fprintf('===========================     THE END [MEGE]   ==============================\n');
        fprintf('===============================================================================\n');
        fprintf('\n\n\n')
        % stop diary
        diary('off');

%         % ZN: need to disable this timeline logging for the multi-GPU
%         % version, to prevent it abort the computation of other sequences
%         %%% PRINT THE TIME LINE (FOR ABORTION SECURITY)
%         logFilePath = fullfile(pwd, 'log', 'record_Disorder_timeline.txt');
%         fileID = fopen(logFilePath, 'a');
%         if fileID ~= -1
%             currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
%             fprintf(fileID, [currentTime,'  ']);
%             fprintf(fileID, 'Close the MATLAB handle: disorder_subcall_bash. Sequence&Patient:%s&%s\n',fileName,patientID);
%             fclose(fileID);
%         end

        %%% UPDATE THE GPU LOG FILE (TO RELEASE THE CURRENT GPU)
        if enable_multiGPU; updateGPUStatus(gpuIndex_activate, 'free'); end
    catch SubcallError
        fprintf('An error occurred in the subcall function!\n');
        if enable_multiGPU; logError(SubcallError,gpuIndex_activate);  end% ZN: catch the error log & log off the current GPU & exit
    end
    
    % exit
    exit;
end

% complexToFloat
% extractGadget

