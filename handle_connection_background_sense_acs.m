function handle_connection_background_sense_acs(Recon_ID)
    
    if nargin < 1 || isempty(Recon_ID); Recon_ID='SENSE_ACS';end % modify this if needed
    
    %% Default settings
    % Fix the X11 related warning and crashes
    opengl('save', 'none'); % stop showing X11 related warning
    
    % Fix for CUBLAS issues
    warning off parallel:gpu:device:DeviceLibsNeedsRecompiling
    try
        gpuArray.eye(2)^2;
    catch; ME
    end
    
    % reset the gpu memory and reset gpu to solve the memory issue 
    gpuArray([]);
    
    %% Write the log
    currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
    timestamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    diary_filename = sprintf('%s/log/Log_%s_%s.txt', pwd, Recon_ID, timestamp);
    diary(diary_filename); 
    fprintf([currentTime,'\n']);
    fprintf('===================== Getting into handle to call custom recon ====================\n');
    fprintf('===================== Performing SENSE reconstruction with ACS ====================\n');
    fprintf('Current folder:\n');
    disp(pwd);
    
    %% Add path
    addpath(genpath('Gadgetron_tools')) 
    addSupportPath(Recon_ID); % Add the path of the custom recon (the path set in Framework_config.xml in the field <support_func_path>)
    
    %% Handling multi-GPU status for parallel computation
    enable_multiGPU = 1; 
    gpuInfo = checkGPUAvailability(); % check GPU resources -> sometimes GPU driver update can cause disconnection
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
        error('No GPU accessable...')
    end
    
    try
        %% Load the saved twix-like raw by the sorting info
        fprintf('=========== Loading the raw data.\n');
        [fileName, nameRef, patientID, save_path, save_time, sorting_info] = read_sortinginfo(Recon_ID); % ZN: the earliest registered data will be loaded for processing
        pause(5); load(nameRef); % 'bufferData','fileName','hdrConnection' are saved in nameRef
        flag_erase = erase_sortinginfo(sorting_info, save_path,patientID,fileName, save_time, Recon_ID); % ZN: the loaded data's sorting info will be wiped from the structure to prevent re-recon; corresponding main_to_sub file is updated
        if flag_erase
            fprintf('=========== Raw data successfully loaded.\n');
        else
            fprintf('!!!WARNING: ERROR IN RAW DATA READING AND SORTING INFO UPDATING!!!')
        end
        
        
        %% Custom reconstruction
        % Put your custom reconstruction here
        
        % e.g., the SENSE reconstruction demo case
        rec=InlineReconPipeline(twix_like);
        
    catch SubcallError % If any errors during the custom recon
        fprintf('An error occurred in the custom recon!\n');
        if enable_multiGPU; logError(SubcallError,gpuIndex_activate);  end% ZN: catch the error log & log off the current GPU & exit (put 0 as the last input if do not wish to exit)
    end
    
    exit; % exit the matlab program to save the RAM
    
end