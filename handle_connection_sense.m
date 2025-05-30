
function [] = handle_connection_sense(connection)
    %%% This is a blank handle
    %%% Just read in Bucket including noise scan (refer to:
    %%% config/configDISORDER_test_Bucket_SWI.xml)
    %%% by Zihan Ning
    
    %% SET ENVIRONMENT
    % ZN: better to be sure that the matlab is launched in
    % /Gadgetron_Parallel_Framework folder to allow file/function
    % recognition
    
    addpath(genpath('Gadgetron_tools'))
    
    %% READ CONNECTION
    fprintf('=========== Reading ismrmrd connection in Bucket.\n');
    
    noise_bucket = {};
    ReadNext_flag = 1;
    loop_counter = 0;
    while ReadNext_flag
        bucketData = connection.next();
        loop_counter = loop_counter + 1;
        if bucketData.data.count == 1 
            % ZN: to make this work, noise lines should awalys be read
            % first and usually line-by-line
            % afterwards, the image & ref data will be read as a whole,
            % which makes bucketData.data.count > 1
            noise_bucket{loop_counter} = bucketData.data.data; % ZN: accumulate noise data
        else % ZN: this is the image & ref bucket, which should be read at the end
            fprintf(':: Finished read noise bucket: %s lines.\n',num2str(loop_counter - 1));
            fprintf(':: Start to read image and ref bucket.\n');
            % current bucketData is the data with kspace lines for imaging
            % and reference lines
            ReadNext_flag = 0; % ZN: this should be the last one to be read, stop looping
        end    
    end
    
    % convert noise lines to 3d matrix
    if ~isempty(noise_bucket) && loop_counter>1
        % ZN: we assume that except the last readin image&ref bucket, the
        % previous lines are all noise lines
         for i = 1:loop_counter-1
             noiseData(:,:,i) = noise_bucket{i};
         end
         fprintf(':: Accumulated noise scan: [%s,%s,%s] \n',num2str(size(noiseData,1)),num2str(size(noiseData,2)),num2str(size(noiseData,3)));
    else
        noiseData = [];
        fprintf('There is no noise lines readin - double check the data or config file! \n');
    end
    
    % ZN: for debugging & developing, save the essential data to prevent
    % resending the kspace data
%     save('radial_test_bucket.mat','connection','bucketData','noiseData','-v7.3')

    %% GENERAL INPUT CONVERTOR (ISMRMRD > TWIX-LIKE RAW)
    % bucket2buffer | matlab version
    twix_like = BucketToBuffer_matlab(bucketData,connection.header,noiseData,0);
    
    %% SAVE RAW AND INFO FOR PARALLEL COMPUTATION
    Recon_ID = 'SENSE_ACS'; % the configurations of each research ID is set within Framework_config.xml, including saved path, related bash, etc
    save_raw_info(twix_like,connection.header,Recon_ID,Recon_ID);
    
    %% CALL SUBCALL FUNCTION IN PARALLEL AND RELEASE SCANNER RECON CHAIN
    debug_mode = 1; % 1: debug mode-run the background function on the interface directly; 0: auto mode-launch a matlab program to run target function on the background
    if debug_mode
%         handle_connection_background_sense_acs
        fprintf('=========== Debug mode \n');
        feval(getReconAlg(Recon_ID)) % the implemented recon pipe can be attached with Recon_ID in the config
    else % call the function by a bash 
        cfg = readScriptConfig(Recon_ID); % get the bash script input by the Recon_ID (info in Framework_config.xml)
        % excute the bash script in parallel
        bashScript = './Bash_script/Parallel_launch_MonitorGPU.sh'; % The bash script 'Parallel_launch_MonitorGPU' will monitor GPU usage status and decide to launch recons or queue up
        cmd = sprintf('%s -l "%s" -f "%s" -r "%s" -w %d -i %d &', ...
              bashScript, ...
              cfg.logFile, ...
              cfg.matlabFunction, ...
              Recon_ID, ...
              cfg.maxWaitTime, ...
              cfg.interval);
        fprintf('=========== Launching custom recon %s in the background... \n', cfg.matlabFunction);
        status = system(cmd);
        if status~=0
            error('%s:LaunchFailed', ...
                  'Failed to start gpu_wait.sh (exit code %s).', bashScript, num2str(status));
        end 
    end
    
    %% REALEASE THE ICE PIPE
    fprintf('=========== Close the saving procedure and return control to ICE. \n');    
    if ~debug_mode; exit; end
end


