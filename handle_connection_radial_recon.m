
function [] = handle_connection_radial_recon(connection)
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
    twix_like = BucketToBuffer_matlab(bucketData,connection.header,noiseData,0,'radial');
    
    %% SAVE RAW AND INFO FOR PARALLEL COMPUTATION
    save_raw_info(twix_like,connection.header,'Sodium_radial','radial_MEGE_volume');
    
    %% CALL SUBCALL FUNCTION IN PARALLEL AND RELEASE SCANNER RECON CHAIN
    
    
    % do nothing and return nothing
    % ZN: please put your own reconstruction code here
    fprintf('=========== Do nothing and shut down the MATLAB.\n');
    
    % ZN: if you'd like to send the image back to the scanner with this
    % handle, the sending out part still need to be developed
    
%     exit
%     system('./kill_all_matlab.sh');
%     exit;
end


