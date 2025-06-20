
function [] = handle_connection_retrieval_manual(connection)
    %%% This is a handler to retrieve the reconstructed images back to
    %%% console 
    %%% Two methods for retrieval available based on this framework:
    %%% - short retrieval dummy scan: its name should contains the target
    %%% scan's name and the planning should be the same as the target scan
    %%% - retro recon
    %%%
    %%% by Zihan Ning
    
    %% SET ENVIRONMENT
    % ZN: better to be sure that the matlab is launched in
    % /Gadgetron_Parallel_Framework folder to allow file/function
    % recognition    
    addpath(genpath('Gadgetron_tools'))
    Recon_ID = 'Manual_retrieval';
    
    %% Write the log
    currentTime = datestr(datetime('now','Format','yyyy-MM-ddd HH:mm:ss:ss.SSS'));
    timestamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    diary_filename = sprintf('%s/log/Log_retrieval_%s_%s.txt', pwd, Recon_ID, timestamp);
    diary(diary_filename); 
    fprintf([currentTime,'\n']);
    fprintf('======================== Retrieving images to the MR console ======================\n');
    disp(pwd);
    
    %% READ CONNECTION
    fprintf('=========== Reading the connection...\n');
    bufferData = connection.next();% for the retrieval scan, we simply use the default gadget from BucketToBuffer
    
    %% LOAD THE IMAGE MANUALLY
    
    %%%%%%%%%%%%%%%%% SET THE BREAKPOINT %%%%%%%%%%%%%%%%%
    % Set a breakpoint here and to load the image_to_be_send file manually
    pause;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% RETRIEVE THE IMAGE BACK TO CONSOLE
    hdr = connection.header;
    % double check the matrix size of the image to be sent and the
    % retrieval dummy scan/retro recon alignes
    fprintf('=========== Checking the alignment of matrix size between retrieval scan and the image.\n');
    target_matrixSize = hdr.encoding.reconSpace.matrixSize;
    target_fov = hdr.encoding.reconSpace.field_of_view;
    [fileToSend, status] = adjust_fileToSend(fileToSend, target_matrixSize, target_fov);

    % send the images
    fprintf('=========== Sending out the images one by one to the console...\n');
    fields = fieldnames(fileToSend);

    for i = 1:numel(fields)
        currentField = fields{i};
        substructures = fileToSend.(currentField);

        if iscell(substructures)
            for j = 1:numel(substructures)
                connection.send(substructures{j});
            end
        end
    end

    % Finish and closing
    fprintf('======================== Images all retreived to the console ======================\n');
    diary_Flag = 1;
    if diary_Flag == 1
        diary off;
    end
    exit;
end


