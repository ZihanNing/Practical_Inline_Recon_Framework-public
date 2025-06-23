
function [] = handle_connection_retrieval_template(connection)
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
    
    %%%%%%%%%%%%%%% HERE NEED TO BE MODIFIED %%%%%%%%%%%%%
    Recon_ID = 'SENSE_exterREF';
    debug_mode = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    %% LOAD THE IMAGE TO BE RETRIEVED
    hdr = connection.header;
    if ~debug_mode
        % For auto mode, the retrieved image is loaded by sorting the
        % Recon_ID and the ID of the retrieval scan/retro-recon
        % The image with the same patient ID & same sequence name/seq_type
        % is loaded automatically 
        % Caution: if you are using retrieval dummy scan, please make sure
        % that the name aligned with the target scan
        
        % sequence name
        fileName = hdr.measurementInformation.protocolName;
        % patient ID
        if isfield(hdr,'subjectInformation') && isfield(hdr.subjectInformation,'patientID')
            patientID = hdr.subjectInformation.patientID;
        else
            fprintf('Fetal error: Patient ID cannot recognize.\n');
            fprintf('=========== Exiting without successful retrieval.\n');
            exit;
        end
        % seq_typ: just in case that the retrieval dummy scan fails to match
        % the name
        [tScanningSequence,~] = searchUPfield(hdr.userParameters,'String','tScanningSequence_zihan'); % sequency type
        [seq_type,~] = recogSeqType(hdr.sequenceParameters.sequence_type,...
            tScanningSequence,1);

        % load the target image to be retrieved
        [fileToSend,status] = load_target_image(Recon_ID, patientID, fileName, seq_typ);
    else
        % In the debug mode, you manually load the images to be retrieved
        % Caution: you should launch matlab program manually and let it to
        % listen to the port (in default, port 23200) as well
        
        fprintf('Manual step: please load the image to be retrieved manually! \n');
        pause; 
    end
    
    %% RETRIEVE THE IMAGE BACK TO CONSOLE
    if status % the image has been found
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
    else
        fprintf('=========== No images found or failed to load. Exit for now...\n');
        exit;
    end
end


