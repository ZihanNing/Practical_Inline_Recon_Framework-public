function logError(errorObject,gpuIndex,exit_flag)

    if nargin<2 || isempty(gpuIndex);gpuIndex=[];end
    if nargin<3 || isempty(exit_flag);exit_flag = 1;end

    % Get current timestamp for unique log file name
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Create a new log file name with the current timestamp
    logFileName = [pwd, '/log/SubcallError_log_', timestamp, '.txt'];
    
    % Open the new log file for writing
    fid = fopen(logFileName, 'w');  % 'w' to create a new file
    if fid == -1
        error('Could not open log file for writing.');
    end
    
    % Get the current timestamp for the log entry
    logTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    % Write error information to the log file
    fprintf(fid, 'Error occurred at: %s\n', logTime);
    fprintf(fid, 'Error message: %s\n', errorObject.message);
    fprintf(fid, 'Error identifier: %s\n', errorObject.identifier);
    fprintf(fid, 'Stack trace:\n');
    
    % Loop through the stack trace and write details to the log
    for k = 1:length(errorObject.stack)
        fprintf(fid, 'File: %s, Line: %d\n', errorObject.stack(k).file, errorObject.stack(k).line);
    end
    
    fprintf(fid, '-----------------------\n');
    fclose(fid);
    
    % Display a message indicating the error has been logged
    disp(['Error has been logged in: ', logFileName]);
    
    % log off from the GPU status & exit
    if ~isempty(gpuIndex) && exit_flag; updateGPUStatus(gpuIndex, 'free'); end % ZN: change the current GPU status to 'free' and exit the matlab
    if exit_flag; exit; end
end
