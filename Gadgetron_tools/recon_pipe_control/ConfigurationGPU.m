function gpuIndex = ConfigurationGPU()
    % Define the log file name
    logFileName = [pwd,'/log/GPU_status_log.txt'];
    
    % Check if the file exists
    if ~isfile(logFileName)
        error('The file %s does not exist.', logFileName);
    end
    
    % Open the file
    fileID = fopen(logFileName, 'r');
    
    % Initialize a variable to store the last line
    lastLine = '';
    
    % Read the file line by line to get the last line
    % ZN: there need to be an initial status in the log file
    while ~feof(fileID)
        currentLine = fgetl(fileID);
        if ischar(currentLine)
            lastLine = currentLine;
        end
    end
    
    % Close the file
    fclose(fileID);
    
    % Split the last line into parts based on ';' separator
    gpuStatusParts = strsplit(lastLine, ';');
    
    % Initialize an array to store free GPU indices
    freeGPUs = [];
    
    % Loop through each part and check for 'free' status
    for i = 1:length(gpuStatusParts)
        gpuStatus = strtrim(gpuStatusParts{i});
        if contains(gpuStatus, 'free')
            % Extract the GPU index
            gpuIndexStr = regexp(gpuStatus, 'GPU(\d+)_free', 'tokens');
            if ~isempty(gpuIndexStr)
                gpuIndex = str2double(gpuIndexStr{1}{1});
                freeGPUs = [freeGPUs, gpuIndex]; %#ok<AGROW>
            end
        end
    end
    
    % Check if there are any free GPUs
    % ZN: this should be taken care by the bash script of activating the
    % subcall function
    if isempty(freeGPUs)
        error('No GPU is free.');
    end
    
    % Sort the free GPUs in ascending order and return the first one
    freeGPUs = sort(freeGPUs);
    gpuIndex = freeGPUs;
end
