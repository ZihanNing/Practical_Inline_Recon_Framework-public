function freeGPUs = ConfigurationGPU()
    % ConfigurationGPU  Return indices of free GPUs based on the most recent non-empty line in the log.
    %
    %   freeGPUs = ConfigurationGPU() reads the file 'log/GPU_status_log.txt'
    %   (in the current working directory), finds the last non-empty line, parses
    %   it for GPU statuses, and returns a sorted vector of all GPU indices that
    %   are marked “free”. If no free GPUs are found, an error is thrown.
    %
    %   This version skips trailing blank lines by walking backward until a
    %   non-empty line is found.
    
    % Define the log file name
    logFileName = fullfile(pwd, 'log', 'GPU_status_log.txt');
    
    % Check if the file exists
    if ~isfile(logFileName)
        error('ConfigurationGPU:NoFile', 'The file %s does not exist.', logFileName);
    end
    
    % Read all lines into a cell array
    fid = fopen(logFileName, 'r');
    allLines = {};
    idx = 1;
    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        allLines{idx} = tline;  %#ok<AGROW>
        idx = idx + 1;
    end
    fclose(fid);
    
    % Find the last non-empty line
    lastLine = '';
    for k = numel(allLines):-1:1
        if ~isempty(strtrim(allLines{k}))
            lastLine = allLines{k};
            break;
        end
    end
    
    if isempty(lastLine)
        error('ConfigurationGPU:NoStatus', 'No non-empty status line found in %s.', logFileName);
    end
    
    % Split the last non-empty line into parts based on ';' separator
    gpuStatusParts = strsplit(lastLine, ';');
    
    % Initialize an array to store free GPU indices
    freeGPUs = [];
    
    % Loop through each part and check for 'free' status
    for i = 1:numel(gpuStatusParts)
        gpuStatus = strtrim(gpuStatusParts{i});
        if contains(gpuStatus, 'free')
            % Extract the GPU index (expects something like "GPU0_free" or "GPU1_free")
            tokens = regexp(gpuStatus, 'GPU(\d+)_free', 'tokens');
            if ~isempty(tokens)
                idxNum = str2double(tokens{1}{1});
                freeGPUs(end+1) = idxNum; %#ok<AGROW>
            end
        end
    end
    
    % If no GPUs are free, throw an error
    if isempty(freeGPUs)
        error('ConfigurationGPU:NoFreeGPU', 'No GPU is free based on the latest log line.');
    end
    
    % Return sorted list of free GPUs
    freeGPUs = sort(freeGPUs);
end
