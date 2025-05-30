function updateGPUStatus(gpuIndex, change_status)
    % Path to log file
    logFileName = fullfile(pwd, 'log', 'GPU_status_log.txt');

    % Validate change_status
    if ~ismember(change_status, {'free','occupied'})
        error('Invalid change_status. Please choose "free" or "occupied".');
    end

    % Read non‐empty lines
    fid = fopen(logFileName,'r');
    if fid==-1, error('Could not open log file.'); end
    raw = textscan(fid,'%s','Delimiter','\n'); fclose(fid);
    lines = raw{1}(~cellfun(@isempty,strtrim(raw{1})));
    if isempty(lines), error('No valid log data found.'); end
    latest = lines{end};
    disp(['Latest Log: ', latest]);

    % Split off timestamp by locating the first ": "
    sep = strfind(latest, ': ');
    if isempty(sep)
        error('Unexpected log format: no “: ” separator found.');
    end
    entriesStr = strtrim(latest(sep(1)+2:end));
    entries    = strtrim(strsplit(entriesStr, ';'));
    nGPUs      = numel(entries);

    % Validate gpuIndex
    if gpuIndex < 1 || gpuIndex > nGPUs
        error('Invalid gpuIndex. This system has %d GPU(s).', nGPUs);
    end

    % Parse entries
    GPUStatus     = cell(1,nGPUs);
    GPUTimestamps = cell(1,nGPUs);
    for k = 1:nGPUs
        tok = regexp(entries{k}, '^GPU(\d+)_(\w+)_([^\s]+)$', 'tokens', 'once');
        if isempty(tok)
            error('Failed to parse GPU entry: %s', entries{k});
        end
        idx = str2double(tok{1});
        GPUStatus{idx}     = tok{2};
        GPUTimestamps{idx} = tok{3};
    end

    % Check if update is needed
    if strcmp(GPUStatus{gpuIndex}, change_status)
        warning('GPU%d is already %s. No update required.', gpuIndex, change_status);
        return;
    end

    % Update only the target GPU
    GPUStatus{gpuIndex}     = change_status;
    GPUTimestamps{gpuIndex} = datestr(now, 'yyyy-mm-dd_HH:MM:SS');

    % Rebuild and append log line
    timestamp  = datestr(now, 'yyyy-mm-dd_HH:MM:SS');
    partsNew   = arrayfun(@(k) sprintf('GPU%d_%s_%s', k, GPUStatus{k}, GPUTimestamps{k}), ...
                          1:nGPUs, 'UniformOutput', false);
    newLine    = [timestamp, ': ', strjoin(partsNew, '; ')];
    fid = fopen(logFileName,'a');
    if fid==-1, error('Could not open log file for appending.'); end
    fprintf(fid, '%s\n', newLine); fclose(fid);

    disp('Log file updated successfully.');
    disp(['New Log: ', newLine]);
end
