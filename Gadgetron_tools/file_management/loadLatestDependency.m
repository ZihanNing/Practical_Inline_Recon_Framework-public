function chosenFullPath = loadLatestDependency(save_path, keyPattern)
% loadLatestDependency  Load the most recent dependency matching a key pattern
%   S = loadLatestDependency(save_path, keyPattern) scans the folder save_path 
%   for files whose names contain keyPattern. If multiple matches are found, it 
%   parses the timestamp (HH-MM-SS) at the end of each filename and picks the file 
%   with the latest timestamp. Assumes filenames end in “_<HH-MM-SS>.<ext>”.
%
%   Inputs:
%     save_path   – full path to the directory containing the dependency files
%     keyPattern  – substring that must appear in the filename (e.g. 'ExterREF_CMS')
%
%   Output:
%     S           – struct returned by MATLAB’s load() on the chosen file
%
%   Example:
%     data = loadLatestDependency('/home/gadgetron/data/sense/', 'ExterREF_CMS');
%
% by Zihan Ning @ King's
% 3 June 2025

    if nargin < 2
        error('loadLatestDependency:MissingInputs', ...
            'Both save_path and keyPattern must be provided.');
    end
    if ~isfolder(save_path)
        error('loadLatestDependency:BadPath', 'save_path must be an existing folder.');
    end

    % List all files that contain keyPattern (any extension)
    pattern = fullfile(save_path, ['*' keyPattern '*']);
    D = dir(pattern);
    if isempty(D)
        error('loadLatestDependency:NoFiles', ...
            'No files matching "*%s*" found in %s.', keyPattern, save_path);
    end

    % Preallocate timestamp array
    nFiles = numel(D);
    timestamps = nan(nFiles,1);

    for k = 1:nFiles
        [~, name, ~] = fileparts(D(k).name);
        parts = strsplit(name, '_');
        % The timestamp is expected to be the last underscore‐delimited part, e.g. "12-11-20"
        ts_str = parts{end};
        if isempty(regexp(ts_str, '^\d{2}-\d{2}-\d{2}$', 'once'))
            error('loadLatestDependency:BadFormat', ...
                'Filename "%s" does not end in a valid HH-MM-SS timestamp.', D(k).name);
        end
        % Convert "HH-MM-SS" → "HH:MM:SS" and then to duration
        ts_colon = strrep(ts_str, '-', ':');  % e.g. "12-11-20" → "12:11:20"
        t = duration(ts_colon, 'InputFormat', 'hh:mm:ss');
        timestamps(k) = hours(t)*3600 + minutes(t)*60 + seconds(t);  % total seconds
    end

    % Find the index of the maximum (latest) timestamp
    [~, idx_max] = max(timestamps);
    chosenFile = D(idx_max).name;
    chosenFullPath = fullfile(save_path, chosenFile);

    fprintf('The latest dependency to be loaded: %s\n', chosenFile);
end
