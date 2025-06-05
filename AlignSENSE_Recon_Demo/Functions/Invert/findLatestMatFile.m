function fileName = findLatestMatFile(SeqName, load_path, seq_type)
%%%%%%%%%%%%5
% usually the dummy scan for image retrieval should be named as 
% high_res_Dummy -> so in this way, the correct seq can be identified
% however, if you mis-named the dummy sequence
% * without Dummy in the name, that's fine, its name should be identical to the
% sequence to be retrieve, so that the right one retrieved -> but good to
% keep the dummy in the name to distinguish the two (if scanner-recon is
% also enabled)
%
% * the name of the dummy is not identical of the original seq; the image
% saved with the seq type in the name will be retreived then, but to note,
% this file could get overwritten if several sequences with the same type
% is implemented and scanned in a single examinatino
%
% * high-res scans scanned twice without changing the name, the one with
% the later timestamp will be retrieved then


    % Check if the SeqName contains 'Dummy' or 'dummy'
    if contains(SeqName, 'Dummy', 'IgnoreCase', true)
        extract_SeqName = regexprep(SeqName, '[_]*Dummy[_]*', ''); % Remove 'Dummy' and surrounding underscores
    else
        extract_SeqName = SeqName;
    end
    
    % Check if the load_path exists
    if ~isfolder(load_path)
        error('The specified load_path does not exist.');
    end
    
    % Construct the search pattern for files
    searchPattern = fullfile(load_path, ['image_to_be_send_highres_' extract_SeqName '*.mat']);
    
    % Get all matching files in the directory
    matFiles = dir(searchPattern);
    
    % Case A: Exactly one file found
    if numel(matFiles) == 1
        fileName = matFiles(1).name;
        return;
    end
    
    % Case B: More than one file found, find the latest based on timestamp
    if numel(matFiles) > 1
        % Extract timestamps from filenames and compare them
        timestamps = [];
        for i = 1:numel(matFiles)
            % Assuming timestamp format is hh-mm-ss within the filename
            [~, fileName, ~] = fileparts(matFiles(i).name);
            timestampStr = regexp(fileName, '\d{2}-\d{2}-\d{2}', 'match', 'once'); % Get timestamp part
            if ~isempty(timestampStr)
                timestamps = [timestamps; datetime(timestampStr, 'InputFormat', 'HH-mm-ss')];
            else
                timestamps = [timestamps; NaT]; % In case no timestamp is found
            end
        end
        
        % Find the index of the latest timestamp
        [~, idx] = max(timestamps);
        fileName = matFiles(idx).name;
        return;
    end
    
    % Case C: No files meet the requirement
    seq_file = ['image_to_be_send_highres_',seq_type,'.mat'];
    fallbackFile = seq_file;
    if isfile(fallbackFile)
        fileName = fallbackFile;
        return;
    end
    
    % If no suitable file is found, report error
    error('No suitable .mat file found in the specified load_path.');
end
