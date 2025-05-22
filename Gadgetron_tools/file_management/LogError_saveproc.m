function LogError_saveproc(SubcallError, ProjectName)
% LogError_saveproc Saves detailed error log to a timestamped text file
%   LogError_saveproc(SubcallError, ProjectName)
%   - SubcallError: caught error object from try-catch
%   - ProjectName: string identifier of the current project

% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 22-May-2025

% Default log folder
log_path = fullfile(pwd, 'log');
if ~isfolder(log_path)
    mkdir(log_path);
end

% Generate timestamp and filename
timestamp = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
filename = ['ErrorDuringSaving_' ProjectName '_' timestamp '.txt'];
filepath = fullfile(log_path, filename);

% Open and write log file
fid = fopen(filepath, 'w');
if fid == -1
    warning('Unable to open log file: %s', filepath);
    return;
end

fprintf(fid, 'Error occurred at: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Error message: %s\n', SubcallError.message);
fprintf(fid, 'Error identifier: %s\n', SubcallError.identifier);
fprintf(fid, 'Stack trace:\n');
for k = 1:length(SubcallError.stack)
    fprintf(fid, 'File: %s, Line: %d\n', ...
        SubcallError.stack(k).file, SubcallError.stack(k).line);
end
fprintf(fid, '-----------------------\n');

fclose(fid);
end
