function [fileToSend_all, status] = load_target_image(Recon_ID, patientID, fileName, seq_type)
% Load image data based on reconstruction ID, patient ID, file name, and sequence type
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 20-Jun-2025

status = 0; % 0-fail; 1-success

% Step 1: Parse XML to find data_saved_path
xmlPath = fullfile(pwd, 'Framework_config.xml');
if ~isfile(xmlPath)
    error('Framework_config.xml not found in the current directory.');
end

try
    xDoc = xmlread(xmlPath);
    allComponents = xDoc.getElementsByTagName('Component');
    foundPath = '';
    for i = 0:allComponents.getLength-1
        comp = allComponents.item(i);
        if strcmp(char(comp.getAttribute('name')), Recon_ID)
            pathNode = comp.getElementsByTagName('data_saved_path').item(0);
            foundPath = char(pathNode.getFirstChild.getData);
            break;
        end
    end
catch
    error('Error parsing XML.');
end

if isempty(foundPath)
    error('No matching Recon_ID (%s) found in Framework_config.xml.', Recon_ID);
end

% Step 2: Search subfolder with patientID and today's date
folders = dir(foundPath);
folders = folders([folders.isdir]);
todayStr = datestr(datetime('today'), 'dd-mmm-yyyy');
latestDate = datetime(0,0,0);
target_folder = '';

for i = 1:length(folders)
    fname = folders(i).name;
    if strcmp(fname, '.') || strcmp(fname, '..')
        continue;
    end
    if contains(fname, patientID)
        % Try extracting the date from folder name
        tokens = regexp(fname, '^(\d{2}-[A-Za-z]{3}-\d{4})', 'tokens');
        if ~isempty(tokens)
            fdate = datetime(tokens{1}{1}, 'InputFormat', 'dd-MMM-yyyy');
            if strcmp(tokens{1}{1}, todayStr)
                target_folder = fullfile(foundPath, fname);
                break;
            elseif fdate > latestDate
                latestDate = fdate;
                target_folder = fullfile(foundPath, fname);
            end
        end
    end
end

if isempty(target_folder)
    error('No folder found for patient ID %s.', patientID);
elseif ~contains(target_folder, todayStr)
    fprintf('No folder found for today. Using latest folder: %s\n', target_folder);
end

% Step 3: Search for image_to_be_send file
files = dir(fullfile(target_folder, sprintf('image_to_be_send*%s*.mat', fileName)));

% Sort by timestamp if multiple
if length(files) > 1
    [~, idx] = sort([files.datenum], 'descend');
    chosenFile = files(idx(1));
elseif length(files) == 1
    chosenFile = files(1);
else
    % Try with seq_type
    files = dir(fullfile(target_folder, sprintf('image_to_be_send*%s*.mat', seq_type)));
    if isempty(files)
        error('No image_to_be_send file found with fileName or seq_type.');
    else
        [~, idx] = sort([files.datenum], 'descend');
        chosenFile = files(idx(1));
    end
end

fprintf('Loading image data from: %s\n', fullfile(chosenFile.folder, chosenFile.name));
try
    load(fullfile(chosenFile.folder, chosenFile.name), 'fileToSend_all');
    if exist('fileToSend_all', 'var')
        status = 1;
        fprintf('=========== Image to be retrieved loaded: \n');
        fprintf('- Recon ID: %s\n', Recon_ID);
        fprintf('- Patient ID: %s\n', patientID);
        fprintf('- Image reconstructed time: %s\n', chosenFile.date);
        fprintf('- Loaded file name: %s\n', chosenFile.name);
        fprintf('- Loaded file path: %s\n', chosenFile.folder);    
    else
        fprintf('Target parameter fileToSend_all not exist in the target file. \n');
    end
catch
    fprintf('Error when loading the target file: %s. \n',fullfile(chosenFile.folder, chosenFile.name));
end

end
