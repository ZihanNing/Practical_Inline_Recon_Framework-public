function [recS,status] = loadExteralREF(folderPath,nameREF)

recS = [];
status = 'NotFound';
% nameREF = '*REF_DATA_EXTRAREF_*.mat';
if isfolder(folderPath)
    % Step 2: Search for .mat file with 'REF_DATA_EXTRAREF_' in the folder
    matFiles = dir(fullfile(folderPath, nameREF));
    
    if ~isempty(matFiles)
        % If .mat files are found, find the one with the latest modification time
        [~, idx] = max([matFiles.datenum]);
        latestFile = matFiles(idx).name;
        
        % Step 3: Load the latest .mat file
        recS = load(fullfile(folderPath, latestFile));
        status = 'Found';
        fprintf('External REF saved previously has been found and will be used for following reconstruction!\n');
        fprintf('Folder: %s \n',folderPath);
        fprintf('File: %s \n',latestFile);
    else
        % Step 4: No matching .mat file found
        disp('No matching .mat file found currently');
    end
else
    % Step A: Folder does not exist, check for 20 attempts to detect folder and file
    count = 0;
    while count < 20
        if isfolder(folderPath)
            matFiles = dir(fullfile(folderPath, nameREF));
            if ~isempty(matFiles)
                [~, idx] = max([matFiles.datenum]);
                latestFile = matFiles(idx).name;
                recS = load(fullfile(folderPath, latestFile));
                status = 'Found';
                fprintf('External REF saved previously has been found and will be used for following reconstruction!\n');
                fprintf('Folder: %s \n',folderPath);
                fprintf('File: %s \n',latestFile);
                break;
            end
        end
        count = count + 1;
        pause(30); % Wait 30 second before retrying
    end
    
    if count == 20
        disp('Folder or file not detected after 20 attempts.');
    end
end
if ~isempty(recS); temp = recS.recS; clear recS; recS = temp; clear temp; end
