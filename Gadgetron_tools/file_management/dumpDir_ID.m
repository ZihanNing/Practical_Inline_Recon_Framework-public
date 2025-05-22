function dirName = dumpDir_ID(patientID, ProjectName)
% dumpDir_ID Creates a data directory based on patientID and component configuration
%   dirName = dumpDir_ID(patientID, componentName)
%   - patientID (optional): string or number representing the patient ID
%   - componentName (optional): name of the component (e.g., 'Default', 'Sodium_radial')

% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 22-May-2025

if nargin < 1; patientID = []; end
if nargin < 2; ProjectName = 'Default'; end

% Get current working directory
current_home_path = pwd;

% Search for Framework_config.xml
configPath = findConfigXML(current_home_path, 'Framework_config.xml');

% Default data path
default_data_path = fullfile(current_home_path, 'data');

if ~isempty(configPath)
    try
        xDoc = xmlread(configPath);
        components = xDoc.getElementsByTagName('Component');
        data_path = default_data_path;  % fallback

        for i = 0:components.getLength-1
            comp = components.item(i);
            nameAttr = char(comp.getAttribute('name'));
            if strcmp(nameAttr, ProjectName)
                data_path = char(comp.getElementsByTagName('data_saved_path').item(0).getTextContent());
                break;
            end
        end
    catch
        warning('Failed to parse config. Using default data path.');
        data_path = default_data_path;
    end
else
    data_path = default_data_path;
end

% Create directory path based on whether patientID is given
dateName = datestr(datetime('today'), 'dd-mmm-yyyy');
if nargin < 1 || isempty(patientID)
    dirName = fullfile(data_path, dateName);
else
    dirName = fullfile(data_path, [dateName '-P' num2str(patientID)]);
end

% Create directory if needed
if ~isfolder(dirName)
    mkdir(dirName);
end
end

function configPath = findConfigXML(baseDir, fileName)
% Recursively search for fileName starting at baseDir
files = dir(fullfile(baseDir, '**', fileName));
if ~isempty(files)
    configPath = fullfile(files(1).folder, files(1).name);
else
    configPath = '';
end
end
