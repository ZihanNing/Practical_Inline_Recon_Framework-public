function sorting_info = generate_sortinginfo(nameRef,save_path,patientID,fileName, save_tag)
% ZN: this function is to generate the main_to_sub matlab file 
% this file is the link between maincall and subcall function, it will let
% subcall function know the information of the to-be-processed data

% fileName: sequence name
% nameRef: the save path of the raw data
% patient ID: subject ID number
% save_path: the saving folder path for this data
% save_time: the timestamp of generating this data group for the current subject & date & sequence

%%% READ THE FILE OR GENERATING A NEW STRUCTURE
sorting_fileName = ['main_to_sub_', save_tag, '.mat'];
if exist(sorting_fileName, 'file') == 2
    % If the file exists, load it
    load(sorting_fileName);
    disp(['File ', sorting_fileName, ' loaded successfully.']);
else
    % If the file does not exist, do nothing
    disp(['File ', sorting_fileName, ' does not exist - generating a new structure']);
    sorting_info = struct('fileName', {}, 'nameRef', {}, 'patient_ID', {}, 'save_path', {}, 'save_time', {});
end

%%% ADDING A NEW GROUP (SUBJECT & DATE & SEQ)
% Define the data for one subject
newCase.fileName = fileName;
newCase.nameRef = nameRef;
newCase.patient_ID = patientID;
newCase.save_path = save_path;
newCase.save_time = datestr(now, 'yyyy-mm-dd-HH-MM-SS'); % generate the current time
% Add the new subject's data to the structure array
sorting_info(end+1) = newCase;


