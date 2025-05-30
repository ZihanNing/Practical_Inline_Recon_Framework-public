function [fileName, nameRef, patientID, save_path, save_time, sorting_info] = read_sortinginfo(save_tag)
% ZN: this function will read the sorting_info from the main_to_sub file
% for each seq type
% The earliest registered subject will be read from the structure to be
% processed

if nargin<1 || isempty(save_tag);save_tag = [];end

sorting_fileName = ['main_to_sub_',save_tag,'.mat'];
try
    load(sorting_fileName); 
catch
    fprintf('Fetal:: Cannot found the file %s contains the sorting info! \n',sorting_fileName);
end

% Checking and loading
if isempty(sorting_info)
    error('There is no registered data to be processed! \n')
else
    numData = numel(sorting_info);
    fprintf('Current number of data %s to be processed: %s. \n', save_tag, num2str(numData));
end

timeArray = zeros(1, length(sorting_info)); % ZN: load the earliest registered data
for i = 1:length(sorting_info)
    timeArray(i) = datenum(sorting_info(i).save_time, 'yyyy-mm-dd-HH-MM-SS');
end
[~, earliestIdx] = min(timeArray);

% Display the subject with the earliest save_time
earliestData = sorting_info(earliestIdx);
disp('Data with the earliest save_time will be processed:');
disp(earliestData);

fileName = earliestData.fileName;
nameRef = earliestData.nameRef;
patientID = earliestData.patient_ID;
save_path = earliestData.save_path;
save_time = earliestData.save_time;



