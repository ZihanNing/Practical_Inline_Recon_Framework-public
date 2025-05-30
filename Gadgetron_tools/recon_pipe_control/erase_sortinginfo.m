function flag_erase = erase_sortinginfo(sorting_info, save_path,patientID,fileName, save_time, seq_type)
% ZN: this function is used for erasing the sorting info of case that has
% already complete processing from the structure

if nargin<6 || isempty(seq_type);seq_type = [];end

% Define the caseToDelete
caseToDelete.patientID = patientID;
caseToDelete.fileName = fileName;
caseToDelete.save_path = save_path;
caseToDelete.save_time = save_time;

% Find the index of the subject(s) that match all three criteria
idx = find(strcmp({sorting_info.patient_ID}, caseToDelete.patientID) & ...
           strcmp({sorting_info.fileName}, caseToDelete.fileName) & ...
           strcmp({sorting_info.save_path}, caseToDelete.save_path) & ...
           strcmp({sorting_info.save_time}, caseToDelete.save_time));

% If any subject matches the criteria, erase it
if ~isempty(idx)
    sorting_info(idx) = [];  % Delete the matching subject(s)
    flag_erase = 1;
    fprintf('Sequence %s from subject %s has been erased from sorting_info since reconstruction complete! \n', fileName, patientID);
else
    flag_erase = 0; 
    fprintf('No group meets criteria for erasing from sorting_info!!! Same sequence with name %s of subject %s might be re-reconstructed! \n', fileName, patientID);
end

% update the sorting_info file
sorting_fileName = ['main_to_sub_',seq_type,'.mat'];
save(sorting_fileName,'sorting_info');
