function new_hdrConnection = check_overwrite(hdrConnection)
% ZN: to prevent overwriting, the sequence name of each data will be
% modified as 'seq_name_timestamp'

% modify the hdrConnection for a new seq name to prevent overwritten
fileName = hdrConnection.measurementInformation.protocolName;
currentTime = datestr(now, 'HH-MM-SS');
new_fileName = [fileName,'_',currentTime]; % ZN: new file name with timestamp after it
new_hdrConnection = hdrConnection;
new_hdrConnection.measurementInformation.protocolName = new_fileName; % replace the seq name
fprintf('To prevent overwriting, the seq name has been modified from %s to %s! \n',fileName, new_fileName);
