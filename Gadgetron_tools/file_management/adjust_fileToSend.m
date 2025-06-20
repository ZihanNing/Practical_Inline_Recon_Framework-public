function [fileToSend, status] = adjust_fileToSend(fileToSend, target_matrixSize, target_fov)
% Adjust fileToSend.Mag and fileToSend.Phs to match target matrix size and FOV
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 20-Jun-2025

% Helper: determine which substructure exists
if isfield(fileToSend, 'Mag') && ~isempty(fileToSend.Mag)
    refData = squeeze(fileToSend.Mag{1}.data);
elseif isfield(fileToSend, 'Phs') && ~isempty(fileToSend.Phs)
    refData = squeeze(fileToSend.Phs{1}.data);
else
    error('Neither Mag nor Phs exists in fileToSend.');
end

actual_matrixSize = size(refData); % Should be [RO, Lin, Par]

if isequal(actual_matrixSize(:), target_matrixSize(:))
    fprintf('Matrix size matched - ready for retrieval. \n');
    status = 'identical';
    return;
else
    fprintf('Mismatch between the matrix size of the retrieval scan and the image to be retrieved detected. \n');
    fprintf('Retrieval scan: [%s,%s,%s] \n',num2str(target_matrixSize(1)),...
        num2str(target_matrixSize(2)),num2str(target_matrixSize(3)));
    fprintf('Image to be retrieved: [%s,%s,%s] \n',num2str(actual_matrixSize(1)),...
        num2str(actual_matrixSize(2)),num2str(actual_matrixSize(3)));
end

% Define handlers
fieldsToAdjust = {'Mag', 'Phs'};
for i = 1:numel(fieldsToAdjust)
    fieldName = fieldsToAdjust{i};
    if isfield(fileToSend, fieldName)
        for j = 1:numel(fileToSend.(fieldName))
            % Original data: [1, RO, Lin, Par]
            origData = squeeze(fileToSend.(fieldName){j}.data);
            origSize = size(origData);

            % Interpolate or crop
            if all(target_matrixSize >= actual_matrixSize)
                % Interpolation to target size
                interpData = interpn(double(origData), ...
                    linspace(1, origSize(1), target_matrixSize(1)), ...
                    linspace(1, origSize(2), target_matrixSize(2)), ...
                    linspace(1, origSize(3), target_matrixSize(3)), 'linear', 0);
                status = 'interpolate';
            else
                % Cropping
                startIdx = floor((actual_matrixSize - target_matrixSize)/2) + 1;
                endIdx = startIdx + target_matrixSize - 1;
                % Clamp indices to actual size
                for d = 1:3
                    if startIdx(d) < 1, startIdx(d) = 1; end
                    if endIdx(d) > actual_matrixSize(d), endIdx(d) = actual_matrixSize(d); end
                end
                interpData = origData(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3));
                status = 'cropped';
            end

            % Reshape back to [1 RO Lin Par]
            fileToSend.(fieldName){j}.data = reshape(interpData, [1, size(interpData)]);
            fileToSend.(fieldName){j}.header.matrix_size = target_matrixSize(:)';
            fileToSend.(fieldName){j}.header.field_of_view = target_fov(:)';
        end
    end
end

switch status
    case 'interpolate'
        fprintf('Interpolated the image to be retrieved to the target matrix size. \n');
    case 'cropped'
        fprintf('Cropped the image to be retrieved to the target matrix size. \n');
end

end