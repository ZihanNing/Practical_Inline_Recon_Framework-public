function idxCell = findMatrixByNonZero(X)
    % cropMatrixByNonZero Crops an N-dimensional matrix by removing all-zero margins.
    %
    %   cropped = cropMatrixByNonZero(X) returns a cropped version of the input
    %   matrix X, removing all leading and trailing slices in each dimension
    %   that contain only zero elements. The function works for matrices of
    %   any number of dimensions.
    %
    %   Example:
    %       % Create a 3D matrix with zero margins
    %       X = zeros(5, 5, 5);
    %       X(2:4, 2:4, 2:4) = 1;
    %       
    %       % Crop the matrix
    %       croppedX = cropMatrixByNonZero(X);
    %
    %       % Display sizes
    %       disp(size(X));       % Outputs: [5 5 5]
    %       disp(size(croppedX));% Outputs: [3 3 3]
    %
    %   Author: [Your Name]
    %   Date: [Date]
    
    % Validate input
    if ~ismatrix(X) && ndims(X) > 2 && ~isvector(X)
        % Do nothing, proceed
    elseif isvector(X)
        % Ensure consistent handling for vectors
        X = reshape(X, [], 1);
    else
        % Handle 2D matrices as is
    end
    
    % Find linear indices of non-zero elements
    linearIdx = find(X);
    
    if isempty(linearIdx)
        % If all elements are zero, return an empty matrix
        cropped = [];
        return;
    end
    
    % Get the size of the input matrix
    sz = size(X);
    dims = length(sz);
    
    % Initialize cell arrays to hold subscript indices
    subs = cell(1, dims);
    
    % Convert linear indices to subscripts
    [subs{:}] = ind2sub(sz, linearIdx);
    
    % Initialize arrays to store minimum and maximum indices for each dimension
    minIdx = zeros(1, dims);
    maxIdx = zeros(1, dims);
    
    for d = 1:dims
        currentSubs = subs{d};
        minIdx(d) = min(currentSubs);
        maxIdx(d) = max(currentSubs);
    end
    
    % Create a cell array of index ranges for cropping
    idxCell = cell(dims,1);
    for d = 1:dims
        sz = size(X,d);
        idxCell{d} = ones(1,sz);
        % pre-pad
        if minIdx(d)>1; idxCell{d}(1:minIdx(d)-1) = 0; end
        % post-pad
        if maxIdx(d) + 1<sz; idxCell{d}(maxIdx(d)+1:sz) = 0; end
    end
end
