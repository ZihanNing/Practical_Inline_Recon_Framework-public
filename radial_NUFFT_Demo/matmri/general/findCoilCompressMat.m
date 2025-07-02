function [A, sv] = findCoilCompressMat(data, Nkeep, doAlignSlices)
%   Find coil compression matrix using method from Zhang et al, doi: 10.1002/mrm.24267
%
%   [A, sv] = findCoilCompressMat(data, Nkeep, doAlignSlices)
%
%   data: Ncoils x NpointsPerSlice x Nslices
%   Nkeep: number of channels to keep
%   doAlignSlices: whether to align slices
%   A:    compression matrices
%   sv:   normalized singular values for each virtual coil for each slice (Ncoil x Nslice)
%
% (c) Corey Baron 2021
% 

if nargin<3 || isempty(doAlignSlices)
    doAlignSlices = 1;
end
if nargin<2 || isempty(Nkeep) || Nkeep < 1 || Nkeep > size(data,1)
    Nkeep = size(data,1);
end
if Nkeep == size(data,1)
    % No point to align slices if all coils are kept, since this will just
    % undo the compression
    doAlignSlices = 0;
end

% Find compression matrix for each slice
A = zeros(Nkeep,size(data,1),size(data,3), 'like', data);
sv = zeros(size(data,1),size(data,3), 'like', data);
for nsl=1:size(data,3)
    [U,S,~] = svd(data(:,:,nsl),'econ');
    sv(:,nsl) = diag(S)/S(1,1);
    Up = U';
    A(:,:,nsl) = Up(1:Nkeep,:);
    if doAlignSlices && nsl>1
        C = A(:,:,nsl)*A(:,:,nsl-1)';
        [U,~,V] = svd(C,'econ');
        P = V*U';
        A(:,:,nsl) = P*A(:,:,nsl);
    end
end


end