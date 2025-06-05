function mask = simple_mask(I,thres, isolated_region_size, enlarge_factor)

%%% This function generate mask for the head based on the root of 
%%% sum-of-square of the coil sensitivity map
%%% 
%%% I - root of sum-of-square of the coil sensitivity map, range from 0 to
%%% 1
%%% thres - the threshold for the hard mask to generate the initial mask,
%%% range from 0 to 1, defualt: 0.3
%%% isolated_region_size - to discard the small isolated region in the
%%% initial mask, defualt: 800
%%% enlarge_factor - expand the mask, defualt 15

if nargin<2;thres = 0.3;end
if nargin<3;isolated_region_size = 800;end
if nargin<4;enlarge_factor = 15;end

ND = length(size(I));
if ND > 3
    error('The image dimension is too large!')
end

% normalize the matrix to [0,1]
min_value = min(I(:));
max_value = max(I(:));
if (min_value ~= 0) || (max_value ~= 1)
    I = (I - min_value) / (max_value - min_value);
end

if ND == 3 % 3D image
    % Choose a suitable initial contour or mask (you may need to customize this)
    maxint =max(max(max(abs(I)))); level = thres * maxint;
    initial_mask = I>level;
    for i = 1:size(initial_mask,3)
        % delete the isolated small region
        initial_mask(:, :, i) = bwareaopen(initial_mask(:, :, i), isolated_region_size);
        % fill the holes
        initial_mask(:,:,i) = imfill(initial_mask(:,:,i),'holes');
    end
    % enlarge the mask a little bit
    se = strel('disk', enlarge_factor); % Adjust the size of the disk based on how much you want to enlarge the mask
    initial_mask = imdilate(initial_mask, se);
%     figure,imshow3D(initial_mask.*I,[])
else % 2D image
    maxint =max(max(abs(I))); level = thres * maxint;
    initial_mask = I>level;
    % delete the isolated small region
    initial_mask = bwareaopen(initial_mask, isolated_region_size);
    % fill the holes
    initial_mask = imfill(initial_mask,'holes');
    % enlarge the mask a little bit
    se = strel('disk', enlarge_factor); % Adjust the size of the disk based on how much you want to enlarge the mask
    initial_mask = imdilate(initial_mask, se);
%     figure,imshow(initial_mask.*I,[])
end

mask = initial_mask;