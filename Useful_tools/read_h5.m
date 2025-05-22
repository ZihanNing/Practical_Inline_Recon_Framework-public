clc
clear all
% close all
home_path = getenv('HOME');
addpath(genpath([home_path,'/matlab/usual_used']))
data_name = 'output_MPRAGE_dummy.h5';
hinfo = hdf5info(data_name);
header = h5read(data_name,hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(3).Name);
img = h5read(data_name,hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name);

fprintf('============ Type the data type =============\n');
disp(['Data type: ',class(img)])

img = double(img);
figure,imshow3D(img,[0 1500])

%% read multi-echo images
clear all
% close all
home_path = getenv('HOME');
addpath(genpath([home_path,'/matlab/usual_used']))
data_name = 'output_SWI_phantom.h5';
hinfo = hdf5info(data_name);
series_num = 1; % MAG or PHS
image_index = 1; % echo 1 & 2
series_names = hinfo.GroupHierarchy.Groups.Groups(1);
% series_names = hinfo.GroupHierarchy.Groups.Groups(series_num);
selected_series_name = series_names.Name;

% Read the images from the selected series
images = h5read(data_name, [selected_series_name '/data']);

fprintf('============ Type the data type =============\n');
disp(['Data type: ',class(images)])
disp(['Data size:',size(images)])

img = double(squeeze(images));
figure,imshow3D(img(:,:,:,1),[])

%% compare the header
fprintf('============ Type the header differences =============\n');
load('header_generic.mat')
struct1 = header;
load('header_debug_new.mat')
struct2 = header;
% compare
com_info = structcmp(struct1,struct2,'Report','on');
% show
disp(com_info);


