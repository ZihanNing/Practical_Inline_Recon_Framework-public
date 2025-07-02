function varargout = bview(varargin)
% BVIEW M-file for bview.fig
%
%   Simple image viewing and ROI drawing GUI.
%
%   (c) Corey Baron, 2010
%
%      BVIEW, by itself, creates a new BVIEW or raises the existing
%      singleton*.
%
%      H = BVIEW returns the handle to a new BVIEW or the handle to
%      the existing singleton*.
%
%      BVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVIEW.M with the given input arguments.
%
%      BVIEW('Property','Value',...) creates a new BVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bview

% Last Modified by GUIDE v2.5 04-Oct-2019 10:57:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bview_OpeningFcn, ...
                   'gui_OutputFcn',  @bview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before bview is made visible.
function bview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bview (see VARARGIN)

% Choose default command line output for bview
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Define the Image Type 
handles.imgopt = 1;
handles.img2 = [];
handles.img3 = []; % 
handles.img4 = []; %
handles.varlength = 1; 
% Populate the listbox
update_listbox(handles)
set(handles.listbox1,'Value',1)
set(handles.listbox2,'Value',1) 

% Define the beginning extraction type (i.e. abs, phase, real, imag)
handles.extractiontype = 1;

% Define the beginning save type (i.e. jpg, tif)
handles.savetype = 1;

% Define the beginning scale type
handles.scaletype = 1;

% Define interpolation factor
handles.interpfactor = 1;

% Define aspect ratio
handles.aspectratio = 1;
handles.aspecttype = 1;

% Set view to be axial
handles.view = 1;

%Set Colourmap to be grayscale 
handles.colour = 1; 

% Define ROI field
handles.ROI = [];
handles.ROI2 = [];
handles.ROI_3d = [];
handles.ROI_3d_par = [];
handles.ROI2_3d = [];
handles.ROI2_3d_par = [];
for i = 3:6 
    handlesname = sprintf('handles.ROI%d_3d',i); 
    eval(sprintf('%s = [];',handlesname)) 
    eval(sprintf('%s_par = [];',handlesname))
end 
%single roi as default 
handles.mroi = 1; 
% Define a colormap
zero_range = 0.1;
nc = 256;
handles.cmap = parula(nc); %winter %hot %jet % summer %hsv %cooper
% nz = round(256*zero_range)/2;
% handles.cmap(nc/2-nz:nc/2+nz,:) = repmat([0 0 0],[2*nz+1 1]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bview wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = bview_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[var1,var2,~,~] = get_var_name(handles); % [var1,~, ~, ~]
handles.img = evalin('base',var1);  
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
imsize = size(handles.img_ex);
handles.mask = ones(imsize(1),imsize(2));
if handles.extractiontype == 5
    handles.mask = ones(imsize(2),imsize(3));
end
% update min and max text boxes.
adjustsliders(handles)
if strcmp(get(handles.edit_min,'String'),'-') || (handles.scaletype == 2)
    array_select = round(get(handles.slider2,'Value'));
    array_select2 = round(get(handles.slider5,'Value'));
    handles.maxvalue = max(max(max(handles.img_ex(:,:,:,array_select,array_select2))));
    handles.minvalue = min(min(min(handles.img_ex(:,:,:,array_select,array_select2))));
    set(handles.edit_min,'String', handles.minvalue)
    set(handles.edit_max,'String', handles.maxvalue)
end

% Initialize 3D ROI if necessary
sz_tmp = size(handles.img_ex);
if length(sz_tmp) < 3
    sz_tmp(3) = 1;
end
if isempty(handles.ROI_3d) %|| (size(handles.ROI_3d) ~= sz_tmp)
    handles.ROI_3d = false(sz_tmp(1:3));
    handles.ROI_3d_par = cell(sz_tmp(3),1);
end
% Enable Multiple/Overlap Options once handles.img is obtained
set(handles.radiobutton_Multiple,'Enable','on') 
set(handles.radiobutton_overlap,'Enable','on') 
% Update handles structure
guidata(hObject, handles);

% Show the image
    showimage(hObject,handles)
% Calculate SNR if ROI's have already been selected
SNR = calcSNR3D(handles);
set(handles.text_SNR,'String',sprintf('%g',SNR))
% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
function adjustsliders(handles)
if (handles.extractiontype == 5) && (size(handles.img_ex,1)==3)
    dim1 = 4;
    dim2 = 5;
    dim3 = 6;
else
    dim1 = 3;
    dim2 = 4;
    dim3 = 5;
end

% Make adjustments for size of image structure
if get(handles.slider1,'Value') > size(handles.img_ex,dim1)
    set(handles.slider1,'Value',size(handles.img_ex,dim1))
end
if get(handles.slider2,'Value') > size(handles.img_ex,dim2)
    set(handles.slider2,'Value',size(handles.img_ex,dim2))
end
if ~iscell(handles.img) && (get(handles.slider5,'Value') > size(handles.img_ex,dim3))
    set(handles.slider5,'Value',size(handles.img_ex,dim3))
end
if  size(handles.img_ex,dim1) > 1
    set(handles.slider1,'Enable','on')
    set(handles.slider1,'Max',size(handles.img_ex,dim1))
    set(handles.slider1,'SliderStep',[1/(size(handles.img_ex,dim1)-1),1/(size(handles.img_ex,dim1)-1)])
else
    set(handles.slider1,'Enable','off')
    set(handles.slider1,'Value',1)
end
if  size(handles.img_ex,dim2) > 1
    set(handles.slider2,'Enable','on')
    set(handles.slider2,'Max',size(handles.img_ex,dim2))
    set(handles.slider2,'SliderStep',[1/(size(handles.img_ex,dim2)-1),1/(size(handles.img_ex,dim2)-1)])
else
    set(handles.slider2,'Enable','off')
    set(handles.slider2,'Value',1)
end
if  size(handles.img_ex,dim3) > 1
    set(handles.slider5,'Enable','on')
    set(handles.slider5,'Max',size(handles.img_ex,dim3))
    set(handles.slider5,'SliderStep',[1/(size(handles.img_ex,dim3)-1),1/(size(handles.img_ex,dim3)-1)])
elseif iscell(handles.img)
    set(handles.slider5,'Enable','on')
    set(handles.slider5,'Max',prod(size(handles.img)))
    set(handles.slider5,'SliderStep',[1/(prod(size(handles.img))-1),1/(prod(size(handles.img))-1)])
else
    set(handles.slider5,'Enable','off')
    set(handles.slider5,'Value',1)
end
set(handles.text1,'String',sprintf('Slice: %d',round(get(handles.slider1,'Value'))))
set(handles.text3,'String',sprintf('Array Dimension: %d',round(get(handles.slider2,'Value'))))
set(handles.text16,'String',sprintf('2nd Array Dimension: %d',round(get(handles.slider5,'Value'))))
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles) % list for the second dialog box
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
[var1,var2,var3,var4] = get_var_name(handles); %([var1 var2 var3 var4] = get_var_name(handles) 
%handles.img = double(evalin('base',var1)); 
if handles.imgopt == 2  
     arrayfun(@cla, findall(0,'type','axes')) 
     for i = 2:4 
        eval(sprintf('handles.img%d = [];',i)) 
     end 
     vars = {var1,var2,var3,var4};  
     vars(cellfun('isempty',vars)) = []; 
     handles.varlength = length(vars);      
for i = 2:length(vars) 
    eval(sprintf('handles.img%d = evalin(''base'',''%s'');',i,vars{i}))
    eval(sprintf('handles.img_ex%d = extraction(handles.slider5,handles.img%d,handles.extractiontype,handles.view);',i,i))
    eval(sprintf('imsize%d = size(handles.img_ex%d);',i,i))
    eval(sprintf('handles.mask%d = ones(imsize%d(1), imsize%d(2));',i,i,i))
    if handles.extractiontype == 5 
        eval(sprintf('handles.mask%d = ones(imsize%d(2),imsize%d(3));',i,i,i))
    end 
end
clear vars var1 var2 var3 var4 
end 
if handles.imgopt == 3 %Overlap 
 handles.img2 = evalin('base',var2);
 handles.img_ex2 = extraction(handles.slider5,handles.img2,handles.extractiontype,handles.view);
 imsize = size(handles.img_ex2);
 handles.mask2 = ones(imsize(1),imsize(2));
 if handles.extractiontype == 5
     handles.mask2 = ones(imsize(2),imsize(3));
 end
end 
% update min and max text boxes.       
adjustsliders(handles)
if strcmp(get(handles.edit9,'String'),'-') || (handles.scaletype == 2) 
    array_select = round(get(handles.slider2,'Value'));
    array_select2 = round(get(handles.slider5,'Value'));
    handles.maxvalue2 = max(max(max(handles.img_ex2(:,:,:,array_select,array_select2)))); % adjust for a new set of min and max depending on the image
    handles.minvalue2 = min(min(min(handles.img_ex2(:,:,:,array_select,array_select2)))); 
    set(handles.edit9,'String', handles.minvalue2)
    set(handles.edit_max2,'String', handles.maxvalue2)
end 
%Check if Array Sizes Match for Multiple Images(Multiple/Overlap options) 
if ~isempty(handles.img2) 
    if (handles.extractiontype == 5) && (size(handles.img_ex,1)==3)
        dim1 = 4;
        dim2 = 5;
        dim3 = 6;
    else
        dim1 = 3;
        dim2 = 4;
        dim3 = 5;
    end
    for i = 2:4
        if eval(sprintf('~isempty(handles.img%d)',i)) 
        if eval(sprintf('size(handles.img_ex,dim1) ~= size(handles.img_ex%d,dim1) | size(handles.img_ex,dim2) ~= size(handles.img_ex%d,dim2) | size(handles.img_ex,dim3) ~= size(handles.img_ex%d,dim3)',i,i,i)) 
            errordlg('Array Dimensions Do Not Match','Array Dimensions Do Not Match','modal') 
        end 
        end 
    end 
end 
%set min and max values if var3 and var4 upon selections
for i = 3:4
    if eval(sprintf('~isempty(handles.img%d)',i)) 
          array_select = round(get(handles.slider2,'Value'));
          array_select2 = round(get(handles.slider5,'Value'));   
          eval(sprintf('handles.maxvalue%d = max(max(max(handles.img_ex%d(:,:,:,array_select,array_select2))));',i,i))
          eval(sprintf('handles.minvalue%d = min(min(min(handles.img_ex%d(:,:,:,array_select,array_select2))));',i,i))
        end 
 end 
%Enable Multi_Image tools depending on number of images formed 
if ~isempty(handles.img3) && (handles.varlength == 3 || handles.varlength == 4) 
    set(handles.multi_image_tools,'Visible','on')
    set(handles.reset_image_3, 'Visible', 'on')
    set(handles.text_minmax3, 'String',sprintf('Min3: %d Max3: %d',handles.minvalue3,handles.maxvalue3))
    set(handles.text_minmax3,'Visible','on')
else    
    set(handles.multi_image_tools,'Visible','off')
    set(handles.reset_image_3, 'Visible', 'off')
    set(handles.text_minmax3,'Visible','off')
end 
if ~isempty(handles.img4) && handles.varlength == 4 
    set(handles.reset_image_4, 'Visible', 'on')
    set(handles.text_minmax4, 'String',sprintf('Min4: %d Max4: %d',handles.minvalue4,handles.maxvalue4))
    set(handles.text_minmax4,'Visible','on')
else 
    set(handles.reset_image_4, 'Visible', 'off')
    set(handles.text_minmax4,'Visible','off')
end 
% Initialize 3D ROI if necessary
sz_tmp = size(handles.img_ex);
if length(sz_tmp) < 3
    sz_tmp(3) = 1;
end
if isempty(handles.ROI_3d) %|| (size(handles.ROI_3d) ~= sz_tmp)
    handles.ROI_3d = false(sz_tmp(1:3));
    handles.ROI_3d_par = cell(sz_tmp(3),1);
end
set(handles.text24,'String','Info')
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)
% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in pushbutton_updatevar.

function pushbutton_updatevar_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_updatevar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1,'Value',1);
update_listbox(handles)
if handles.imgopt == 2 || handles.imgopt == 3 
    set(handles.listbox2,'Value',1)
    update_listbox(handles)
end 
function update_listbox(handles)
% handles    structure with handles and user data (see GUIDATA)
% Updates the listbox to match the current workspace
vars = evalin('base','who');
set(handles.listbox1,'String',vars)
if handles.imgopt == 2 || handles.imgopt == 3 
    vars = evalin('base','who'); 
    set(handles.listbox2,'String',vars)
end 

function [var1,var2,var3,var4] = get_var_name(handles) % [var1,var2,~,~] = get_var_names(handles) 
% Returns the names of the image variable to plot
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');

if length(index_selected) > 1 
    errordlg('You must select one variable','Incorrect Selection','modal')
elseif length(index_selected) == 2 
    var1 = list_entries{index_selected(1)}; 
    var2 = list_entries{index_selected(2)}; 
elseif length(index_selected) == 1
    var1 = list_entries{index_selected(1)};
    var2 = [];
    var3 = []; 
    var4 = []; 
end
if handles.imgopt == 2 || handles.imgopt == 3 
    list_entries2 = get(handles.listbox2,'String'); 
    index_selected2 = get(handles.listbox2,'Value'); 
if (length(index_selected2) > 1 && handles.imgopt == 3) || length(index_selected2) > 3
    errordlg('You must select one variable','Incorrect Selection','modal')
elseif length(index_selected2) == 3 
    var2 = list_entries2{index_selected2(2)};   
    var3 = list_entries2{index_selected2(3)};   
    var4 = list_entries2{index_selected2(1)};   
elseif length(index_selected2) == 2 
    var2 = list_entries2{index_selected2(1)}; 
    var3 = list_entries2{index_selected2(2)}; 
    var4 = []; 
elseif length(index_selected2) == 1 
    var2 = list_entries2{index_selected2(1)};
    var3 = []; 
    var4 = [];  
end
end 
clear list_entries list_entries2 index_selected2 index_selected
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showimage(hObject,handles)
set(handles.text1,'String',sprintf('Slice: %d',round(get(handles.slider1,'Value'))))
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles) % Slice Selection 
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function output = extraction(slider5,input,type,view)
if iscell(input)
  array_select2 = round(get(slider5,'Value'));
  input = input{array_select2};
end
% Apply contrast change
if type == 1
    output = abs(input);
elseif type == 2
    output = angle(input);
elseif type == 3
    if (size(input,1) ~= 6) && (size(input,1) ~= 3)
        errordlg('You must use DTI data','Incorrect Selection','modal')
        output = input;
    else
        output = squeeze(input(1,:,:,:,:,:) + input(2,:,:,:,:,:) + input(3,:,:,:,:,:));
    end
elseif type == 4
    if size(input,1) ~= 3
        errordlg('You must use DTI eigenvalues','Incorrect Selection','modal')
        output = input;
    else
        trace = ( input(1,:,:,:,:,:) + input(2,:,:,:,:,:) + input(3,:,:,:,:,:) )/3;
        output = squeeze( sqrt(3/2)*...
            sqrt( (input(1,:,:,:,:,:)-trace).^2 +...
                (input(2,:,:,:,:,:)-trace).^2 +...
                (input(3,:,:,:,:,:)-trace).^2 )./...
            sqrt( input(1,:,:,:,:,:).^2 + input(2,:,:,:,:,:).^2 + input(3,:,:,:,:,:).^2 ));
        clear trace;
    end
elseif type == 5
    output = abs(input);
    imsize = size(output);
    if imsize(end) == 3
        parray = circshift(1:length(imsize),[0 1]);
        output = permute(output,parray);
        temp = output(1,:,:,:,:);
        output(1,:,:,:,:) = output(2,:,:,:,:);
        output(2,:,:,:,:) = temp;
        clear temp
    end
elseif type == 6
    output = log(1+0.00001*abs(input));
elseif type == 7
    output = real(input);
elseif type == 8
    output = imag(input);
else
    output = input;
end
% Apply view
if type~=5
    if view == 2
        output = permute(input,[3 1 2 4 5 6]);
        output = flipdim(output,1);
    elseif view == 3
        output = permute(input,[3 2 1 4 5 6]);
        output = flipdim(output,1);
    end
end

function output = FA_calc(i1,i2,i3,i4,i5,i6)
D = [i1,i4,i5;...
    i4,i2,i6;...
    i5,i6,i3];
D = eig(D);
meanD = sum(D)/3;
output = sqrt(3/2)*sqrt((D(1)-meanD)^2 + (D(2)-meanD)^2 + (D(3)-meanD)^2)/...
    sqrt(D(1)^2+D(2)^2+D(3)^2);
clear meanD D

function showimage(hObject,handles)
% Retrieve relevant parameters for slice and array selection
slice_select = round(get(handles.slider1,'Value'));
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
% wnd = get(handles.slider3,'Value');
% lvl = get(handles.slider4,'Value');

if iscell(handles.img)
  array_select2 = 1;
  handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
  for i = 2:4
     if eval(sprintf('~isempty(handles.img%d)',i))   
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5,handles.img%d,handles.extractiontype,handles.view);',i,i))
     end 
  end 
end

if (handles.extractiontype == 5) && (size(handles.img_ex,1)==3)
%     imsize = size(handles.img_ex);
%     parray = circshift(1:length(imsize),[0 1]);
%     img_tmp = handles.img_ex(:,:,slice_select,array_select,array_select2);
    imsize = size(squeeze(handles.img_ex(1,:,:,slice_select,array_select,array_select2)));
    if min( imsize == size(handles.mask) )
        image( permute( handles.img_ex(:,:,:,slice_select,array_select,array_select2), [2 3 1] ).*...
            repmat(handles.mask, [1 1 3] ))
    else
        image( permute( handles.img_ex(:,:,:,slice_select,array_select,array_select2), [2 3 1] ) )
    end
    axis image
else
% Scale the image
    maxlevel = handles.maxvalue;
    minlevel = handles.minvalue;
    for i = 2:4
        if eval(sprintf('~isempty(handles.img%d)',i)) %handles.imgopt == 2 || handles.imgopt == 3  
            eval(sprintf('maxlevel%d = handles.maxvalue%d;',i,i)) % maxlevel2 = handles.maxvalue2; 
            eval(sprintf('minlevel%d = handles.minvalue%d;',i,i))  %minlevel2 = handles.minvalue2;
        end 
    end 
    % Apply the mask, if present
    if min(size(handles.img_ex(:,:,slice_select,array_select,array_select2)) == size(handles.mask) )
        img = handles.img_ex(:,:,slice_select,array_select,array_select2).*handles.mask;
        for i = 2:4 
            if eval(sprintf('~isempty(handles.img%d)',i)) && eval(sprintf('min(size(handles.img_ex%d(:,:,slice_select,array_select,array_select2)) == size(handles.mask%d))',i,i))  
                eval(sprintf('img%d = handles.img_ex%d(:,:,slice_select,array_select,array_select2).*handles.mask%d;',i,i,i))  
            end 
        end 
    else
        img = handles.img_ex(:,:,slice_select,array_select,array_select2);
        for i = 2:4 
            if eval(sprintf('~isempty(handles.img%d)',i))  
                eval(sprintf('img%d = handles.img_ex%d(:,:,slice_select,array_select,array_select2);',i,i)) 
            end 
        end 
    end

    % Apply the aspectrato
    if handles.aspecttype == 1
        % Here handles.aspectratio gives the in-plane aspect ratio of the voxels
        mfact = handles.aspectratio;
        for i = 2:4 
        if eval(sprintf('~isempty(handles.img%d)',i))   
            eval(sprintf('mfact%d = handles.aspectratio;',i))  
        end 
        end 
    else
        % Here handles.aspectratio gives the in-plane aspect ratio of the FOV
        mfact = handles.aspectratio*size(handles.img_ex,2)/size(handles.img_ex,1);
        for i = 2:4 
        if eval(sprintf('~isempty(handles.img%d)',i))    
            eval(sprintf('mfact%d = handles.aspectratio*size(handles.img_ex%d,2)/size(handles.img_ex%d,1);',i,i,i)) 
        end 
        end 
    end
    if mfact < 1
        newmatrix = [size(handles.img_ex,1) round(size(handles.img_ex,2)/mfact)];
    else
        newmatrix = [round(size(handles.img_ex,1)*mfact) size(handles.img_ex,2)];
    end
    for i = 2:4
        if  eval(sprintf('~isempty(handles.img%d)',i))
            if eval(sprintf('mfact%d < 1',i))  
                eval(sprintf('newmatrix%d = [size(handles.img_ex%d,1) round(size(handles.img_ex%d,2)/mfact)];',i,i,i,i))
            else
                eval(sprintf('newmatrix%d = [round(size(handles.img_ex%d,1)*mfact%d) size(handles.img_ex%d,2)];',i,i,i,i))
            end
                eval(sprintf('img%d = imresize(img%d,newmatrix%d,''nearest'');',i,i,i))
        end
    end 
    img = imresize(img,newmatrix,'nearest');
    
% Scale data
    maxthresh = img <= maxlevel;
    minthresh = img >= minlevel;
    img = ((img.*maxthresh.*minthresh + (1-maxthresh)*maxlevel +...
        (1-minthresh)*minlevel) - minlevel) / (maxlevel-minlevel);
    for i = 2:4 
        if eval(sprintf('~isempty(handles.img%d)',i)) 
            eval(sprintf( 'maxthresh%d = img%d <= maxlevel%d;',i,i,i)) 
            eval(sprintf('minthresh%d = img%d >= minlevel%d;',i,i,i)) 
            eval(sprintf('img%d = ((img%d.*maxthresh%d.*minthresh%d + (1-maxthresh%d)*maxlevel%d + (1-minthresh%d)*minlevel%d) - minlevel%d) / (maxlevel%d-minlevel%d);',i,i,i,i,i,i,i,i,i,i,i))
        end 
    end 
% show the image
    if ~isreal(img)
        img = abs(img);
        for i = 2:4 
            if eval(sprintf('~isempty(handles.img%d)',i)) 
               eval(sprintf('img%d = abs(img%d);',i,i))  
            end 
        end 
    end
%Plotting the Images     
    if handles.colour == 2  %Colour selected or Real data type
%         imshow(img,handles.cmap,'InitialMagnification','fit')
        if handles.imgopt == 3 % Overlap
            set(gca,'Position',[0.1562,0.079,0.625,0.8180]) %possibly use reset instead
            I = imshow(img,'Colormap',handles.cmap,'InitialMagnification','fit'); 
            red = cat(3,ones(size(img)),zeros(size(img)),zeros(size(img)));
            hold on 
            h = imshow(red,'InitialMagnification','fit'); 
            set(h,'AlphaData',img2); %Transparency Map 
            hold off 
        elseif handles.imgopt == 2 %Multiple 
            subplot(round((handles.varlength/4)*2),round(handles.varlength/2)*2,1)
                imshow(img,'Colormap',handles.cmap,'InitialMagnification','fit');
                set(gca,'Position',[0.1562+((1-rem(1,2))*0.3125), 0.079+((round(1/5))*0.490),0.3125,0.8180/round(handles.varlength/2)])
                showroi(handles,slice_select,img)
           for i = 2:4
               if eval(sprintf('~isempty(handles.img%d)',i)) 
               subplot(round((handles.varlength/4)*2),round(handles.varlength/2)*2,i)
                  eval(sprintf('imshow(img%d,''Colormap'',handles.cmap,''InitialMagnification'',''fit'')',i))
                  set(gca,'Position',[0.1562+((1-rem(i,2))*0.3125), 0.079+((round(i/5))*0.490),0.3125,0.8180/round(handles.varlength/2)])
                  showroi(handles,slice_select,img) 
               end 
           end       
        else %Single 
           set(gca,'Position',[0.1562,0.079,0.625,0.8180])
           imshow(img,'Colormap',handles.cmap,'InitialMagnification','fit')
        end 
    else %Grayscale
        if handles.imgopt == 3 %Overlap
            set(gca,'Position',[0.1562,0.079,0.625,0.8180])
            I = imshow(img,[0,1],'InitialMagnification','fit');   
            red = cat(3,ones(size(img)),zeros(size(img)), zeros(size(img)));
            hold on 
            h = imshow(red,'InitialMagnification','fit'); 
            set(h,'AlphaData',img2);
            hold off 
        elseif handles.imgopt == 2  %Multiple
            subplot(round((handles.varlength/4)*2),round(handles.varlength/2)*2,1)
                imshow(img,[0 1],'InitialMagnification','fit');
                set(gca,'Position',[0.1562+((1-rem(1,2))*0.3125), 0.079+((round(1/5))*0.490),0.3125,0.8180/round(handles.varlength/2)])
                showroi(handles,slice_select,img)
           for i = 2:4
               if eval(sprintf('~isempty(handles.img%d)',i)) 
                 subplot(round((handles.varlength/4)*2),round(handles.varlength/2)*2,i)
                  eval(sprintf('imshow(img%d,[0 1],''InitialMagnification'',''fit'')',i))
                  set(gca,'Position',[0.1562+((1-rem(i,2))*0.3125), 0.079+((round(i/5))*0.490),0.3125,0.8180/round(handles.varlength/2)])
                  showroi(handles,slice_select,img) 
               end 
           end 
        else  %Single (Default Image) 
            set(gca,'Position',[0.1562,0.079,0.625,0.8180])
            imshow(img,[0 1],'InitialMagnification','fit')
%         imshow(img,[0 1],'InitialMagnification','fit')
        end 
    end
%Show volume ROIs if present
        showroi(handles,slice_select,img)
        clear img clear img2 img3 img4 
end

function showroi(handles,slice_select,img) %shows all ROI Selected
    sz = size(handles.ROI_3d); % ROI (default)  
    if (slice_select <= length(handles.ROI_3d_par)) && ~isempty(handles.ROI_3d_par{slice_select}) && isequal(size(img),sz(1:2)) 
      hold on   
           for n=1:length(handles.ROI_3d_par{slice_select})
             plot(handles.ROI_3d_par{slice_select}{n}(:,1),...
             handles.ROI_3d_par{slice_select}{n}(:,2), 'y')
           end
       hold off 
    end
    for k = 2:6 %Additional ROIs 
        roiname = sprintf('handles.ROI%d_3d',k);
        colours = {'y','m','g','r','b','c'}; 
        eval(sprintf('sz = size(%s);',roiname))
        if eval(sprintf('(slice_select <= length(%s_par))',roiname)) && eval(sprintf('~isempty(%s_par{slice_select}) && isequal(size(img),sz(1:2))',roiname)) 
         hold on   
           for n=1:eval(sprintf('length(%s_par{slice_select})',roiname))
             eval(sprintf('plot(%s_par{slice_select}{n}(:,1),%s_par{slice_select}{n}(:,2),''%s'')',roiname,roiname,colours{k}))
           end
         hold off 
        end
     end 
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles) %Array Dimension 1 
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showimage(hObject,handles)
set(handles.text3,'String',sprintf('Array Dimension: %d',round(get(handles.slider2,'Value'))))
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showimage(hObject,handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showimage(hObject,handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton3. 
function pushbutton3_Callback(hObject, eventdata, handles) % Reset Min and Max scales
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[minval,maxval,~,~] = resetscale(handles); %Callback function on line ... 1118
handles.maxvalue = maxval;
handles.minvalue = minval; 
set(handles.edit_min,'String', handles.minvalue)
set(handles.edit_max,'String', handles.maxvalue)
%Show the Image 
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_ave_allsl. (Removed Before 2014 Edition) 
function pushbutton_ave_allsl_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ave_allsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.ROI)
    errordlg('Select an ROI first','Incorrect Selection','modal')
% elseif size(handles.img_ex,4) == 1
%     slice_select = round(get(handles.slider1,'Value'));
%     array_select = round(get(handles.slider2,'Value'));
%     array_select2 = round(get(handles.slider5,'Value'));
%     numelements  = length(find(handles.ROI));
%     ave = squeeze(sum(sum(handles.img_ex(:,:,slice_select,array_select,array_select2).*...
%         handles.ROI,1),2)) / numelements
else
%    slice_select = round(get(handles.slider1,'Value'));
%     mxarray      = round(get(handles.slider2,'Max'));
%     array_select2 = round(get(handles.slider5,'Value'));
%     numelements  = length(find(handles.ROI));
%     ave = squeeze(sum(sum(handles.img_ex(:,:,slice_select,:,array_select2).*...
%         repmat(handles.ROI,[1 1 1 mxarray 1]),1),2)) / numelements
    %array_select = round(get(handles.slider2,'Value'));
    mxarray      = round(get(handles.slider1,'Max'));
    mxarray2      = round(get(handles.slider2,'Max'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    ave = squeeze(sum(sum(handles.img_ex(:,:,:,:,array_select2).*...
        repmat(handles.ROI,[1 1 mxarray mxarray2 1]),1),2)) / numelements
    save ave ave
    clear ave
end

% --- Executes on button press in pushbutton_DT_allsl. (Removed Before 2014 Edition)
function pushbutton_DT_allsl_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DT_allsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.ROI)
    errordlg('Select an ROI first','Incorrect Selection','modal')
elseif (size(handles.img,1) ~= 6) && (size(handles.img,1) ~= 3)
    errordlg('You must use DTI data','Incorrect Selection','modal')
elseif size(handles.img,5) == 1
    slice_select = round(get(handles.slider1,'Value'));
    array_select = round(get(handles.slider2,'Value'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    ave = squeeze(sum(sum(handles.img(:,:,:,slice_select,array_select,array_select2).*...
        repmat(reshape(handles.ROI,[1 size(handles.ROI,1) size(handles.ROI,2)]),...
        [size(handles.img,1) 1 1]),2),3)) / numelements
else
    slice_select = round(get(handles.slider1,'Value'));
    mxarray      = round(get(handles.slider2,'Max'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    ave = squeeze(sum(sum(handles.img(:,:,:,slice_select,:,array_select2).*...
        repmat(reshape(handles.ROI,[1 size(handles.ROI,1) size(handles.ROI,2)]),...
        [size(handles.img,1) 1 1 1 mxarray 1]),2),3)) / numelements
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles) % The Aspect Ratio
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.interpfactor = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement. 
function slider5_Callback(hObject, eventdata, handles) % Second Array Dimension 
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showimage(hObject,handles)
set(handles.text16,'String',sprintf('2nd Array Dimension: %d',round(get(handles.slider5,'Value'))))
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton_std_allsl. (%Removed)
function pushbutton_std_allsl_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_std_allsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.ROI)
    errordlg('Select an ROI first','Incorrect Selection','modal')
% elseif size(handles.img_ex,4) == 1
%     slice_select = round(get(handles.slider1,'Value'));
%     array_select = round(get(handles.slider2,'Value'));
%     array_select2 = round(get(handles.slider5,'Value'));
%     numelements  = length(find(handles.ROI));
%     ave = squeeze(sum(sum(handles.img_ex(:,:,slice_select,array_select,array_select2).*...
%         handles.ROI,1),2)) / numelements
else
%     slice_select = round(get(handles.slider1,'Value'));
%     mxarray      = round(get(handles.slider2,'Max'));
%     array_select2 = round(get(handles.slider5,'Value'));
%     numelements  = length(find(handles.ROI));
%     ave = squeeze(sum(sum(handles.img_ex(:,:,slice_select,:,array_select2).*...
%         repmat(handles.ROI,[1 1 1 mxarray 1]),1),2)) / numelements

    %array_select = round(get(handles.slider2,'Value'));
    mxarray      = round(get(handles.slider1,'Max'));
    mxarray2      = round(get(handles.slider2,'Max'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    SD = zeros(size(handles.img_ex,3),size(handles.img_ex,4));
    for i=1:mxarray
        for j=1:mxarray2
            temp = handles.img_ex(:,:,i,j,array_select2);
            SD(i,j) = std(temp(handles.ROI));
%             SD = squeeze(sum(sum(handles.img_ex(:,:,:,:,array_select2).*...
%                 repmat(handles.ROI,[1 1 mxarray mxarray2 1]),1),2)) / numelements
        end
    end
    SD
    save SD SD
    clear temp SD
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aspectratio = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.minvalue = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Show image
showimage(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxvalue = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Show image
showimage(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel3. (Auto-Scale) 
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_nvr'
        % Code for when radiobutton_nvr is selected.
        handles.scaletype = 1;
    case 'radiobutton_newselect'
        % Code for when radiobutton_newselect is selected.
        handles.scaletype = 2;
end
% Update handles structure
guidata(hObject, handles);

function [minval,maxval,minval2,maxval2] = resetscale(handles)
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
if iscell(handles.img)
  array_select2 = 1;
end
maxval = max(max(max(handles.img_ex(:,:,:,array_select,array_select2))));
minval = min(min(min(handles.img_ex(:,:,:,array_select,array_select2))));
maxval2 = '-'; %maintains state when handles is empty
minval2 ='-';  
if ~isempty(handles.img2) 
maxval2 =  max(max(max(handles.img_ex2(:,:,:,array_select,array_select2)))); 
minval2 = min(min(min(handles.img_ex2(:,:,:,array_select,array_select2)))); 
end 
    
% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(hObject,'String','FOV')
    handles.aspecttype = 0;
else
    set(hObject,'String','voxel')
    handles.aspecttype = 1;
end
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton1

% --- Executes on button press in pushbutton21. (Add ROI)
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.mroi == 2 
roi_select = round(get(handles.roi_selection_slider,'Value')); 
else 
roi_select = 1; 
end 
if handles.imgopt == 2 
    set(handles.text24,'String','Click on an Image to Select Axes for ROI')
    w = waitforbuttonpress; 
    if w == 0 
        set(handles.text24,'String','Select ROI') 
    else
        errordlg('Axes Not Selected!','Incorrect Selection','modal')
    end 
end 
slice_select = round(get(handles.slider1,'Value'));
[roi_sub,xi,yi] = roipoly;
sz = size(handles.img_ex);
if length(sz)<3
    sz = [sz,1];
end
szR = size(handles.ROI_3d);
if length(szR)<3
    szR = [szR,1];
end
if roi_select == 1  
if ~isequal(sz(1:3),szR)
    handles.ROI_3d = false(sz(1:3));
    handles.ROI_3d_par = cell(sz(3),1);
end 
    next = length(handles.ROI_3d_par{slice_select}) + 1;
    handles.ROI_3d_par{slice_select}{next} = [xi,yi];
    handles.ROI_3d(:,:,slice_select) = or(roi_sub,handles.ROI_3d(:,:,slice_select));
    set(handles.text24,'String','Info')
else
    roiname = sprintf('handles.ROI%d_3d',roi_select); 
    if eval(sprintf('~isequal(sz(1:3),size(%s))',roiname))
    eval(sprintf('%s = false(sz(1:3));',roiname))
    eval(sprintf('%s_par = cell(sz(3),1);',roiname))
    end
    eval(sprintf('next = length(%s_par{slice_select}) + 1;',roiname))
    eval(sprintf('%s_par{slice_select}{next} = [xi,yi];',roiname))
    eval(sprintf('%s(:,:,slice_select) = or(roi_sub,%s(:,:,slice_select));',roiname,roiname))
    set(handles.text24,'String','ROI1:Yellow, ROI2:Magenta, ROI3:Green, ROI4:Red, ROI5:Blue, ROI6:Cyan')
end 
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)
clear roi_sub xi yi sz slice_select

% --- Executes on button press in pushbutton28. (Add ROI2) 
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.imgopt == 2 %if Multiple Image option 
    set(handles.text24,'String','Click on an Image to Select Axes for ROI')
    w = waitforbuttonpress; 
    if w == 0 
        set(handles.text24,'String','Select ROI') 
    else
        errordlg('Axes Not Selected!','Incorrect Selection','modal')
    end 
end
slice_select = round(get(handles.slider1,'Value'));
[roi_sub,xi,yi] = roipoly;
sz = size(handles.img_ex);
if length(sz)<3
    sz = [sz,1];
end
szR = size(handles.ROI2_3d);
if length(szR)<3
    szR = [szR,1];
end
if ~isequal(sz(1:3),szR)
    handles.ROI2_3d = false(sz(1:3));
    handles.ROI2_3d_par = cell(sz(3),1);
end
next = length(handles.ROI2_3d_par{slice_select}) + 1;
handles.ROI2_3d_par{slice_select}{next} = [xi,yi];
handles.ROI2_3d(:,:,slice_select) = or(roi_sub,handles.ROI2_3d(:,:,slice_select));
set(handles.text24,'String','Info')
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)
clear roi_sub xi yi sz slice_select

% --- Executes on button press in pushbutton22. (Clear Slice)
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slice_select = round(get(handles.slider1,'Value'));
handles.ROI_3d(:,:,slice_select) = false(size(handles.ROI_3d,1),size(handles.ROI_3d,2));
handles.ROI_3d_par{slice_select} = [];
for k = 2:6 %Clears all additional ROI 
    roiname = sprintf('handles.ROI%d_3d',k); 
    eval(sprintf('%s(:,:,slice_select) = false(size(%s,1),size(%s,2));',roiname,roiname,roiname)); 
    eval(sprintf('%s_par{slice_select} = [];',roiname))
end 
set(handles.text24,'String','Remember to Clear All for Deleted Additional ROI to Prevent Unwanted Saving')
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

% --- Executes on button press in pushbutton23. (Disp. Value)
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.mroi == 2 
    roi_select = round(get(handles.roi_selection_slider,'Value')); 
else 
    roi_select = 1; 
end 
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
if roi_select == 1
if ~isempty(handles.ROI_3d)
    if handles.extractiontype ~= 5
        temp = handles.img_ex(:,:,:,array_select,array_select2);
        temp = temp(handles.ROI_3d);
        val = mean(temp);
        sd = std(temp);
        set(handles.text23,'String',sprintf('%.4g (%.4g)',val,sd))
    end 
end
else % Display Values for Additional ROIs if Selected 
if eval(sprintf('~isempty(handles.ROI%d_3d)',roi_select))
    if handles.extractiontype ~= 5
        temp = handles.img_ex(:,:,:,array_select,array_select2);
        eval(sprintf('temp = temp(handles.ROI%d_3d);',roi_select))
        val = mean(temp);
        sd = std(temp);
        set(handles.text23,'String',sprintf('%.4g (%.4g)',val,sd))
    end 
end
end 
clear temp val sd

% --- Executes on button press in pushbutton24. (Save ROI)
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scfile = strtrim(get(handles.edit8,'String')); %subject file  
svname = strtrim(get(handles.edit7,'String')); %ROI File 
pvals = handles.ROI_3d_par;
eval(sprintf('%s = handles.ROI_3d;',svname));
if eval(sprintf('isfile(''%s.mat'');',scfile)) %already existing subject file
    eval(sprintf('save(''%s.mat'',''%s'',''pvals'');',svname,svname));
    eval(sprintf('%s = load(''%s.mat'');',svname,svname)); 
    eval(sprintf('save(''%s.mat'',''%s'',''-append'');',scfile,svname)); 
else %new subject file
    eval(sprintf('save(''%s.mat'',''%s'',''pvals'');',svname,svname));
    eval(sprintf('%s = load(''%s.mat'');',svname,svname));
    eval(sprintf('save(''%s.mat'',''%s'');',scfile,svname));
end 
eval(sprintf('delete(''%s.mat'');',svname)); 
clear(svname)
clear svname pvals scfile 
%Save Additional ROIs if present/~isempty 
for k = 2:6 
roiname = sprintf('handles.ROI%d_3d',k); 
if eval(sprintf('~isempty(%s)',roiname)) && eval(sprintf('~isempty(%s_par)',roiname))
    scfile = strtrim(get(handles.edit8,'String'));
    svname = inputdlg(sprintf('Enter ROI%d Name:',k),'Name the ROI');
    svname = char(strtrim(svname)); 
    eval(sprintf('pvals = %s_par;',roiname))
    eval(sprintf('%s = %s;',svname,roiname));
        if eval(sprintf('isfile(''%s.mat'');',scfile)) %already existing subject file
            eval(sprintf('save(''%s.mat'',''%s'',''pvals'');',svname,svname));
            eval(sprintf('%s = load(''%s.mat'');',svname,svname)); 
            eval(sprintf('save(''%s.mat'',''%s'',''-append'');',scfile,svname)); 
        else %new subject file
            eval(sprintf('save(''%s.mat'',''%s'',''pvals'');',svname,svname));
            eval(sprintf('%s = load(''%s.mat'');',svname,svname));
            eval(sprintf('save(''%s.mat'',''%s'');',scfile,svname));
        end 
    eval(sprintf('delete(''%s.mat'');',svname)); 
    clear(svname)
    clear svname pvals scfile 
end 
end 

function edit7_Callback(hObject, eventdata, handles) %(ROI name)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles) % Change the ROI name  
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton25 (Clear All).
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROI_3d = false(size(handles.ROI_3d));
handles.ROI_3d_par = cell(size(handles.ROI_3d,3));
%Clear Additional ROIs (All Slices) 
for k = 2:6 
    roiname = sprintf('handles.ROI%d_3d',k); 
    eval(sprintf('%s = false(size(%s));',roiname,roiname)); 
    eval(sprintf('%s_par = cell(size(%s,3));',roiname,roiname))
    eval(sprintf('%s = [];',roiname)) 
    eval(sprintf('%s_par = [];',roiname)) 
end 
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

% --- Executes on button press in pushbutton26. (Load ROI)
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scfile = get(handles.edit8,'String'); % Scan File (.mat file) 
svname = get(handles.edit7,'String'); % ROI name (Struct) 
load(sprintf('%s.mat',scfile))
if handles.mroi == 2 
    roi_select = round(get(handles.roi_selection_slider,'Value')); 
else 
    roi_select = 1; 
end 
if roi_select  == 1
    handles.ROI_3d = eval(sprintf('%s.%s',svname,svname));
    handles.ROI_3d_par = eval(sprintf('%s.pvals',svname));
else 
    eval(sprintf('handles.ROI%d_3d = %s.%s;',roi_select,svname,svname))
    eval(sprintf('handles.ROI%d_3d_par = %s.pvals;',roi_select,svname))
end 
clear(svname)
clear scfile svname pvals
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

% --- Executes on button press in pushbutton27. (Lesion)
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.imgopt == 2 
    set(handles.text24,'String','Click on an Image to Select Axes for Lesion')
    w = waitforbuttonpress; 
    if w == 0 
        set(handles.text24,'String','Click on the Brightest part of lesion then nearby healthly tissue') 
    else
        errordlg('Axes Not Selected!','Incorrect Selection','modal')
    end 
end 
slice_select = round(get(handles.slider1,'Value'));
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
img = handles.img_ex(:,:,slice_select,array_select,array_select2);
% First user clicks on the lesion and a reference healthy tissue point
set(handles.text24,'String','Click on the brightest part of lesion, then nearby healthy tissue.')
[y,x,button] = ginput(2);
x = round(x);
y = round(y);
refval = (img(x(1),y(1)) + img(x(2),y(2)))/2;
% Search through neighboring points in a raster fashion
mask = false(size(img));
mask(x(1),y(1)) = true;
ntotal = 1;
newpnts_nxt = [x(1),y(1)];
while ntotal > 0
    newmsk = false(size(img));
    newpnts = [];
    n3 = 1;
    for n0 = 1:ntotal
        for n1 = newpnts_nxt(n0,1)-1:newpnts_nxt(n0,1)+1
            for n2 = newpnts_nxt(n0,2)-1:newpnts_nxt(n0,2)+1
                if ~newmsk(n1,n2) && ~mask(n1,n2) && (img(n1,n2) >= refval)
                    newmsk(n1,n2) = true;
                    newpnts(n3,:) = [n1,n2];
                    n3 = n3+1;
                end
            end
        end
    end
    ntotal = n3-1;
    newpnts_nxt = newpnts;
    mask = or(mask,newmsk);
end
% Add to the overall 3D ROI
pvals = polysort(mask2poly(mask));
sz = size(handles.img_ex);
%Assign Lesion Mask to a ROI Handles (ROI1 default) 
if handles.mroi == 2
    roi_select = round(get(handles.roi_selection_slider,'Value')); 
else 
    roi_select = 1; %Default
end 
if roi_select == 1
if ~isequal(sz(1:3),size(handles.ROI_3d))
    handles.ROI_3d = false(sz(1:3));
    handles.ROI_3d_par = cell(sz(3),1);
end
next = length(handles.ROI_3d_par{slice_select}) + 1;
handles.ROI_3d_par{slice_select}{next} = pvals;
handles.ROI_3d(:,:,slice_select) = or(mask,handles.ROI_3d(:,:,slice_select));
else %if Multiple ROI is Selected 
    roiname = sprintf('handles.ROI%d_3d',roi_select); 
    if eval(sprintf('~isequal(sz(1:3),size(%s))',roiname))
       eval(sprintf(' %s= false(sz(1:3));',roiname))
        eval(sprintf('%s_par = cell(sz(3),1);',roiname)) 
    end
        eval(sprintf('next = length(%s_par{slice_select}) + 1;',roiname)) 
        eval(sprintf('%s_par{slice_select}{next} = pvals;',roiname))
        eval(sprintf(' %s(:,:,slice_select) = or(mask,%s(:,:,slice_select));',roiname,roiname))
end 
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_ROI.
function pushbutton_ROI_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton22.
function pushbutton22_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton21.
function pushbutton21_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ROI1_Callback(hObject, eventdata, handles) % ROI 1 in Menu Bar  
% hObject    handle to ROI1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function defROI1_Callback(hObject, eventdata, handles) % Define ROI 1 
% hObject    handle to defROI1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
handles.ROI = roipoly;
% Update handles structure
guidata(hObject, handles);

function saveROI1_Callback(hObject, eventdata, handles) % save ROI 1 (2D) 
% hObject    handle to saveROI1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = handles.ROI;
save roi roi
clear roi

function loadROI1_Callback(hObject, eventdata, handles) %load 2D-ROI 
% hObject    handle to loadROI1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call loadfile dialog
user_response = loadfile_gui('Title','Load ROI');
if ~isempty(user_response)
    handles.ROI = user_response;
end
% Update handles structure
guidata(hObject, handles);

function defROI2_Callback(hObject, eventdata, handles) %Define 2D ROI 2 
% hObject    handle to defROI2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROI2 = roipoly;
% Update handles structure
guidata(hObject, handles);

function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask = handles.ROI;
% Update handles structure
guidata(hObject, handles);
showimage(hObject,handles)

function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imsize = size(handles.img_ex);
handles.mask = ones(imsize(1),imsize(2));
% Update handles structure
guidata(hObject, handles);
showimage(hObject,handles)

function Untitled_11_Callback(hObject, eventdata, handles) %2D DTI 
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.ROI)
    errordlg('Select an ROI first','Incorrect Selection','modal')
else
    slice_select = round(get(handles.slider1,'Value'));
    array_select = round(get(handles.slider2,'Value'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    temp = handles.img_ex(:,:,slice_select,array_select,array_select2);
    ave = mean(temp(handles.ROI));
%     ave = squeeze(sum(sum(handles.img_ex(:,:,slice_select,array_select,array_select2).*...
%         handles.ROI,1),2)) / numelements;
    display(ave)
    clear temp ave
end

function Untitled_12_Callback(hObject, eventdata, handles) %Average for 2D ROI 
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.ROI)
    errordlg('Select an ROI first','Incorrect Selection','modal')
elseif (size(handles.img,1) ~= 6) && (size(handles.img,1) ~= 3)
    errordlg('You must use DTI data','Incorrect Selection','modal')
else
    slice_select = round(get(handles.slider1,'Value'));
    array_select = round(get(handles.slider2,'Value'));
    array_select2 = round(get(handles.slider5,'Value'));
    numelements  = length(find(handles.ROI));
    ave = squeeze(sum(sum(handles.img(:,:,:,slice_select,array_select,array_select2).*...
        repmat(reshape(handles.ROI,[1 size(handles.ROI,1) size(handles.ROI,2)]),...
        [size(handles.img,1) 1 1]),2),3)) / numelements
end

function Untitled_13_Callback(hObject, eventdata, handles) %2D SNR Calculations 
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SNR = calcSNR(handles);
set(handles.text_SNR,'String',sprintf('%g',SNR))

function output = calcSNR(handles)
% Retrieve relevant parameters
slice_select = round(get(handles.slider1,'Value'));
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));

if isempty(handles.ROI) || isempty(handles.ROI2)
    output = 0;
else
    if handles.extractiontype ~= 5
        temp = handles.img_ex(:,:,slice_select,array_select,array_select2);
        output = mean(temp(handles.ROI))/std(temp(handles.ROI2));
        clear temp
    else
        output = 0;
    end
end

function defROI1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to defROI1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton29. (SNR for 3-D ROI) 
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SNR = calcSNR3D(handles);
set(handles.text_SNR,'String',sprintf('%g',SNR))

function output = calcSNR3D(handles) %Callbacked in line 173 
% Retrieve relevant parameters
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));

if isempty(handles.ROI_3d) || isempty(handles.ROI2_3d)
    output = 0;
else
    if handles.extractiontype ~= 5
        temp = handles.img_ex(:,:,:,array_select,array_select2);
        output = mean(temp(handles.ROI_3d))/std(temp(handles.ROI2_3d));
        clear temp 
    else
        output = 0;
    end
end

function edit8_Callback(hObject, eventdata, handles) %Edit the Scan File Name
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Untitled_16_Callback(hObject, eventdata, handles) % Legacy Load ROI function 
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function LoadROIV2_Callback(hObject, eventdata, handles) %Loads ROI from .mat file with ROIname only
% hObject    handle to LoadROIV2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
svname = get(handles.edit7,'String');
load(sprintf('%s.mat',svname))
handles.ROI_3d = eval(svname);
handles.ROI_3d_par = pvals;
clear(svname)
clear svname pvals
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function saveROIV2_Callback(hObject, eventdata, handles) %Saves ROI into .mat file with ROIname only 
% hObject    handle to saveROIV2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
svname = strtrim(get(handles.edit7,'String'));
pvals = handles.ROI_3d_par;
eval(sprintf('%s = handles.ROI_3d;',svname));
eval(sprintf('save(''%s.mat'',''%s'',''pvals'');',svname,svname));
clear(svname)
clear svname pvals

% --- Executes when selected object is changed in uibuttongroup3. %image
% options
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_single' %Disable Functions and Set back to ROI 
        handles.imgopt = 1;
        for i = 2:4 
            eval(sprintf('handles.img%d = [];',i)) 
            eval(sprintf('handles.img_ex%d = [];',i)) 
        end 
        arrayfun(@cla, findall(0,'type','axes'))
        set(handles.listbox2,'Enable','off','Visible','off')
        set(handles.text32,'Visible','off') 
        set(handles.edit9,'Visible','off','Enable','off') 
        set(handles.text33,'Visible','off')
        set(handles.edit_max2,'Enable','off','Visible','off') 
        set(handles.reset_button2, 'Enable','off', 'Visible', 'off')
        set(handles.reset_image_3, 'Visible', 'off')
        set(handles.text_minmax3,'Visible','off')
        set(handles.reset_image_4, 'Visible', 'off')
        set(handles.text_minmax4,'Visible','off')
        set(handles.multi_image_tools,'Visible','off')
    case 'radiobutton_Multiple' %Enable Functions to add img2 to plots
        handles.imgopt = 2; 
        update_listbox(handles)
        set(handles.listbox2,'Enable','on','Visible','on')
        set(handles.text32,'Visible','on') 
        set(handles.edit9,'Visible','on','Enable','on') 
        set(handles.text33,'Visible','on')
        set(handles.edit_max2,'Enable','on','Visible','on') 
        set(handles.reset_button2, 'Enable','on', 'Visible', 'on')
        set(handles.text24,'String','Select a Second Variable')
    case 'radiobutton_overlap' %Enable Functions to add img2 to plots
        handles.imgopt = 3;
        arrayfun(@cla, findall(0,'type','axes'))
        update_listbox(handles)
        set(handles.listbox2,'Enable','on','Visible','on')
        set(handles.text32,'Visible','on') 
        set(handles.edit9,'Visible','on','Enable','on')
        set(handles.text33,'Visible','on')
        set(handles.edit_max2,'Enable','on','Visible','on') 
        set(handles.reset_button2, 'Enable','on', 'Visible', 'on')
        set(handles.text24,'String','Select a Second Variable')
        set(handles.reset_image_3, 'Visible', 'off')
        set(handles.text_minmax3,'Visible','off')
        set(handles.reset_image_4, 'Visible', 'off')
        set(handles.text_minmax4,'Visible','off')
        set(handles.multi_image_tools,'Visible','off')
end
% Update handles structure
guidata(hObject, handles);

function edit9_Callback(hObject, eventdata, handles)%Min2 
% hObject    handle to edit_max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_max2 as text
%        str2double(get(hObject,'String')) returns contents of edit_max2 as a double
handles.minvalue2 = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Show image
showimage(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles) % Min2 
% hObject    handle to edit_max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_max2_Callback(hObject, eventdata, handles) % Max2 Callback Function 
% hObject    handle to edit_max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_max2 as text
%        str2double(get(hObject,'String')) returns contents of edit_max2 as a double
handles.maxvalue2 = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Show image
showimage(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_max2_CreateFcn(hObject, eventdata, handles) %Max2 Create Function 
% hObject    handle to edit_max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clear_roi_1_Callback(hObject, eventdata, handles) % Clear ROI1 only 
% hObject    handle to clear_roi_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slice_select = round(get(handles.slider1,'Value'));
handles.ROI_3d(:,:,slice_select) = false(size(handles.ROI_3d,1),size(handles.ROI_3d,2));
handles.ROI_3d_par{slice_select} = [];
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function clear_roi_2_Callback(hObject, eventdata, handles) % Clear ROI2 Only
% hObject    handle to clear_roi_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slice_select = round(get(handles.slider1,'Value'));
handles.ROI2_3d(:,:,slice_select) = false(size(handles.ROI2_3d,1), size(handles.ROI2_3d,2)); 
handles.ROI2_3d_par{slice_select} = []; 
set(handles.text24,'String','Remember to Clear All for Deleted Additional ROI to Prevent Unwanted Saving')
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function current_roi_Callback(hObject, eventdata, handles) %Clear Current ROI 
% hObject    handle to current_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_select = round(get(handles.roi_selection_slider,'Value')); 
slice_select = round(get(handles.slider1,'Value'));
if roi_select == 1
    handles.ROI_3d(:,:,slice_select) = false(size(handles.ROI_3d,1),size(handles.ROI_3d,2));
    handles.ROI_3d_par{slice_select} = [];
else 
    roiname = sprintf('handles.ROI%d_3d',roi_select); 
    eval(sprintf('%s(:,:,slice_select) = false(size(%s,1), size(%s,2));',roiname,roiname,roiname)); 
    eval(sprintf('%s_par{slice_select} = [];',roiname)) 
end 
set(handles.text24,'String','Remember to Clear All for Deleted Additional ROI to Prevent Unwanted Saving')
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function Clear_slice_Callback(hObject, eventdata, handles) % Clear Slice Context Menu
% hObject    handle to Clear_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function save_image_Callback(hObject, eventdata, handles) %Save image menu bar 
% hObject    handle to save_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function file_type_Callback(hObject, eventdata, handles)
% hObject    handle to file_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function save_jpg_Callback(hObject, eventdata, handles) % jpg 
% hObject    handle to save_jpg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.savetype = 1;
% Update handles structure
guidata(hObject, handles);
%callback save function 
save_image(handles); 

function save_image(handles) %function to callback in all save_images 
%Prompt Image Name 
imagename = inputdlg('Enter Image Name:','Save Image'); 
imagename = char(strtrim(imagename)); 
% Set interpolation factor
factor = handles.interpfactor;
% Retrieve relevant parameters
slice_select = round(get(handles.slider1,'Value'));
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
% wnd = get(handles.slider3,'Value');
% lvl = get(handles.slider4,'Value');
% Apply window and level
if handles.extractiontype == 5
    imsize = size(squeeze(handles.img_ex(1,:,:,slice_select,array_select,array_select2)));
    if min( imsize == size(handles.mask) )
        img = permute( handles.img_ex(:,:,:,slice_select,array_select,array_select2), [2 3 1] ).*...
            repmat(handles.mask, [1 1 3] );
    else
        img = permute( handles.img_ex(:,:,:,slice_select,array_select,array_select2), [2 3 1] );
    end
%     % Interpolate image
%     if factor>1
%         [X,Y] = meshgrid(1:1:size(imagedata,1));
%         [XI,YI] = meshgrid(1:1/factor:size(imagedata,1));
%         imagedata_a = zeros(length(XI),length(YI),3);
%         for n=1:3
%             imagedata_a(:,:,n) = interp2(X,Y,imagedata(:,:,n),XI,YI,'spline');
%         end
%         imagedata = imagedata_a;
%     end
else
    im_max = max(max(max(handles.img_ex(:,:,:,array_select,array_select2))));
    im_min = min(min(min(handles.img_ex(:,:,:,array_select,array_select2))));
%     maxlevel = 0.5*wnd*(im_max-im_min) + lvl*(im_max+im_min);
%     minlevel = 2*lvl*(im_max+im_min)-maxlevel;
    % Use min and max fields instead
    maxlevel = handles.maxvalue;
    minlevel = handles.minvalue;

    % Apply mask
    if min( size(handles.img_ex(:,:,slice_select,array_select,array_select2)) == size(handles.mask) )
        img = handles.img_ex(:,:,slice_select,array_select,array_select2).*handles.mask;
    else
        img = handles.img_ex(:,:,slice_select,array_select,array_select2);
    end
    % Scale data
    maxthresh = img <= maxlevel;
    minthresh = img >= minlevel;
    img = ((img.*maxthresh.*minthresh + (1-maxthresh)*maxlevel +...
        (1-minthresh)*minlevel) - minlevel) / (maxlevel-minlevel);
    % Apply the aspectrato
    if handles.aspecttype == 1
        % Here handles.aspectratio gives the in-plane aspect ratio of the
        % voxels
        mfact = handles.aspectratio;
    else
        % Here handles.aspectratio gives the in-plane aspect ratio of the
        % FOV
        mfact = handles.aspectratio*size(handles.img_ex,2)/size(handles.img_ex,1);
    end
    if mfact < 1
        newmatrix = [size(handles.img_ex,1) round(size(handles.img_ex,2)/mfact)];
    else
        newmatrix = [round(size(handles.img_ex,1)*mfact) size(handles.img_ex,2)];
    end
    img = imresize(img,newmatrix,'nearest');
end
% Convert the image to uint8
img = uint8(round(img*255));
% Upsample image to be nearly 256 pixels wide (better viewing on computers)
targetres = 256;
if round(targetres/size(img,2)) > 1
    img = imresize(img, round(targetres/size(img,2)),'nearest');
end
% Use colormap when saving when real flag or colourmap is selected
if handles.extractiontype == 7 || handles.colour == 2 
  usecolormap = 1;
else
  usecolormap = 0;
end
% Save the image
if usecolormap == 1
    map = colormap(handles.cmap);
%     map(1,:) = 0;
    if handles.savetype == 1
        imwrite(img,map,sprintf('%s.jpg',imagename),'jpg','Quality',100)
    elseif handles.savetype == 2
        imwrite(img,map,sprintf('%s.tif',imagename),'tif')
    elseif handles.savetype == 3
        imwrite(img,map,sprintf('%s.bmp',imagename),'bmp')
    end
else
    if handles.savetype == 1
        imwrite(img,sprintf('%s.jpg',imagename),'jpg','Quality',100)
    elseif handles.savetype == 2
        imwrite(img,sprintf('%s.tif',imagename),'tif')
    elseif handles.savetype == 3
        imwrite(img,sprintf('%s.bmp',imagename),'bmp')
    end
end
clear imagedata img

function save_tif_Callback(hObject, eventdata, handles)
% hObject    handle to save_tif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.savetype = 2; 
% Update handles structure
guidata(hObject, handles);
%callback save function 
save_image(handles); 

function save_bmp_Callback(hObject, eventdata, handles)
% hObject    handle to save_bmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.savetype = 3;
% Update handles structure
guidata(hObject, handles);
%callback save function 
save_image(handles) 

% --- Executes on button press in reset_button2. (%Reset Button) 
function reset_button2_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,~,minval2,maxval2] = resetscale(handles); 
handles.maxvalue2 =  maxval2;  
handles.minvalue2 = minval2; 
set(handles.edit9,'String', handles.minvalue2)
set(handles.edit_max2,'String', handles.maxvalue2)
%Show the Image 
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function View_options_Callback(hObject, eventdata, handles)
% hObject    handle to View_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function axial_view_Callback(hObject, eventdata, handles) %Axial View 
% hObject    handle to axial_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.view = 1;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
if eval(sprintf('~isempty(handles.img%d)',i))
    eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function sagittal_view_Callback(hObject, eventdata, handles) % Sagittal View 
% hObject    handle to sagittal_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.view = 2;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
if eval(sprintf('~isempty(handles.img%d)',i))
    eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function coronal_view_Callback(hObject, eventdata, handles) %Coronal View 
% hObject    handle to coronal_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.view = 3;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function data_type_Callback(hObject, eventdata, handles) %MenuBar --> List Data Types 
% hObject    handle to data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function abs_type_Callback(hObject, eventdata, handles)
% hObject    handle to abs_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 1;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function phase_type_Callback(hObject, eventdata, handles)
% hObject    handle to phase_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.extractiontype = 2;
 handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function real_type_Callback(hObject, eventdata, handles)
% hObject    handle to real_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 7;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function imag_type_Callback(hObject, eventdata, handles)
% hObject    handle to imag_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 8;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end  
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function trace_type_Callback(hObject, eventdata, handles)
% hObject    handle to trace_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 3;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end  
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function FA_type_Callback(hObject, eventdata, handles)
% hObject    handle to FA_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 4;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function orient_type_Callback(hObject, eventdata, handles) %(Orientation)
% hObject    handle to orient_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 5;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end  
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function log_type_Callback(hObject, eventdata, handles)
% hObject    handle to log_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.extractiontype = 6;
handles.img_ex = extraction(handles.slider5,handles.img,handles.extractiontype,handles.view);
for i = 2:4 
    if eval(sprintf('~isempty(handles.img%d)',i))
        eval(sprintf('handles.img_ex%d = extraction(handles.slider5, handles.img%d,handles.extractiontype,handles.view);',i,i))  
    end 
end 
% Show the image
adjustsliders(handles)
showimage(hObject,handles)
% Update handles structure
guidata(hObject, handles);

function colour_options_Callback(hObject, eventdata, handles) 
% hObject    handle to colour_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function grayscale_opt_Callback(hObject, eventdata, handles) % Grayscale 
% hObject    handle to grayscale_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colour = 1;
% Update handles structure
guidata(hObject, handles);
%show image 
showimage(hObject,handles)

function colour_opt_Callback(hObject, eventdata, handles) %Parula Colourmap
% hObject    handle to colour_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colour = 2;
% Update handles structure
guidata(hObject, handles);
%show image 
showimage(hObject,handles)

function clear_all_roi_Callback(hObject, eventdata, handles) %Clears All ROI1 only 
% hObject    handle to clear_all_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROI_3d = false(size(handles.ROI_3d));
handles.ROI_3d_par = cell(size(handles.ROI_3d,3));
handles.ROI_3d = [];
handles.ROI_3d_par = [];
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function clear_all_roi2_Callback(hObject, eventdata, handles) %Clears All ROI2 only 
% hObject    handle to clear_all_roi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROI2_3d = false(size(handles.ROI2_3d));
handles.ROI2_3d_par = cell(size(handles.ROI2_3d,3));
handles.ROI2_3d = [];
handles.ROI2_3d_par = [];
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function clear_all_current_Callback(hObject, eventdata, handles) %All Selected ROI Only 
% hObject    handle to clear_all_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_select = round(get(handles.roi_selection_slider,'Value')); 
if roi_select == 1 
    handles.ROI_3d = false(size(handles.ROI_3d));
    handles.ROI_3d_par = cell(size(handles.ROI_3d,3));
else
    eval(sprintf('handles.ROI%d_3d = false(size(handles.ROI%d_3d));',roi_select,roi_select))
    eval(sprintf('handles.ROI%d_3d_par = cell(size(handles.ROI%d_3d,3));', roi_select,roi_select))
    eval(sprintf('handles.ROI%d_3d = [];',roi_select))
    eval(sprintf('handles.ROI%d_3d_par = [];',roi_select))
end 
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function Clear_all_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function disp_val_img2_Callback(hObject, eventdata, handles) %Displays the end of the Image
% hObject    handle to disp_val_img2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.mroi == 2 
    roi_select = round(get(handles.roi_selection_slider,'Value')); 
else 
    roi_select = 1; 
end 
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
if ~isempty(handles.img2)
if roi_select == 1
if ~isempty(handles.ROI_3d)
    if handles.extractiontype ~= 5
        temp = handles.img_ex2(:,:,:,array_select,array_select2);
        temp = temp(handles.ROI_3d);
        val = mean(temp);
        sd = std(temp);
        set(handles.text23,'String',sprintf('%.4g (%.4g)',val,sd))
    end
end
else
  if eval(sprintf('~isempty(handles.ROI%d_3d)',roi_select))
    if handles.extractiontype ~= 5
        temp = handles.img_ex2(:,:,:,array_select,array_select2);
        eval(sprintf('temp = temp(handles.ROI%d_3d);',roi_select))
        val = mean(temp);
        sd = std(temp);
        set(handles.text23,'String',sprintf('%.4g (%.4g)',val,sd))
    end
  end
end 
clear temp val sd
else 
    errordlg('Second Variable Must be Selected','Incorrect Selection','modal') 
end 
%Calculations for Image 2
function Disp_val2_Callback(hObject, eventdata, handles) %context menu link 
% hObject    handle to Disp_val2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function SNR_img_2_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_img_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.img2)
SNR = calcSNR3D2(handles);
set(handles.text_SNR,'String',sprintf('%g',SNR))
else
    errordlg('Must Select Second Variable','Incorrect Selection','modal')
end 
function output = calcSNR3D2(handles)
% Retrieve relevant parameters
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));

if isempty(handles.ROI_3d) || isempty(handles.ROI2_3d)
    output = 0;
else
    if handles.extractiontype ~= 5
        temp = handles.img_ex2(:,:,:,array_select,array_select2); %change to img_ex2 
        output = mean(temp(handles.ROI_3d))/std(temp(handles.ROI2_3d));
        clear temp 
    else
        output = 0;
    end
end

function multi_ROI_add_Callback(hObject, eventdata, handles)
% hObject    handle to multi_ROI_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mroi = 2; %Enable Multiple ROI Functions on for all related Pushbuttons/Context Menus
set(handles.roi_selection_slider, 'Enable','on','Visible','on') 
set(handles.roi_selection_display,'Visible','on')
set(handles.current_roi,'Visible','on')
set(handles.clear_all_current, 'Visible','on')
set(handles.text24,'String','ROI1:Yellow, ROI2:Magenta, ROI3:Green, ROI4:Red, ROI5:Blue, ROI6:Cyan')
%update handles information
guidata(hObject, handles);

function Multiple_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Multiple_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Single_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Single_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mroi = 1; %default (single ROI) 
set(handles.roi_selection_slider, 'Enable','off','Visible','off')
set(handles.roi_selection_display,'Visible','off')
set(handles.current_roi,'Visible','off')
set(handles.clear_all_current, 'Visible','off')
set(handles.text24,'String','Info')
guidata(hObject, handles);

% --- Executes on slider movement.
function roi_selection_slider_Callback(hObject, eventdata, handles)
% hObject    handle to roi_selection_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.roi_selection_display,'String',sprintf('ROI %d\n',round(get(handles.roi_selection_slider,'Value'))))
% --- Executes during object creation, after setting all properties.
function roi_selection_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_selection_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function multi_image_tools_Callback(hObject, eventdata, handles)
% hObject    handle to multi_image_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function reset_image_3_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function reset_image_4_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function reset_image4_manual_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image4_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minimax = inputdlg({'Enter Min value', 'Enter Max Value'},'Input Min and Max Values',...
    [1 16; 1 16],{sprintf('%d',handles.minvalue4), sprintf('%d',handles.maxvalue4)}); 
handles.minvalue4 = str2double(char(minimax{1})); 
handles.maxvalue4 = str2double(char(minimax{2})); 
set(handles.text_minmax4, 'String',sprintf('Min4: %d Max4: %d',handles.minvalue4,handles.maxvalue4))
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function reset_image4_auto_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image4_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
handles.maxvalue4= max(max(max(handles.img_ex4(:,:,:,array_select,array_select2))));
handles.minvalue4 = min(min(min(handles.img_ex4(:,:,:,array_select,array_select2))));
set(handles.text_minmax4, 'String',sprintf('Min4: %d Max4: %d',handles.minvalue4,handles.maxvalue4))
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function reset_image_3_manual_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image_3_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minimax = inputdlg({'Enter Min value', 'Enter Max Value'},'Input Min and Max Values', ...
    [1 16; 1 16],{sprintf('%d',handles.minvalue3), sprintf('%d',handles.maxvalue3)}); 
handles.minvalue3 = str2double(char(minimax{1})); 
handles.maxvalue3 = str2double(char(minimax{2}));
set(handles.text_minmax3, 'String',sprintf('Min3: %d Max3: %d',handles.minvalue3,handles.maxvalue3))
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)

function reset_img3_auto_Callback(hObject, eventdata, handles)
% hObject    handle to reset_img3_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
array_select = round(get(handles.slider2,'Value'));
array_select2 = round(get(handles.slider5,'Value'));
handles.maxvalue3 = max(max(max(handles.img_ex3(:,:,:,array_select,array_select2))));
handles.minvalue3 = min(min(min(handles.img_ex3(:,:,:,array_select,array_select2))));
set(handles.text_minmax3, 'String',sprintf('Min3: %d Max3: %d',handles.minvalue3,handles.maxvalue3))
% Update handles structure
guidata(hObject, handles);
% Show the image
showimage(hObject,handles)
