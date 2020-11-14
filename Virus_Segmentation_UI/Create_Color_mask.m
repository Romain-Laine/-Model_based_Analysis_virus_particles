function varargout = Create_Color_mask(varargin)
% CREATE_COLOR_MASK MATLAB code for Create_Color_mask.fig
%      CREATE_COLOR_MASK, by itself, creates a new CREATE_COLOR_MASK or raises the existing
%      singleton*.
%
%      H = CREATE_COLOR_MASK returns the handle to a new CREATE_COLOR_MASK or the handle to
%      the existing singleton*.
%
%      CREATE_COLOR_MASK('CA LLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATE_COLOR_MASK.M with the given input arguments.
%
%      CREATE_COLOR_MASK('Property','Value',...) creates a new CREATE_COLOR_MASK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Create_Color_mask_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Create_Color_mask_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Create_Color_mask

% Last Modified by GUIDE v2.5 12-Jun-2014 11:41:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Create_Color_mask_OpeningFcn, ...
                   'gui_OutputFcn',  @Create_Color_mask_OutputFcn, ...
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


% --- Executes just before Create_Color_mask is made visible.
function Create_Color_mask_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Create_Color_mask (see VARARGIN)
handles.output = hObject;

set(handles.Saturation_min, 'String', '0');
set(handles.Saturation_max, 'String', '0.1');
set(handles.Sensitivity, 'String', '1');

set(handles.Gamma, 'String', '0.5');
set(handles.BG_remove, 'Value', 1);
set(handles.Invert_mask, 'Value', 0);

set(handles.x_shift, 'String', 0);
set(handles.y_shift, 'String', 0);

set(handles.n_dilate, 'String', 20);
handles.n_dilate = 20;

% For the old camera system CamPixelSize = 158.4;
% 157.2 nm / pixel was obtained from measurement
set(handles.CamPixelSize, 'String', '157.2');
handles.CamPixelSize = 157.2;

handles.x_shift = 0;
handles.y_shift = 0;

handles.Saturation_min = 0;
handles.Saturation_max = 0.1;
handles.Gamma = 0.5;
handles.Sensitivity = 1;

handles.BG_remove = 1;
handles.Invert_mask = 0;

StartPath = varargin{1}(1);
ColorImage_in = varargin{1}(2);
ColorImage_in = cell2mat(ColorImage_in);

GreyImage_in = varargin{1}(3);
GreyImage_in = cell2mat(GreyImage_in);

PixelSize = varargin{1}(4);
PixelSize = cell2mat(PixelSize);

[FileName,PathName] = uigetfile(strcat(StartPath,'*.tif'),'Select a file :');
disp(['Mask loaded: ',PathName,FileName]);

ImInfo = imfinfo([PathName,FileName],'tif');     % Extract file headers and info
n_files = numel(ImInfo);         % Number of images in the tif

x_res = ImInfo(1).Width;
y_res = ImInfo(1).Height;

SumIm = double(zeros(x_res,y_res));

if n_files > 1
    h = waitbar(0,'Please wait while it calculates the average image ...');
    for i = 1:n_files
        waitbar(i/n_files)
        SumIm = SumIm + double(imread([PathName,FileName],'tif',i,'Info',ImInfo));
    end
    close(h);
else
    SumIm = double(imread([PathName,FileName],'tif'));
end

SumIm = SumIm/n_files;

figure('Color','white','name','Mask creation - Intensity image');
imshow(SumIm,[]);
title 'Intensity image'

handles.SumIm = SumIm;
handles.PathName = PathName;
handles.FileName = FileName;
handles.hFig1 = '';
handles.hFig2 = '';

handles.ColorImage_in = ColorImage_in;
handles.GreyImage_in = GreyImage_in;

handles.ColorImage_out = ColorImage_in;
handles.GreyImage_out = GreyImage_in;

handles.PixelSize = PixelSize;

% Update handles structure
guidata(hObject, handles);
uiwait(hObject);


% UIWAIT makes Create_Color_mask wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Create_Color_mask_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% ColorImage_out = handles.ColorImage_out;
% GreyImage_out = handles.GreyImage_out;
% 
% size(cat(3,ColorImage_out, GreyImage_out))
% varargout{1} = cat(3,ColorImage_out, GreyImage_out);
Mask_resized = handles.Mask_resized;
varargout{1} = Mask_resized;

close(hObject);


% --- Executes on button press in BG_remove.
function BG_remove_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to BG_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BG_remove
BG_remove = get(hObject, 'Value');

handles.BG_remove = BG_remove;
guidata(hObject,handles)


% --- Executes on button press in Invert_mask.
function Invert_mask_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Invert_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Invert_mask
Invert_mask = get(hObject, 'Value');

% Save the new density value
handles.Invert_mask = Invert_mask;
guidata(hObject,handles)




% --- Executes on button press in Create_mask.
function Create_mask_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Create_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Gaussian filter %

BG_remove = handles.BG_remove;

Invert_mask = handles.Invert_mask;
SumIm = handles.SumIm;
Gamma = handles.Gamma;
hFig1 = handles.hFig1;
hFig2 = handles.hFig2;
PixelSize = handles.PixelSize;
CamPixelSize = handles.CamPixelSize;    % Camera pixel size (nm)

x_shift = handles.x_shift;
y_shift = handles.y_shift;

PathName = handles.PathName;
FileName = handles.FileName;

Saturation_min = handles.Saturation_min;
Saturation_max = handles.Saturation_max;
Sensitivity = handles.Sensitivity;

ColorImage_in = handles.ColorImage_in;
% assignin('base','SumIm',SumIm);

%% Gaussian filter
hsize = 6;
sigma = 1;
filter = fspecial('gaussian', hsize, sigma);
I_n = imfilter(SumIm, filter,'replicate','same');



%% Average filter %%
disp('Remove background?');
disp(BG_remove);

if BG_remove == 1;
    hsize = 50;
    filter = fspecial('average', hsize);
    BG_image = imfilter(I_n, filter,'replicate','same');
    I_n = I_n - BG_image;
end


% Normalise to 0-1
I_n = (I_n-min(I_n(:)))/(max(I_n(:)) - min(I_n(:)));

% Adjust the levels
I = imadjust(I_n,[Saturation_min Saturation_max],[],Gamma);         


% Convert to mask
Level = graythresh(I);
Level = max(Level/Sensitivity,0);
Level = min(Level/Sensitivity,1);
disp(['Threshold level: ',num2str(Level)])
Mask = im2bw(I,Level);
Mask = bwmorph(Mask,'clean');


% Creating edges and Color mask
Mask_edges = bwmorph(Mask,'remove');
ImColor = cat(3,I_n,I_n,I_n + Mask_edges);

if isempty(hFig1) == 0
    close(hFig1);
end

hFig1 = figure('Color','white');
imshow(ImColor)
set(hFig1,'name',[PathName,FileName]);


% Resizing the mask 

% This resizing does not take into account the borders introduced by
% rapidSTORM --> just rescaling the image tyo fit that of the SR image

% numrows = size(ColorImage_in,1);
% numcols = size(ColorImage_in,2);
% Mask_resized =  imresize(Mask, [numrows numcols],'bicubic');


% This resizing takes into account the borders added by rapidSTORM
% However it needs to know the SR pixel size and the camera pixel size

numrows = round(size(SumIm,1)*CamPixelSize/PixelSize);
numcols = round(size(SumIm,2)*CamPixelSize/PixelSize);

Mask_resized = zeros(size(ColorImage_in,1),size(ColorImage_in,2));
Mask_resized(1:numrows, 1:numcols) = imresize(Mask, [numrows numcols],'bicubic');


% Smoothing the edges of the mask through opening 
n_smoothing = 10;
se = strel('disk',max(1,round(abs(n_smoothing))));
Mask_resized = imopen(Mask_resized, se);


% Shifting
disp(['x shift: ',num2str(x_shift),' - y shift: ',num2str(y_shift)]);
Mask_resized = circshift(Mask_resized, [x_shift y_shift]);

if x_shift > 0;
    Mask_resized(1:x_shift, :) = 0;
end
if x_shift < 0;
    Mask_resized(end+1+x_shift:end, :) = 0;
end

if y_shift > 0;
    Mask_resized(:,1:y_shift) = 0;
end
if y_shift < 0;
    Mask_resized(:,end+1+y_shift:end) = 0;
end


% Dilating the mask
disp(['Dilation factor: ',num2str(handles.n_dilate),' pixels']);
se = strel('disk',max(1,round(abs(handles.n_dilate))));
Mask_resized = imdilate(Mask_resized, se);

Mask_resized_edges = bwmorph(Mask_resized,'remove');

% Make the edge mask thicker - for display purposes
n_dilate_display = 10;
se = strel('disk',max(1,round(abs(n_dilate_display))));
Mask_resized_edges = imdilate(Mask_resized_edges, se);

Blank_image = zeros(size(ColorImage_in,1),size(ColorImage_in,2));

if isempty(hFig2) == 0
    close(hFig2);
end

hFig2 = figure('Color','white');
imshow(ColorImage_in + uint8(256*cat(3,Blank_image,Blank_image, Mask_resized_edges)))
set(hFig2,'name',[PathName,FileName])

%% Calculate the output image %%

% Make Mask negative
if Invert_mask == 1
    Mask_resized = (-1)*(Mask_resized - 1);
end

handles.Mask_resized = Mask_resized;
handles.hFig1 = hFig1;
handles.hFig2 = hFig2;

guidata(hObject,handles);




% --- Executes on button press in Ok_button.
function Ok_button_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to Ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume();



% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mask_resized = handles.Mask_resized;
Mask_resized = ones(size(Mask_resized));
handles.Mask_resized = Mask_resized;

guidata(hObject, handles);
uiresume();


function Saturation_min_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Saturation_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Saturation_min as text
%        str2double(get(hObject,'String')) returns contents of Saturation_min as a double
Saturation_min = str2double(get(hObject, 'String'));
if isnan(Saturation_min)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Saturation_min = Saturation_min;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Saturation_min_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Saturation_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Saturation_max_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Saturation_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Saturation_max as text
%        str2double(get(hObject,'String')) returns contents of Saturation_max as a double
Saturation_max = str2double(get(hObject, 'String'));
if isnan(Saturation_max)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Saturation_max = Saturation_max;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function Saturation_max_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Saturation_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gamma_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gamma as text
%        str2double(get(hObject,'String')) returns contents of Gamma as a double
Gamma = str2double(get(hObject, 'String'));
if isnan(Gamma)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Gamma = Gamma;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Gamma_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sensitivity_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Sensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sensitivity as text
%        str2double(get(hObject,'String')) returns contents of Sensitivity as a double
Sensitivity = str2double(get(hObject, 'String'));
if isnan(Sensitivity)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Sensitivity = Sensitivity;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Sensitivity_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Sensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_shift_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to x_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_shift as text
%        str2double(get(hObject,'String')) returns contents of x_shift as a double
x_shift = str2double(get(hObject, 'String'));
if isnan(x_shift)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.x_shift = x_shift;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function x_shift_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to x_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_shift_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to y_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_shift as text
%        str2double(get(hObject,'String')) returns contents of y_shift as a double
y_shift = str2double(get(hObject, 'String'));
if isnan(y_shift)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.y_shift = y_shift;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function y_shift_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to y_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CamPixelSize_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to CamPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CamPixelSize as text
%        str2double(get(hObject,'String')) returns contents of CamPixelSize as a double
CamPixelSize = str2double(get(hObject, 'String'));
if isnan(CamPixelSize)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.CamPixelSize = CamPixelSize;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function CamPixelSize_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to CamPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_dilate_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to n_dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_dilate as text
%        str2double(get(hObject,'String')) returns contents of n_dilate as a double
n_dilate = str2double(get(hObject, 'String'));
if isnan(n_dilate)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.n_dilate = n_dilate;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function n_dilate_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to n_dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Invert_mask_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to Invert_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function BG_remove_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to BG_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



