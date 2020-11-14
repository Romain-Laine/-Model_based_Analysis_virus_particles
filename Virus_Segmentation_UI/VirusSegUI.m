function varargout = VirusSegUI(varargin)
% VIRUSSEGUI MATLAB code for VirusSegUI.fig
%      VIRUSSEGUI, by itself, creates a new VIRUSSEGUI or raises the existing
%      singleton*.
%
%      H = VIRUSSEGUI returns the handle to a new VIRUSSEGUI or the handle to
%      the existing singleton*.
%
%      VIRUSSEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIRUSSEGUI.M with the given input arguments.
%
%      VIRUSSEGUI('Property','Value',...) creates a new VIRUSSEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VirusSegUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VirusSegUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VirusSegUI

% Last Modified by GUIDE v2.5 28-May-2014 13:16:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VirusSegUI_OpeningFcn, ...
    'gui_OutputFcn',  @VirusSegUI_OutputFcn, ...
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


% --- Executes just before VirusSegUI is made visible.
function VirusSegUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VirusSegUI (see VARARGIN)
% Choose default command line output for VirusSegUI
handles.output = hObject;

set(handles.Min_size, 'String', '50');
set(handles.Max_size, 'String', '20000');
set(handles.Min_ecc, 'String', '0');
set(handles.Max_ecc, 'String', '1');
set(handles.Sensitivity, 'String', '1');
set(handles.PixelSize, 'String', '10');
set(handles.Scale_bar_size, 'String', '100');
set(handles.MaxPanelSize, 'String', '15');
set(handles.ADC_cutoff, 'String', '0');
set(handles.Radius_cutoff, 'String', '150');
set(handles.Frame_min, 'String', '0');
set(handles.Frame_max, 'String', '50000');
set(handles.I_cutoff, 'String', '1');
set(handles.Gamma, 'String', '1');

handles.Gamma = 1;
handles.I_cutoff = 1;
handles.Frame_min = 0;
handles.Frame_max = 50000;
handles.Radius_cutoff = 150;
handles.ADC_cutoff = 0;
handles.MaxPanelSize = 15;
handles.Min_size = 50;
handles.Max_size = 20000;
handles.Min_ecc = 0;
handles.Max_ecc = 1;
handles.PixelSize = 10;
handles.Scale_bar_size = 100;
handles.Sensitivity = 1;

warning('off','images:initSize:adjustingMag');    % Switches off the warning relatied to the display of large images
DefaultPath = 'E:\STORM data\2014_06_03 Radius analysis\';
handles.DefaultPath = DefaultPath;

PathName = DefaultPath;
handles.PathName = PathName;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VirusSegUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VirusSegUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in Load.
function Load_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% profile on
PixelSize = handles.PixelSize;

DefaultPath = handles.DefaultPath;
[FileName,PathName,FilterIndex] = uigetfile({'*.png'; '*.txt'},'Select a file (*.png or *.txt)', DefaultPath);

if FileName ~= 0
    % Clear the command window each time you load a new image
    clc
    DefaultPath = PathName;
    
    if FilterIndex == 1
        ColorImage = imread([PathName, FileName]);
        LocInfo = {};
        GreyImage = double(rgb2gray(ColorImage));
        Original_GreyImage = GreyImage;
    elseif FilterIndex == 2
        %[Original_GreyImage, LocInfo] = SRimage_LocFile( [PathName, FileName], PixelSize);
        LocInfo = Read_LocFile([PathName, FileName]);
        Original_GreyImage = Loc2SRimage(LocInfo, PixelSize);
        GreyImage = imadjust(Original_GreyImage,[0 handles.I_cutoff],[0 1],handles.Gamma);
        ColorImage = Grey2Color(GreyImage);
        
    end
    
    disp(['Loaded: ',PathName,FileName]);
    
    figure('Color','white','name',[PathName, FileName]);
    imshow(ColorImage)
    
    handles.PathName = PathName;
    handles.FileName = FileName;
    
    handles.Original_GreyImage = Original_GreyImage;
    handles.GreyImage = GreyImage;
    handles.ColorImage = ColorImage;
    handles.LocInfo = LocInfo;
    handles.DefaultPath = DefaultPath;
    
    guidata(hObject, handles);
end


function Min_size_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Min_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_size as text
%        str2double(get(hObject,'String')) returns contents of Min_size as a double
Min_size = str2double(get(hObject, 'String'));
if isnan(Min_size)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Min_size = Min_size;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Min_size_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Min_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_size_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Max_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_size as text
%        str2double(get(hObject,'String')) returns contents of Max_size as a double

Max_size = str2double(get(hObject, 'String'));
if isnan(Max_size)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Max_size = Max_size;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Max_size_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Max_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Min_ecc_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Min_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_ecc as text
%        str2double(get(hObject,'String')) returns contents of Min_ecc as a double
Min_ecc = str2double(get(hObject, 'String'));
if isnan(Min_ecc)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Min_ecc = Min_ecc;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Min_ecc_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Min_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Max_ecc_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Max_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_ecc as text
%        str2double(get(hObject,'String')) returns contents of Max_ecc as a double
Max_ecc = str2double(get(hObject, 'String'));
if isnan(Max_ecc)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Max_ecc = Max_ecc;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Max_ecc_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Max_ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Sort.
function Sort_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Sensitivity = handles.Sensitivity;
PathName = handles.PathName;
FileName = handles.FileName;
MinSize = handles.Min_size;
MaxSize = handles.Max_size;
EccMin = handles.Min_ecc;
EccMax = handles.Max_ecc;

Image = handles.GreyImage;
Image = OtsuThreshold(Image,Sensitivity);
LogicalIm = logical(Image);

figure('Color','white','name',['Logical image - ',PathName, FileName]);
imshow(LogicalIm);
title(['Sensitivity: ',num2str(Sensitivity)]);

ParamObject = regionprops(LogicalIm,'Area','Centroid','EquivDiameter','PixelList','Eccentricity');

BigObjects = find([ParamObject.Area] >= MinSize & [ParamObject.Area] <= MaxSize & [ParamObject.Eccentricity] < EccMax & [ParamObject.Eccentricity] >= EccMin);
SmallObjects = find([ParamObject.Area] < MinSize & [ParamObject.Eccentricity] < EccMax);

ImMask = DispAllObjects(LogicalIm, ParamObject, BigObjects);
ImPanel = ObjectPanel(handles.ColorImage, ParamObject, BigObjects, handles.PixelSize, handles.MaxPanelSize, handles.Scale_bar_size);

Centroids = cat(1, ParamObject(BigObjects).Centroid);

disp(' ');
disp('Sorting parameters:');
disp('------------------');
disp(['Sensitivity: ',num2str(Sensitivity)]);
disp(['Size: ',num2str(MinSize),' - ',num2str(MaxSize)])
disp(['Eccentricity: ', num2str(EccMin),' - ',num2str(EccMax)])
disp('------------------')
disp(['Total number of objects: ',num2str(size(Centroids,1))]);


handles.Centroids = Centroids;
handles.LogicalIm = LogicalIm;
handles.ImMask = ImMask;
handles.ParamObject = ParamObject;
handles.BigObjects = BigObjects;
handles.SmallObjects = SmallObjects;
handles.ImPanel = ImPanel;

guidata(hObject,handles)

% --- Executes on button press in Close_all.
function Close_all_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to Close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keep_guis;
disp('All figures closed.');


% --- Executes on button press in Review.
function Review_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Review (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ParamObject = handles.ParamObject;
BigObjects = handles.BigObjects;

HistSize(ParamObject(BigObjects));
HistParam(ParamObject(BigObjects));
DispInfo(ParamObject,BigObjects);


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


function PixelSize_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to PixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PixelSize as text
%        str2double(get(hObject,'String')) returns contents of PixelSize as a double
PixelSize = str2double(get(hObject, 'String'));
if isnan(PixelSize)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.PixelSize = PixelSize;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PixelSize_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to PixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_bar_size_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Scale_bar_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scale_bar_size as text
%        str2double(get(hObject,'String')) returns contents of Scale_bar_size as a double
Scale_bar_size = str2double(get(hObject, 'String'));
if isnan(Scale_bar_size)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Scale_bar_size = Scale_bar_size;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Scale_bar_size_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Scale_bar_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Show_mask.
function Show_mask_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Show_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ParamObject = handles.ParamObject;
BigObjects = handles.BigObjects;
ColorImage = handles.ColorImage;
Centroids = handles.Centroids;
LogicalIm = handles.LogicalIm;

figure('Color','white','name','Color image');
imshow(ColorImage)
title(['Total number of Objects: ',num2str(max(size(BigObjects)))])


if size(BigObjects,2) == size(ParamObject,2)
    ImMask = LogicalIm;
else
    ImMask = DispAllObjects(LogicalIm, ParamObject, BigObjects);
end

ObjectPanel(ColorImage, ParamObject, BigObjects, handles.PixelSize, handles.MaxPanelSize, handles.Scale_bar_size);
set(gcf,'name','Particle panel');

ColorImage = ShowMask_on_Image(ImMask, ColorImage);
ColorImage = AddCentroid2Image(Centroids, ColorImage);

figure('Color','white','name','Color image - with mask');
imshow(ColorImage)
title(['Total number of Objects: ',num2str(max(size(BigObjects)))])

% ObjectPanel(256*ImMask, ParamObject, BigObjects, handles.PixelSize, handles.MaxPanelSize, Scale_bar_size);
ImPanelMask = ObjectPanel(ColorImage, ParamObject, BigObjects, handles.PixelSize, handles.MaxPanelSize, handles.Scale_bar_size);
set(gcf,'name','Particle panel - with mask');

handles.ImPanelMask = ImPanelMask;
guidata(hObject,handles)




% --- Executes on button press in Save.
function Save_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is Save particle panel

PathName = handles.PathName;
PixelSize = handles.PixelSize;
Centroids = handles.Centroids;

% Density = max(size(handles.BigObjects))/(handles.PixelSize^2*size(handles.ImMask,1)*size(handles.ImMask,2));
% disp(['Density of particles: ', num2str(Density*1000000),' particles / um^2']);
% DensityBg = max(size(handles.SmallObjects ))/(handles.PixelSize^2*size(handles.ImMask,1)*size(handles.ImMask,2));
% disp(['Density of background objects: ', num2str(DensityBg*1000000),' objects / um^2 (e.g. ', num2str(round(DensityBg/Density)),'/1)']);

Filename = strsplit(handles.FileName,'.');
Panel_Filename = [handles.PathName,'Particles_',Filename{1},'.png'];
imwrite(handles.ImPanel,Panel_Filename,'png');

disp('----------------');
disp('Files saved as:');
disp(['Saved: ',Panel_Filename]);

PanelMask_Filename = [handles.PathName,'Particles_',Filename{1},'_Mask.png'];
imwrite(handles.ImPanelMask,PanelMask_Filename,'png');
disp(['Saved: ',PanelMask_Filename]);

% Converting the centroids into nm
Centroids = (Centroids - 1)*PixelSize;
Filename = strsplit(handles.FileName,'.');
Centroids_filename = strcat(PathName,'Centroids_',Filename{1},'.txt');
Header = {'X','Y'};

dlmwrite(Centroids_filename,Header,'delimiter','\t','newline','pc');
dlmwrite(Centroids_filename,Centroids,'-append','delimiter','\t','newline','pc','precision','%6.1f');
disp(['Saved: ',Centroids_filename]);



function MaxPanelSize_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to MaxPanelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxPanelSize as text
%        str2double(get(hObject,'String')) returns contents of MaxPanelSize as a double
MaxPanelSize = str2double(get(hObject, 'String'));
if isnan(MaxPanelSize)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.MaxPanelSize = MaxPanelSize;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function MaxPanelSize_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to MaxPanelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Select_particles.
function Select_particles_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Select_particles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

BigObjects = handles.BigObjects;
ParamObject = handles.ParamObject;


TableData = Delete_GUI(size(BigObjects,2));
BigObjects(TableData) = [];

handles.BigObjects = BigObjects ;

Centroids = cat(1, ParamObject(BigObjects).Centroid);
handles.Centroids = Centroids;

ImPanel = ObjectPanel(handles.ColorImage, handles.ParamObject, BigObjects, handles.PixelSize, handles.MaxPanelSize, handles.Scale_bar_size);
handles.ImPanel = ImPanel ;

guidata(hObject,handles)


% --- Executes on button press in Review_loc_file.
function Review_loc_file_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Review_loc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LocInfo = handles.LocInfo;

if isempty(LocInfo) == 1
    disp('No localization file loaded.')
    PathName = handles.PathName;
    [FileName_Loc,PathName_Loc,~] = uigetfile( '*.txt','Select a file (*.txt)',PathName);
    disp(' ');
    disp('Read localization file:');
    disp([PathName_Loc, FileName_Loc]);
    LocInfo = Read_LocFile( [PathName_Loc,FileName_Loc] );
    handles.PathName_Loc = PathName_Loc;
    handles.FileName_Loc = FileName_Loc;
else
    PathName_Loc = handles.PathName_Loc;
    FileName_Loc = handles.FileName_Loc;
end

hFigure_LocRev = Review_loc_file(LocInfo);
set(hFigure_LocRev,'name',[PathName_Loc, FileName_Loc]);

Filename = strsplit(FileName_Loc,'.');
LocReview_filename = strcat(PathName_Loc,'LocReview_',Filename{1},'.png');
print(hFigure_LocRev, '-dpng', LocReview_filename);
disp(['Saved: ',LocReview_filename])

handles.LocInfo = LocInfo;
guidata(hObject,handles)


% --- Executes on button press in Centre_mass.
function Centre_mass_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Centre_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LocInfo = handles.LocInfo;

if isempty(LocInfo) == 1
    disp('No localization file loaded.')
    PathName = handles.PathName;
    [FileName_Loc,PathName_Loc,~] = uigetfile( '*.txt','Select a file (*.txt)',PathName);
    disp(' ');
    disp('Read localization file:');
    disp([PathName_Loc, FileName_Loc]);
    LocInfo = Read_LocFile( [PathName_Loc,FileName_Loc] );
    handles.PathName_Loc = PathName_Loc;
    handles.FileName_Loc = FileName_Loc;
end

Centroids = handles.Centroids;
ADC_cutoff = handles.ADC_cutoff;
PixelSize = handles.PixelSize;
Frame_min = handles.Frame_min;
Frame_max = handles.Frame_max;

Frame = LocInfo(:,3);
ADC = LocInfo(:,4);
LocInfo(Frame < Frame_min | Frame > Frame_max | ADC < ADC_cutoff, :) = [];

ParticleInfo = Centre_mass(LocInfo, Centroids, PixelSize, handles.Radius_cutoff);
hFig_CMass = COM_display(ParticleInfo, handles.Scale_bar_size, handles.Gamma, handles.I_cutoff);

handles.hFig_CMass = hFig_CMass;
handles.ParticleInfo = ParticleInfo;

disp(' ');
disp('Centre of mass analysis parameters:');
disp(['ADC cut-off: ',num2str(ADC_cutoff)]);
disp(['Frame:',num2str(Frame_min),' - ',num2str(Frame_max)]);

guidata(hObject,handles)



function ADC_cutoff_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to ADC_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ADC_cutoff as text
%        str2double(get(hObject,'String')) returns contents of ADC_cutoff as a double
ADC_cutoff = str2double(get(hObject, 'String'));
if isnan(ADC_cutoff)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.ADC_cutoff = ADC_cutoff;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ADC_cutoff_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to ADC_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Radius_cutoff_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Radius_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Radius_cutoff as text
%        str2double(get(hObject,'String')) returns contents of Radius_cutoff as a double
Radius_cutoff = str2double(get(hObject, 'String'));
if isnan(Radius_cutoff)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Radius_cutoff = Radius_cutoff;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Radius_cutoff_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Radius_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frame_min_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Frame_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frame_min as text
%        str2double(get(hObject,'String')) returns contents of Frame_min as a double
Frame_min = str2double(get(hObject, 'String'));
if isnan(Frame_min)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Frame_min = Frame_min;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Frame_min_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Frame_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frame_max_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Frame_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frame_max as text
%        str2double(get(hObject,'String')) returns contents of Frame_max as a double
Frame_max = str2double(get(hObject, 'String'));
if isnan(Frame_max)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.Frame_max = Frame_max;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Frame_max_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Frame_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_centre_mass.
function Save_centre_mass_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Save_centre_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hFig_CMass = handles.hFig_CMass;
PathName = handles.PathName;
ParticleInfo = handles.ParticleInfo;

Filename = strsplit(handles.FileName,'.');
CM_filename = strcat(PathName,'COM_',Filename{1},'.png');
print(hFig_CMass, '-dpng', CM_filename)
disp(['Saved: ',CM_filename])

PI_filename = strcat(PathName,'ParticleInfo_',Filename{1},'.txt');
Header = {'X','Y','f','A','n'};

dlmwrite(PI_filename,Header,'delimiter','\t','newline','pc');
dlmwrite(PI_filename,ParticleInfo,'-append','delimiter','\t','newline','pc');
disp(['Saved: ',PI_filename]);



% --- Executes on button press in Load_mask.
function Load_mask_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PathName = handles.PathName;
ColorImage = handles.ColorImage;
GreyImage = handles.GreyImage;
PixelSize = handles.PixelSize;

Mask = Create_Color_mask({PathName, ColorImage, GreyImage, PixelSize});

ColorImage = ColorImage.*uint8(cat(3,Mask,Mask,Mask));
GreyImage = GreyImage.*double(Mask);
disp('Done.');
handles.ColorImage = ColorImage;
handles.GreyImage = GreyImage;

guidata(hObject,handles)


% --- Executes on button press in Multiple_centre_mass.
function Multiple_centre_mass_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Multiple_centre_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PathName = handles.PathName;
[ParticleInfo, FolderName] = Multiple_COM(PathName);
% assignin('base','ParticleInfo',ParticleInfo);

Radius = sqrt(ParticleInfo(:,1).^2 + ParticleInfo(:,2).^2);
ADC_cutoff = handles.ADC_cutoff;
Frame_min = handles.Frame_min;
Frame_max = handles.Frame_max;

Frame = ParticleInfo(:,3);
ADC = ParticleInfo(:,4);
ParticleInfo(Frame < Frame_min | Frame > Frame_max | ADC < ADC_cutoff | Radius > handles.Radius_cutoff, :) = [];
% assignin('base','ParticleInfo2',ParticleInfo);

hFig_CMass = COM_display(ParticleInfo, handles.Scale_bar_size, handles.Gamma, handles.I_cutoff);

button = questdlg('Save Multiple COM analysis?','Save Multiple COM','Yes','No','No');
if strcmp(button,'Yes')
    CM_filename = strcat(FolderName,'\COM_multiple.png');
    print(hFig_CMass, '-dpng', CM_filename)
    disp(['Saved: ',CM_filename])
end


% --- Executes on button press in Adjust_image.
function Adjust_image_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Adjust_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Original_GreyImage = handles.Original_GreyImage;

GreyImage = imadjust(Original_GreyImage,[0 handles.I_cutoff],[0 1],handles.Gamma);
ColorImage = Grey2Color(GreyImage);

handles.GreyImage = GreyImage;
handles.ColorImage = ColorImage;

h2 = figure('Color','white','name','Color image - adjusted');
imshow(ColorImage)
set(h2,'name',[handles.PathName, handles.FileName])

guidata(hObject,handles)



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



function I_cutoff_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to I_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I_cutoff as text
%        str2double(get(hObject,'String')) returns contents of I_cutoff as a double
I_cutoff = str2double(get(hObject, 'String'));
if isnan(I_cutoff)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.I_cutoff = I_cutoff;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function I_cutoff_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to I_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load_loc_file.
function Load_loc_file_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Load_loc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PathName = handles.PathName;
[FileName,PathName,~] = uigetfile( '*.txt','Select a file (*.txt)',PathName);

if FileName ~= 0
    disp(' ');
    disp('Read localization file:');
    disp([PathName, FileName]);
    LocInfo = Read_LocFile( [PathName,FileName] );
    
    handles.PathName_Loc = PathName;
    handles.FileName_Loc = FileName;
    
    handles.LocInfo = LocInfo;
    guidata(hObject,handles)
end


% --- Executes on button press in Load_centroid.
function Load_centroid_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Load_centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PathName = handles.PathName;
PixelSize = handles.PixelSize;
% Scale_bar_size = handles.Scale_bar_size;
% ColorImage = handles.ColorImage;

[FileName_Centroid,PathName_Centroid,~] = uigetfile( '*.txt','Select a file (*.txt)',PathName);

if FileName_Centroid ~= 0
    disp(' ');
    disp('Read localization file containing centroids:');
    disp([PathName_Centroid, FileName_Centroid]);
    Centroids = dlmread([PathName_Centroid,FileName_Centroid],'\t',1,0);
    disp(['Number of centroids loaded: ',num2str(size(Centroids,1))]);
    
% % This line should be used when loading centroids from localization file
% % obtained from mTurquoise localization
%     CentroidInfo = Read_LocFile( [PathName_Centroid,FileName_Centroid] );
%     Centroids = CentroidInfo(:,1:2);
%     
%     % Linear shifts
%     xy_shift = Centroid_shift_GUI();
%     x_shift = xy_shift(1);
%     y_shift = xy_shift(2);
%     
%     Centroids(:,1) = Centroids(:,1) + x_shift;
%     Centroids(:,2) = Centroids(:,2) + y_shift;
%     disp('Linear shift applied to centroids:');
%     disp(['x: ',num2str(x_shift),' nm - y: ',num2str(y_shift),' nm']);
%     

% Converting the centroids in nm into pixel units
    Centroids = Centroids / (PixelSize) + 1;
%     
%     DispListObject_Centroids(ColorImage, Centroids, PixelSize, Scale_bar_size);
%     TableData = Delete_GUI(size(Centroids,1));
%     Centroids(TableData,:) = [];
%     DispListObject_Centroids(ColorImage, Centroids, PixelSize, Scale_bar_size);
%     TableData = Delete_GUI(size(Centroids,1));
%     Centroids(TableData,:) = [];
%     DispListObject_Centroids(ColorImage, Centroids, PixelSize, Scale_bar_size);
    
    handles.Centroids = Centroids;
    guidata(hObject,handles)
end


% --- Executes on button press in Hough_tr.
function Hough_tr_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Hough_tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PixelSize = handles.PixelSize;
ColorImage = handles.ColorImage;
MaxPanelSize = handles.MaxPanelSize;
Scale_bar_size = handles.Scale_bar_size;

ColorImage_raw = ColorImage;
disp('------------------');
disp('Hough transform:');


% Apply some Gaussian smoothing
Smoothing = 1;
if Smoothing == 1
    disp('Gaussian smoothing applied.');
    GaussFilter = fspecial('gaussian',[5 5],0.75);
    ColorImage = imfilter(ColorImage,GaussFilter,'same');
    figure('Color','white','name','Hough transform - Smoothing');
    subplot(1,2,1)
    imshow(ColorImage_raw,[]);
    title 'Raw image'
    subplot(1,2,2)
    imshow(ColorImage,[]);
    title 'With Gaussian filter'
elseif Smoothing == 0
    disp('No smoothing applied.');
    figure('Color','white','name','Hough transform - Smoothing');
    imshow(ColorImage,[]);
    title 'No smoothing applied'
end

% Create a dialog window to get the parameters from user
prompt = {'Radius min (nm):','Radius max (nm):','Sensitivity:','Edge threshold:'};
dlg_title = 'Hough transform parameters';
num_lines = 1;
def = {'70','200','0.95','0.05'};
Hough_params = inputdlg(prompt,dlg_title,num_lines,def);
% assignin('base','Answer',answer);

% Good initial values are as following:
% Radius_range = [70 150]; % nm
% Sensitivity = 0.95;
% EdgeThresh = 0.05;

if isempty(Hough_params) == 0
    Radius_range = [str2double(Hough_params{1}), str2double(Hough_params{2})];
    Sensitivity = str2double(Hough_params{3});
    EdgeThresh = str2double(Hough_params{4});
        
    % Convert in pixels
    allowed_radius = round(Radius_range/PixelSize); % in pixels
    
    %% Perform the Hough transform
    
    Method = 'PhaseCode';
    % Method = 'TwoStage';
    
    disp(['Allowed radius range: ',num2str(allowed_radius), ' (in pixels)']);
    disp(['Sensitivity: ',num2str(Sensitivity),' - Edge threshold: ', num2str(EdgeThresh)]);
    disp('Computing Hough transform...');
    tic
    [centers, radii, circle_strength] = imfindcircles(ColorImage,allowed_radius,'ObjectPolarity', 'bright', 'Sensitivity',Sensitivity,'Method',Method,'EdgeThreshold',EdgeThresh);
    toc
    
    % Total number of circles found by the Hough transform
    n_circles = size(radii,1);
    disp(['Number of circles found: ',num2str(n_circles)]);
    % assignin('base', 'circle_strength', circle_strength);
    
    % Display a graph to adjust the sensitivity (=1-circle_strength)
    figure('Color','white','name','Hough transform - Sensitivity');
    plot(1-circle_strength,0:(n_circles-1),'r+')
    ylabel 'Number of circles found'
    xlabel 'Sensitivity'
    grid on
    
    % figure('Color','white');
    % hist(1-circle_strength,50)
    % xlabel 'Sensitivity'
    % ylabel 'Occurences'
    
    % Create the image with the detected parameters on
%     Image_Hough = ColorImage;
    Image_Hough_raw = ColorImage_raw;
    
    % This one is time-consumming --> god to improve if possible
    tic
    [ ImCircle, ImCentres, LogicalIm ] = CreateCircleImage(size(ColorImage), centers, radii);
    toc
    
    % Display mask image (detected circles)
    figure('Color','white','name','Hough transform - Binary mask image');
    imshow(LogicalIm,[]);
    
    % Calculate RGB images for display purposes
    Image_Hough_raw(:,:,3) = Image_Hough_raw(:,:,3) + uint8(255*ImCircle);
%     Image_Hough(:,:,3) = Image_Hough(:,:,3) + uint8(255*ImCircle);
    
    Image_Hough_raw(:,:,2) = Image_Hough_raw(:,:,2) + uint8(255*ImCentres);
%     Image_Hough(:,:,2) = Image_Hough(:,:,2) + uint8(255*ImCentres);
    
%     figure('Color','white');
%     subplot(2,2,1)
%     imshow(ColorImage_raw,[]);
%     title 'Raw image'
%     subplot(2,2,2)
%     imshow(ColorImage,[]);
%     title 'With Gaussian filter'
%     subplot(2,2,3)
%     imshow(Image_Hough_raw);
%     title 'Hough circle detection (raw)'
%     subplot(2,2,4)
%     imshow(Image_Hough);
%     title 'Hough circle detection (smoothed)'
    
    
    %% Create the ParamObject variables as identical to that obtained from object detection function
    
    if n_circles == 0
        disp('No circles were found by Hough transform.')
    else
        ParamObject(n_circles).Centroid = [0 0];
        ParamObject(n_circles).EquivDiameter = 0;
        BigObjects = 1:n_circles;
        for i = 1:n_circles
            ParamObject(i).Centroid = centers(i,:);
            ParamObject(i).EquivDiameter = 2*radii(i);
        end
        
        ObjectPanel(Image_Hough_raw, ParamObject, BigObjects, PixelSize, MaxPanelSize, Scale_bar_size);
        set(gcf,'name','Particle panel - with mask')
        ImPanel = ObjectPanel(ColorImage_raw, ParamObject, BigObjects, PixelSize, MaxPanelSize, Scale_bar_size);
        Centroids = cat(1, ParamObject(BigObjects).Centroid);
        
        handles.ParamObject = ParamObject;
        handles.ImPanel = ImPanel;
        handles.ParamObject = ParamObject;
        handles.BigObjects = BigObjects;
        handles.Centroids = Centroids ;
        handles.LogicalIm = LogicalIm;
        
        guidata(hObject,handles)
    end
    
    % End of the if isempty(Hough_params) == 0
end







