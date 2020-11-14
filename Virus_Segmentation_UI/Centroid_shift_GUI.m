function varargout = Centroid_shift_GUI(varargin)
% CENTROID_SHIFT_GUI MATLAB code for Centroid_shift_GUI.fig
%      CENTROID_SHIFT_GUI, by itself, creates a new CENTROID_SHIFT_GUI or raises the existing
%      singleton*.
%
%      H = CENTROID_SHIFT_GUI returns the handle to a new CENTROID_SHIFT_GUI or the handle to
%      the existing singleton*.
%
%      CENTROID_SHIFT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CENTROID_SHIFT_GUI.M with the given input arguments.
%
%      CENTROID_SHIFT_GUI('Property','Value',...) creates a new CENTROID_SHIFT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Centroid_shift_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Centroid_shift_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Centroid_shift_GUI

% Last Modified by GUIDE v2.5 10-Jan-2014 11:53:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Centroid_shift_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Centroid_shift_GUI_OutputFcn, ...
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


% --- Executes just before Centroid_shift_GUI is made visible.
function Centroid_shift_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Centroid_shift_GUI (see VARARGIN)

% Choose default command line output for Centroid_shift_GUI
handles.output = hObject;

set(handles.x_shift, 'String', '0');
set(handles.y_shift, 'String', '0');
handles.x_shift = 0;
handles.y_shift = 0;

% Update handles structure
guidata(hObject, handles);
uiwait(hObject);

% UIWAIT makes Centroid_shift_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Centroid_shift_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
xy_shift = [handles.x_shift, handles.y_shift];
varargout{1} = xy_shift;

close(hObject);


function x_shift_Callback(hObject, eventdata, handles)
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
function x_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_shift_Callback(hObject, eventdata, handles)
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
function y_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume();

