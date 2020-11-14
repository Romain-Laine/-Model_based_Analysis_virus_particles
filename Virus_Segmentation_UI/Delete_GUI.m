function varargout = Delete_GUI(varargin)
% DELETE_GUI MATLAB code for Delete_GUI.fig
%      DELETE_GUI, by itself, creates a new DELETE_GUI or raises the existing
%      singleton*.
%
%      H = DELETE_GUI returns the handle to a new DELETE_GUI or the handle to
%      the existing singleton*.
%
%      DELETE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DELETE_GUI.M with the given input arguments.
%
%      DELETE_GUI('Property','Value',...) creates a new DELETE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Delete_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Delete_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Delete_GUI

% Last Modified by GUIDE v2.5 29-Oct-2013 12:27:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Delete_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Delete_GUI_OutputFcn, ...
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


% --- Executes just before Delete_GUI is made visible.
function Delete_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Delete_GUI (see VARARGIN)

% Choose default command line output for Delete_GUI
handles.output = hObject;


disp('Remove objects from the panel:')
n_object = varargin{1};

n = ceil(sqrt(n_object));
handles.n = n;

Pixel_width = 30;
Pixel_height = 20;
ColWidth = Pixel_width*ones(1,n);
ColWidth = num2cell(ColWidth);

set(handles.UI_table,'ColumnWidth',ColWidth);

rect = [0 0 Pixel_width*(n+1)+50 Pixel_height*n + 100];
set(handles.figure1,'Units','pixels',...
    'Position',rect);

Data_array = false(n,n);
Data_cell = num2cell(Data_array);
columnformat = cell(1,n);
columnformat(:) = {'logical'};
columneditable =  true(1,n);

rect = [0 0 1 1];
set(handles.UI_table,'Units','normalized',...
    'Position',rect,...
    'Data',Data_cell,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', columneditable);

set(handles.Ok_button,'Units','pixels',...
    'Position',[30 30 100 50]);

TableData = zeros(n,n);
handles.TableData = TableData;

handles.n_object = n_object;

% Update handles structure
guidata(hObject, handles);
uiwait();

% UIWAIT makes Delete_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Delete_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargout = handles.Data;

TableData = handles.TableData;
n_object = handles.n_object;
n = handles.n;

TableData = reshape(TableData',1,n^2);
if max(size(TableData)) > n_object
    TableData((n_object+1):end) = [];
end

TableData = logical(TableData);
varargout = {TableData};

close(hObject);


% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in Ok_button.
function Ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to Ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TableData = get(handles.UI_table, 'Data');
handles.TableData = cell2mat(TableData);
disp('Ok');
disp(handles.TableData);

guidata(hObject, handles);
uiresume();



% --- Executes during object creation, after setting all properties.
function UI_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UI_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in UI_table.
function UI_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to UI_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
