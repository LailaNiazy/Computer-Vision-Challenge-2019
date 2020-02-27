function varargout = challengeGUI(varargin)
% GUITEST MATLAB code for guiTest.fig
%      GUITEST, by itself, creates a new GUITEST or raises the existing
%      singleton*.
%
%      H = GUITEST returns the handle to a new GUITEST or the handle to
%      the existing singleton*.
%
%      GUITEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITEST.M with the given input arguments.
%
%      GUITEST('Property','Value',...) creates a new GUITEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiTest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiTest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiTest

% Last Modified by GUIDE v2.5 25-Aug-2018 12:26:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiTest_OpeningFcn, ...
                   'gui_OutputFcn',  @guiTest_OutputFcn, ...
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


% --- Executes just before guiTest is made visible.
function guiTest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiTest (see VARARGIN)

% Choose default command line output for guiTest
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiTest wait for user response (see UIRESUME)
% uiwait(handles.figure1);
addpath(genpath('RectifKitE'));
addpath(genpath('Functions'));
addpath(genpath('DisparityMap'));
addpath(genpath('Rendering_view_synthesis'));
I = imread('img/PlatzhalterOben.png');
J = imread('img/PlatzhalterOben.png');
K = imread('img/PlatzhalterUnten.jpg');
axes(handles.LeftPicture);
imshow(I);
axes(handles.RightPicture);
imshow(J);
axes(handles.result);
imshow(K);
handles.ancle = 0.5;
handles.kalib = zeros(3,3);
handles.depthMap = 0;
handles.LeftPic = I;
handles.RightPic = J;
guidata(hObject,handles)


% --- Outputs from this function are returned to the command line.
function varargout = guiTest_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
output = free_viewpoint(handles.LeftPic,handles.RightPic,handles.ancle, handles.kalib, handles.depthMap);
axes(handles.result);
elapsed_time = toc;
imshow(output);
textLabel = sprintf('Task has perfomed in %f seconds', elapsed_time);
set(handles.performance, 'String', textLabel);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%textLabel = sprintf('Variable C = %f', C);
label = sprintf('p = %g', get(hObject,'Value'));
set(handles.winkel, 'String', label);
handles.ancle = get(hObject,'Value');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(hObject,'Value')
    I = imread('img/L1.JPG');
    axes(handles.LeftPicture);
    imshow(I);
    handles.LeftPic = I;
    load('Kalibrierungsmatrix.mat');
    handles.depthMap = 1;
    handles.kalib = K1_opt;
    guidata(hObject,handles)
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(hObject,'Value')
    I = imread('img/R1.JPG');
    axes(handles.RightPicture);
    imshow(I);
    handles.RightPic = I;
    guidata(hObject,handles)
end

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
if get(hObject,'Value')
    I = imread('img/L2.JPG');
    axes(handles.LeftPicture);
    imshow(I);
    handles.LeftPic = I;
    handles.depthMap = 2;
    load('Kalibrierungsmatrix.mat');
    handles.kalib = K2_opt;
    guidata(hObject,handles)
end

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4

if get(hObject,'Value')
    I = imread('img/R2.JPG');
    axes(handles.RightPicture);
    imshow(I);
    handles.RightPic = I;
    guidata(hObject,handles)
end
