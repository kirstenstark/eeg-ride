function varargout = RIDE_GUI(varargin)
% RIDE_GUI MATLAB code for RIDE_GUI.fig
%      RIDE_GUI, by itself, creates a new RIDE_GUI or raises the existing
%      singleton*.
%
%      H = RIDE_GUI returns the handle to a new RIDE_GUI or the handle to
%      the existing singleton*.
%
%      RIDE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIDE_GUI.M with the given input arguments.
%
%      RIDE_GUI('Property','Value',...) creates a new RIDE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RIDE_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RIDE_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RIDE_GUI

% Last Modified by GUIDE v2.5 09-Jan-2015 12:27:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RIDE_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RIDE_GUI_OutputFcn, ...
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


% --- Executes just before RIDE_GUI is made visible.
function RIDE_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RIDE_GUI (see VARARGIN)

% Choose default command line output for RIDE_GUI
handles.output = hObject;
if ~isempty(varargin) handles.results = varargin{1}; end
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes RIDE_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RIDE_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_ride.
function run_ride_Callback(hObject, eventdata, handles)
% hObject    handle to run_ride (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp(size(handles.data));
cfg = RIDE_cfg(handles.cfg);

temp = find(strcmpi(cfg.comp.latency,'rt')); if ~isempty(temp) cfg.comp.latency{temp} = handles.rt;end

assignin('base','data',handles.data);
assignin('base','cfg',cfg);
handles.results = RIDE_call(handles.data,cfg);disp(isfield(handles,'results'));
RIDE_plot(handles.results,{'erp','s','c','r'},handles.chan1,handles.axes1);
axis(handles.axes1,'tight');
xlabel(handles.axes1,'time after stimulus (ms)');
ylabel(handles.axes1,'potential');
guidata(hObject, handles);

assignin('base','results',handles.results);


function chan_Callback(hObject, eventdata, handles)
% hObject    handle to chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chan as text
%        str2double(get(hObject,'String')) returns contents of chan as a double
handles.chan1 = str2double(get(hObject,'String'));
guidata(hObject, handles);
if isfield(handles,'results')
    if handles.plot_item == 1
        RIDE_plot(handles.results,{'erp',handles.cfg.comp.name{:}},handles.chan1,handles.axes1);
    end
    if handles.plot_item == 2
        RIDE_plot(handles.results,{'erp','erp_new'},handles.chan1,handles.axes1);
    end
    axis(handles.axes1,'tight');
    xlabel(handles.axes1,'time after stimulus (ms)');
    ylabel(handles.axes1,'potential');
end


% --- Executes during object creation, after setting all properties.
function chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.chan1 = str2double(get(hObject,'String'));%disp(str2double(get(hObject,'String')));
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over chan.
function chan_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function samp_interval_Callback(hObject, eventdata, handles)
% hObject    handle to samp_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samp_interval as text
%        str2double(get(hObject,'String')) returns contents of samp_interval as a double
handles.cfg.samp_interval = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function samp_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samp_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epoch_twd_Callback(hObject, eventdata, handles)
% hObject    handle to epoch_twd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch_twd as text
%        str2double(get(hObject,'String')) returns contents of epoch_twd as a double
temp = textscan(get(hObject,'String'),'%f'); 
handles.cfg.epoch_twd = temp{1}';
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function epoch_twd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch_twd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name1_Callback(hObject, eventdata, handles)
% hObject    handle to name1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name1 as text
%        str2double(get(hObject,'String')) returns contents of name1 as a double



% --- Executes during object creation, after setting all properties.
function name1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.cfg.comp.name{1} = get(hObject,'String');
guidata(hObject, handles);


% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('load');%disp(size(data));
handles.data = data;if exist('rt','var') handles.rt = rt;end
guidata(hObject,handles);disp('load data success');


% --- Executes on button press in load_res.
function load_res_Callback(hObject, eventdata, handles)
% hObject    handle to load_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('load');handles.results = results;guidata(hObject,handles);


function high_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to high_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_cutoff as text
%        str2double(get(hObject,'String')) returns contents of high_cutoff as a double



% --- Executes during object creation, after setting all properties.
function high_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.cfg.high_cutoff = str2double(get(hObject,'String'));
guidata(hObject, handles);


function re_samp_Callback(hObject, eventdata, handles)
% hObject    handle to re_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of re_samp as text
%        str2double(get(hObject,'String')) returns contents of re_samp as a double
handles.cfg.re_samp = str2double(get(hObject,'String'));%disp(handles.cfg.re_samp);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function re_samp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to re_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latency_search_Callback(hObject, eventdata, handles)
% hObject    handle to latency_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency_search as text
%        str2double(get(hObject,'String')) returns contents of latency_search as a double



% --- Executes during object creation, after setting all properties.
function latency_search_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.cfg.latency_search = get(hObject,'String');
guidata(hObject, handles);


function dur1_Callback(hObject, eventdata, handles)
% hObject    handle to dur1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dur1 as text
%        str2double(get(hObject,'String')) returns contents of dur1 as a double
handles.cfg.dur{1} = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function dur1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dur2_Callback(hObject, eventdata, handles)
% hObject    handle to dur2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dur2 as text
%        str2double(get(hObject,'String')) returns contents of dur2 as a double
handles.cfg.dur{2} = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function dur2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dur3_Callback(hObject, eventdata, handles)
% hObject    handle to dur3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dur3 as text
%        str2double(get(hObject,'String')) returns contents of dur3 as a double
handles.cfg.dur{3} = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function dur3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dur4_Callback(hObject, eventdata, handles)
% hObject    handle to dur4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dur4 as text
%        str2double(get(hObject,'String')) returns contents of dur4 as a double
handles.cfg.dur{4} = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function dur4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function twd1_Callback(hObject, eventdata, handles)
% hObject    handle to twd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of twd1 as text
%        str2double(get(hObject,'String')) returns contents of twd1 as a double
temp = textscan(get(hObject,'String'),'%f');
handles.cfg.comp.twd{1} = temp{1}';
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function twd1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to twd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latency1_Callback(hObject, eventdata, handles)
% hObject    handle to latency1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency1 as text
%        str2double(get(hObject,'String')) returns contents of latency1 as a double

if strcmpi(get(hObject,'String'),'0') handles.cfg.comp.latency{1} = 0;end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function latency1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if strcmpi(get(hObject,'String'),'0') handles.cfg.comp.latency{1} = 0;end
guidata(hObject, handles);




function name2_Callback(hObject, eventdata, handles)
% hObject    handle to name2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name2 as text
%        str2double(get(hObject,'String')) returns contents of name2 as a double
handles.cfg.comp.name{2} = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function name2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function twd2_Callback(hObject, eventdata, handles)
% hObject    handle to twd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of twd2 as text
%        str2double(get(hObject,'String')) returns contents of twd2 as a double
temp = textscan(get(hObject,'String'),'%f');
handles.cfg.comp.twd{2} = temp{1}'; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function twd2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to twd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latency2_Callback(hObject, eventdata, handles)
% hObject    handle to latency2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency2 as text
%        str2double(get(hObject,'String')) returns contents of latency2 as a double
temp = get(hObject,'String'); 
if strcmpi(temp,'rt') 
    handles.cfg.comp.latency{2} = temp;
    if isfield(handles,'rt') assignin('base','rt',handles.rt);end
else handles.cfg.comp.latency{2} = 'unknown';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function latency2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name3_Callback(hObject, eventdata, handles)
% hObject    handle to name3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name3 as text
%        str2double(get(hObject,'String')) returns contents of name3 as a double
handles.cfg.comp.name{3} = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function name3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function twd3_Callback(hObject, eventdata, handles)
% hObject    handle to twd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of twd3 as text
%        str2double(get(hObject,'String')) returns contents of twd3 as a double
temp = textscan(get(hObject,'String'),'%f');
handles.cfg.comp.twd{3} = temp{1}'; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function twd3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to twd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latency3_Callback(hObject, eventdata, handles)
% hObject    handle to latency3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency3 as text
%        str2double(get(hObject,'String')) returns contents of latency3 as a double
temp = get(hObject,'String'); 
if strcmpi(temp,'rt') 
    handles.cfg.comp.latency{3} = temp;
    if isfield(handles,'rt') assignin('base','rt',handles.rt);end
else handles.cfg.comp.latency{3} = 'unknown';
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function latency3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name4_Callback(hObject, eventdata, handles)
% hObject    handle to name4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name4 as text
%        str2double(get(hObject,'String')) returns contents of name4 as a double
handles.cfg.comp.name{4} = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function name4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function twd4_Callback(hObject, eventdata, handles)
% hObject    handle to twd4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of twd4 as text
%        str2double(get(hObject,'String')) returns contents of twd4 as a double
temp = textscan(get(hObject,'String'),'%f');
handles.cfg.comp.twd{4} = temp{1}'; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function twd4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to twd4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latency4_Callback(hObject, eventdata, handles)
% hObject    handle to latency4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency4 as text
%        str2double(get(hObject,'String')) returns contents of latency4 as a double
temp = get(hObject,'String'); 
if strcmpi(temp,'rt') 
    handles.cfg.comp.latency{4} = temp;
    if isfield(handles,'rt') assignin('base','rt',handles.rt);end
else handles.cfg.comp.latency{4} = 'unknown';
end

% --- Executes during object creation, after setting all properties.
function latency4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pop_out.
function pop_out_Callback(hObject, eventdata, handles)
% hObject    handle to pop_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = figure;temp2=copyobj(handles.axes1,temp);set(temp2,'Units','normalized');set(temp2,'position',[0.1,0.1,0.8,0.8]);


% --- Executes on button press in save_res.
function save_res_Callback(hObject, eventdata, handles)
% hObject    handle to save_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'results')
    results = handles.results;uisave('results');
end


% --- Executes on button press in batch_proc.
function batch_proc_Callback(hObject, eventdata, handles)
% hObject    handle to batch_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'cfg')
    [a,b] = uigetfile('*.mat','MultiSelect','on');rt_i = find(strcmpi(handles.cfg.comp.latency,'rt'));
    for j = 1:length(a)
        load([b,a{j}]); if ~isempty(rt_i) handles.cfg.comp.latency{rt_i} = rt;end
        disp(handles.cfg.comp.latency{2});
        cfg = RIDE_cfg(handles.cfg);
        results = RIDE_call(data,cfg);
        save([b,'results_',a{j}],'results');
    end
end
    


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
handles.plot_item = get(hObject,'Value');guidata(hObject,handles);

    if handles.plot_item == 1
        RIDE_plot(handles.results,{'erp',handles.cfg.comp.name{:}},handles.chan1,handles.axes1);
    end
    if handles.plot_item == 2
        RIDE_plot(handles.results,{'erp','erp_new'},handles.chan1,handles.axes1);
    end
    axis(handles.axes1,'tight');
    xlabel(handles.axes1,'time after stimulus (ms)');
    ylabel(handles.axes1,'potential');
    
    
% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.plot_item = get(hObject,'Value');guidata(hObject,handles);
