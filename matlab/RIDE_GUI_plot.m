function varargout = RIDE_GUI_plot(varargin)
% RIDE_GUI_plot MATLAB code for RIDE_GUI_plot.fig
%      RIDE_GUI_plot, by itself, creates a new RIDE_GUI_plot or raises the existing
%      singleton*.
%
%      H = RIDE_GUI_plot returns the handle to a new RIDE_GUI_plot or the handle to
%      the existing singleton*.
%
%      RIDE_GUI_plot('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIDE_GUI_plot.M with the given input arguments.
%
%      RIDE_GUI_plot('Property','Value',...) creates a new RIDE_GUI_plot or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RIDE_GUI_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RIDE_GUI_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RIDE_GUI_plot

% Last Modified by GUIDE v2.5 24-Jan-2015 18:36:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RIDE_GUI_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @RIDE_GUI_plot_OutputFcn, ...
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


% --- Executes just before RIDE_GUI_plot is made visible.
function RIDE_GUI_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RIDE_GUI_plot (see VARARGIN)

% Choose default command line output for RIDE_GUI_plot
handles.output = hObject;
if ~isempty(varargin) 
    handles.results = varargin{1}.ride.results; 
    handles.cfg = varargin{1}.ride.results.cfg; 
end
% Update handles structure
guidata(hObject, handles);

RIDE_plot(handles.results,{'erp',handles.results.cfg.comp.name{:}},handles.chan1,handles.axes1);
axis(handles.axes1,'tight');
xlabel(handles.axes1,'time after stimulus (ms)');
ylabel(handles.axes1,'potential');



% UIWAIT makes RIDE_GUI_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RIDE_GUI_plot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




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





% --- Executes on button press in pop_out.
function pop_out_Callback(hObject, eventdata, handles)
% hObject    handle to pop_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = figure;temp2=copyobj(handles.axes1,temp);set(temp2,'Units','normalized');set(temp2,'position',[0.1,0.1,0.8,0.8]);
  


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

