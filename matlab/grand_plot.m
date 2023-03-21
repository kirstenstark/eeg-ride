function varargout = grand_plot(varargin)
% grand_plot MATLAB code for grand_plot.fig
%      grand_plot, by itself, creates a new grand_plot or raises the existing
%      singleton*.
%
%      H = grand_plot returns the handle to a new grand_plot or the handle to
%      the existing singleton*.
%
%      grand_plot('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in grand_plot.M with the given input arguments.
%
%      grand_plot('Property','Value',...) creates a new grand_plot or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before grand_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to grand_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help grand_plot

% Last Modified by GUIDE v2.5 20-Jan-2015 17:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grand_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @grand_plot_OutputFcn, ...
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


% --- Executes just before grand_plot is made visible.
function grand_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to grand_plot (see VARARGIN)

% Choose default command line output for grand_plot
handles.output = hObject;
if ~isempty(varargin) handles.results = varargin{1}; end
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes grand_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = grand_plot_OutputFcn(hObject, eventdata, handles) 
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



% --- Executes on button press in load_data1.
function load_data1_Callback(hObject, eventdata, handles)
% hObject    handle to load_data1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



    [handles.a{1},handles.b{1}] = uigetfile('*.mat','MultiSelect','on');rt = [];
    if ~iscell(handles.a{1}) handles.a{1} = {handles.a{1}};end
    for j = 1:length(handles.a{1})
        load([handles.b{1},handles.a{1}{j}]);
        if j==1 
            handles.cfg = results.cfg;
            comp_item = {'erp',handles.cfg.comp.name{:}};
            for jj = 1:length(handles.cfg.comp.name)
                comp_item{end+1} = [handles.cfg.comp.name{jj},'_sl'];
            end
            comp_item{end+1} = 'erp_new';
        end
    end
    
    set(handles.object_popup,'String',comp_item);
    handles.comp_item = comp_item;
    guidata(hObject,handles);



% --- Executes on button press in pop_out.
function pop_out_Callback(hObject, eventdata, handles)
% hObject    handle to pop_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = figure;temp2=copyobj(handles.axes1,temp);set(temp2,'Units','normalized');set(temp2,'position',[0.1,0.1,0.8,0.8]);



% --- Executes on button press in load_data2.
function load_data2_Callback(hObject, eventdata, handles)
% hObject    handle to load_data2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



    [handles.a{2},handles.b{2}] = uigetfile('*.mat','MultiSelect','on');rt = [];
    if ~iscell(handles.a{2}) handles.a{2} = {handles.a{2}};end
    for j = 1:length(handles.a{2})
        load([handles.b{2},handles.a{2}{j}]);
        if j==1 
            handles.cfg = results.cfg;
            comp_item = {'erp',handles.cfg.comp.name{:}};
            for jj = 1:length(handles.cfg.comp.name)
                comp_item{end+1} = [handles.cfg.comp.name{jj},'_sl'];
            end
            comp_item{end+1} = 'erp_new';
        end
    end
    
    set(handles.object_popup,'String',comp_item);
    handles.comp_item = comp_item;
    guidata(hObject,handles);




% --- Executes on selection change in items2plot.
function items2plot_Callback(hObject, eventdata, handles)
% hObject    handle to items2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns items2plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from items2plot
handles.plot_item = get(hObject,'Value');guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.
function items2plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to items2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.object_popup = hObject;guidata(hObject,handles);


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);hold off;
if ~isfield(handles,'plot_item') handles.plot_item = 1;end

if isfield(handles,'early_t1') handles = rmfield(handles,'early_t1');end
if isfield(handles,'early_t2') handles = rmfield(handles,'early_t2');end

if isfield(handles,'a')
    colors = {'b','r'};cons = {'cond1','cond2'};
    for p = 1:length(handles.a)
        if ~isempty(handles.a{p})
            a = handles.a{p};b = handles.b{p};
            rt = [];
            for j = 1:length(a)
                load([b,a{j}]);
                comp_item = handles.comp_item;
                for c = 1:length(comp_item)
                    eval(['comp_',comp_item{c},'(:,:,j) = results.',comp_item{c},';']);
                    if ~isempty(find(strcmpi(comp_item,'r'))) rt(j) = mean(results.cfg.comp.latency{end});end
                end
            end
            
            %re_syn R to RT
            if ~isempty(find(strcmpi(comp_item,'r'))) comp_r = move3(comp_r,round((median(rt)-rt)/handles.cfg.samp_interval));end
            
            
            %re_syn ERP to RT
            if strcmpi(comp_item(handles.plot_item),'erp_new') && ~isempty(find(strcmpi(comp_item,'r')))
                comp_erp_new = zeros(size(comp_s));
                for c = 1:length(handles.cfg.comp.name)
                    eval(['comp_erp_new = comp_erp_new + comp_',handles.cfg.comp.name{c},';']);
                end
            end
            
            
            for c = 1:length(comp_item)
                eval(['handles.comp_',comp_item{c},' = mean(comp_',comp_item{c},',3);']);
            end
            
            eval(['temp = handles.comp_',comp_item{handles.plot_item},'(:,handles.chan1);']);
            
            axes(handles.axes1);
            if ~strcmpi(comp_item(handles.plot_item),'r')
                plot(handles.axes1,linspace(handles.cfg.epoch_twd(1),handles.cfg.epoch_twd(2),size(handles.comp_erp,1)),temp,'color',colors{p});hold on;axis tight;
                xlabel('time after stimulus (ms)');ylabel('potential (\muV))');
            else
                plot(handles.axes1,linspace(handles.cfg.epoch_twd(1),handles.cfg.epoch_twd(2),size(handles.comp_erp,1)),temp,'color',colors{p});hold on;axis tight;
                xlabel('time after RT (ms)');ylabel('potential (\muV))');
            end
            
            
            
%             if ~isempty(find(strcmpi(comp_item,'r')))
%                 plot(handles.axes1,linspace(handles.cfg.epoch_twd(1),handles.cfg.epoch_twd(2),size(handles.comp_erp,1)),temp,'color',colors{p});hold on;axis tight;
%                 xlabel('time after stimulus (ms)');ylabel('potential (\muV))');
%             else
%                 plot(handles.axes1,linspace(handles.cfg.epoch_twd(1),handles.cfg.epoch_twd(2),size(handles.comp_erp,1)),temp,'color',colors{p});hold on;axis tight;
%                 xlabel('time after RT (ms)');ylabel('potential (\muV))');
%             end
            
            for c = 1:length(comp_item)
                eval(['handles.cps{p}{c} = comp_',comp_item{c},';']);
            end
        end
            
    end
    legend(handles.axes1,cons(1:p));hold off;
    guidata(hObject,handles);
end



% --- Executes on button press in clear_button.
function clear_button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'a') handles = rmfield(handles,'a');end
if isfield(handles,'b') handles = rmfield(handles,'b');end

guidata(hObject,handles);

plot(handles.axes1,0);hold off;



function ttest_twd_Callback(hObject, eventdata, handles)
% hObject    handle to ttest_twd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttest_twd as text
%        str2double(get(hObject,'String')) returns contents of ttest_twd as a double


temp = textscan(get(hObject,'String'),'%f'); 
handles.t_twd = temp{1}';
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function ttest_twd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttest_twd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ttest_button.
function ttest_button_Callback(hObject, eventdata, handles)
% hObject    handle to ttest_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temp_twd = round((handles.t_twd-handles.cfg.epoch_twd(1))/handles.cfg.samp_interval);


    colors = {'b','r'};cons = {'cond1','cond2'};
    for p = 1:length(handles.a)

        temp(:,p) = squeeze(mean(handles.cps{p}{handles.plot_item}(temp_twd,handles.chan1,:)));
        
    end
    
    [a,b,c,d] = ttest(temp(:,1)-temp(:,2));
    x = get(handles.axes1,'xlim');y = get(handles.axes1,'ylim');
    axes(handles.axes1);
    
    x1 = mean(x);y1 = y(1) + (y(2)-y(1))/5;
    x2 = mean(x);y2 = y(1) + (y(2)-y(1))/10;
    
    
    
    if isfield(handles,'early_t1')
        set(handles.early_t1,'Position',[x1,y1],'String',['t=',num2str(d.tstat)]);
    else
        handles.early_t1 = text(x1,y1,['t=',num2str(d.tstat)]);
    end
     if isfield(handles,'early_t2')
        set(handles.early_t2,'Position',[x2,y2],'String',['p=',num2str(b)]);
     else
         handles.early_t2 = text(x2,y2,['p=',num2str(b)]);
     end
    
     
     
   axes(handles.axes1);
   temp = findobj(gca,'type','image');for j = 1:length(temp) delete(temp(j));end
   hold on;h=imagesc(linspace(handles.t_twd(1),handles.t_twd(2),100),linspace(y(1),y(2),100),zeros(100,100));alpha(0.4);
   uistack(h,'bottom');
    
    guidata(hObject,handles);
    
    
        





function chan2_Callback(hObject, eventdata, handles)
% hObject    handle to chan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chan2 as text
%        str2double(get(hObject,'String')) returns contents of chan2 as a double
handles.chan2 = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function chan2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.chan2 = str2double(get(hObject,'String'));%disp(str2double(get(hObject,'String')));
guidata(hObject, handles);


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in items2plot2.
function items2plot2_Callback(hObject, eventdata, handles)
% hObject    handle to items2plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns items2plot2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from items2plot2
handles.plot_item2 = get(hObject,'Value');guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function items2plot2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to items2plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.object_popup2 = hObject;guidata(hObject,handles);


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in unp_ttest.
function unp_ttest_Callback(hObject, eventdata, handles)
% hObject    handle to unp_ttest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


temp_twd = round((handles.t_twd-handles.cfg.epoch_twd(1))/handles.cfg.samp_interval);


    colors = {'b','r'};cons = {'cond1','cond2'};
    for p = 1:length(handles.a)

        temp(:,p) = squeeze(mean(handles.cps{p}{handles.plot_item}(temp_twd,handles.chan1,:)));
        
    end
    
    [a,b,c,d] = ttest2(temp(:,1),temp(:,2));
    x = get(handles.axes1,'xlim');y = get(handles.axes1,'ylim');
    axes(handles.axes1);
    
    x1 = mean(x);y1 = y(1) + (y(2)-y(1))/5;
    x2 = mean(x);y2 = y(1) + (y(2)-y(1))/10;
    
    
    
    if isfield(handles,'early_t1')
        set(handles.early_t1,'Position',[x1,y1],'String',['t=',num2str(d.tstat)]);
    else
        handles.early_t1 = text(x1,y1,['t=',num2str(d.tstat)]);
    end
     if isfield(handles,'early_t2')
        set(handles.early_t2,'Position',[x2,y2],'String',['p=',num2str(b)]);
     else
         handles.early_t2 = text(x2,y2,['p=',num2str(b)]);
    end
    
    
    
    guidata(hObject,handles);
