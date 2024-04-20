function varargout = fg_absolute(varargin)
% FG_ABSOLUTE MATLAB code for fg_absolute.fig
%      FG_ABSOLUTE, by itself, creates a new FG_ABSOLUTE or raises the existing
%      singleton*.
%
%      H = FG_ABSOLUTE returns the handle to a new FG_ABSOLUTE or the handle to
%      the existing singleton*.
%
%      FG_ABSOLUTE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_ABSOLUTE.M with the given input arguments.
%
%      FG_ABSOLUTE('Property','Value',...) creates a new FG_ABSOLUTE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_absolute_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_absolute_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_absolute

% Last Modified by GUIDE v2.5 25-Aug-2017 14:26:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_absolute_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_absolute_OutputFcn, ...
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


% --- Executes just before fg_absolute is made visible.
function fg_absolute_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_absolute (see VARARGIN)
    
    % load the current directory project's database
    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    load(fullfile(vars.working_folder,'database.mat'))
    
    if ~isempty(database.benchmarks)
        set(handles.lstAllBench,'string',sort({database.benchmarks.name}));
    end
    
    if ~isempty(database.agravbench)
        set(handles.lstAbsBenchmarks,'string',sort({database.agravbench.name}));
    end

    % Choose default command line output for fg_absolute
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % save the database to the local handle
    hObj = guidata(hObject);
    hObj.database = database;
    hObj.working_folder = vars.working_folder;
    
    hObj.Visualization = CssVisualization([], database.benchmarks, handles.axMap, [.5,.5,.5], [1 0 1], [.2,.2,.2], [1 0 0]);
    hObj.Visualization = hObj.Visualization.PlotBenchmarks();
    
    guidata(hObject, hObj);
    
    % UIWAIT makes fg_absolute wait for user response (see UIRESUME)
    % uiwait(handles.fg_absolute);


% --- Outputs from this function are returned to the command line.
function varargout = fg_absolute_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


% --- Executes on selection change in lstAbsBenchmarks.
function lstAbsBenchmarks_Callback(hObject, eventdata, handles)
    % hObject    handle to lstAbsBenchmarks (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns lstAbsBenchmarks contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from lstAbsBenchmarks

    selection = get(hObject,{'string','value'});
    
    vars = guidata(hObject);
    database = vars.database;
    
    % erase any selection
    % erase_selection(database);
    
    selection = selection{1}(selection{2});
    
    for i = 1:length(selection)
        benchmark = CssBenchmark.ReturnBenchmark(database.agravbench,selection{i});

        if ~isempty(benchmark)
            % select_benchmark(benchmark)
            hObj = guidata(hObject);

            hObj.Visualization = hObj.Visualization.SelectBenchmark(benchmark);

            guidata(hObject, hObj);
    
            set(handles.txtAbsGrav,'String',num2str(benchmark.absolute_g))
            set(handles.txtUncert,'String',num2str(benchmark.uncertainty))
        end
        
    end
    
    set(handles.lstAllBench,'Value',[])
    
    

% --- Executes during object creation, after setting all properties.
function lstAbsBenchmarks_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lstAbsBenchmarks (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function txtAbsGrav_Callback(hObject, eventdata, handles)
    % hObject    handle to txtAbsGrav (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtAbsGrav as text
    %        str2double(get(hObject,'String')) returns contents of txtAbsGrav as a double


% --- Executes during object creation, after setting all properties.
function txtAbsGrav_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to txtAbsGrav (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function txtUncert_Callback(hObject, eventdata, handles)
    % hObject    handle to txtUncert (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtUncert as text
    %        str2double(get(hObject,'String')) returns contents of txtUncert as a double


% --- Executes during object creation, after setting all properties.
function txtUncert_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to txtUncert (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    vars = guidata(hObject);
    database = vars.database;
    save(fullfile(vars.working_folder,'database.mat'), 'database')
    

% --- Executes on button press in cmdCancel.
function cmdCancel_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdCancel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in lstAllBench.
function lstAllBench_Callback(hObject, eventdata, handles)
    % hObject    handle to lstAllBench (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns lstAllBench contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from lstAllBench

    selection = get(hObject,{'string','value'});
    
    vars = guidata(hObject);
    database = vars.database;
    
    selection = selection{1}(selection{2});
    
    for i = 1:length(selection)
        benchmark = CssBenchmark.ReturnBenchmark(database.benchmarks,selection{i});

        if ~isempty(benchmark)
            hObj = guidata(hObject);
            hObj.Visualization = hObj.Visualization.SelectBenchmark(benchmark);
            guidata(hObject, hObj);
        end
    end
    
    

% --- Executes during object creation, after setting all properties.
function lstAllBench_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstAllBench (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdAdd.
function cmdAdd_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdAdd (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % button to add a benchmark to the absolute gravity list
    vars = guidata(hObject);
    database = vars.database;
    
    selection = get(handles.lstAllBench,{'string','value'});
    selection = selection{1}(selection{2});
    
    if ~isempty(selection)
    
        AGravBenchmarks = database.agravbench;
        Benchmarks      = database.benchmarks;
        
        for i = 1:length(selection)
            % check if the selected benchmarks are already in the list
            if ~isempty(CssBenchmark.ReturnBenchmark(AGravBenchmarks, selection(i)))
                % in list, do not add
                continue
            else
                % not in list, add
                benchmark = CssBenchmark.ReturnBenchmark(Benchmarks, selection(i));
                benchmark.absolute_g  = 977452.188;
                benchmark.uncertainty = 0.002;
                database.agravbench = [database.agravbench; benchmark];
            end
        end
        
        set(handles.lstAbsBenchmarks,'string',sort({database.agravbench.name}));
        
        hObj = guidata(hObject);
        hObj.database = database;
        guidata(hObject, hObj);
    
    end

% --- Executes on button press in cmdRemove.
function cmdRemove_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function erase_selection(database)
    if ~isempty(database.benchmarks)
        lat = [database.benchmarks.lat];
        lon = [database.benchmarks.lon];

        m_plot(lon, lat, 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'w')
    end
    
    if ~isempty(database.agravbench)
        lat = [database.agravbench.lat];
        lon = [database.agravbench.lon];
    
        m_plot(lon, lat, '^r', 'MarkerSize', 6, 'MarkerFaceColor', 'c')
    end

function select_benchmark(benchmark)

    if isempty(benchmark.absolute_g)
        m_plot(benchmark.lon, benchmark.lat, 'ob', 'MarkerSize', 3, 'MarkerFaceColor', 'w')
    else
        m_plot(benchmark.lon, benchmark.lat, '^b', 'MarkerSize', 6, 'MarkerFaceColor', 'c')
    end

% --- Executes on button press in cmdSelectMap.
function cmdSelectMap_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSelectMap (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    vars = guidata(hObject);
    
    % request user's input
    [x, y] = ginput(1);
    
    vars.Visualization = vars.Visualization.PickBenchmark(x, y);
    
    if ~isempty(vars.Visualization.selected_benchmark)
        benchmarks = sort({vars.database.benchmarks.name});
        
        index = find(ismember(benchmarks,{vars.Visualization.benchmarks(vars.Visualization.selected_benchmark).name}));
        
        set(handles.lstAllBench,'Value',index)
    end
    
    % save data
    guidata(hObject, vars);
    

% --- Executes on button press in cmdSaveValue.
function cmdSaveValue_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSaveValue (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    selection = get(handles.lstAbsBenchmarks,{'string','value'});
    
    vars = guidata(hObject);
    database = vars.database;

    selection = selection{1}(selection{2});
    
    if ~isempty(selection)
        % save the values to the current benchmark
        benchmark = CssBenchmark.ReturnBenchmark(database.agravbench,selection{1});

        S = get(handles.txtAbsGrav,'String');
        
        if ~all(ismember(S, '0123456789+-.eEdD'))
            warndlg('Only numbers are accepted in the absolute gravity field')
        else
            benchmark.absolute_g = str2double(S);
        end
        
        S = get(handles.txtUncert,'String');
        
        if ~all(ismember(S, '0123456789+-.eEdD'))
            warndlg('Only numbers are accepted in the uncertainty field')
        else
            benchmark.uncertainty = str2double(S);
        end
    end
    
    hObj = guidata(hObject);
    hObj.database = database;
    guidata(hObject, hObj);
    
