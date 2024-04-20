function varargout = fg_find_benchmark(varargin)
% FG_FIND_BENCHMARK MATLAB code for fg_find_benchmark.fig
%      FG_FIND_BENCHMARK, by itself, creates a new FG_FIND_BENCHMARK or raises the existing
%      singleton*.
%
%      H = FG_FIND_BENCHMARK returns the handle to a new FG_FIND_BENCHMARK or the handle to
%      the existing singleton*.
%
%      FG_FIND_BENCHMARK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_FIND_BENCHMARK.M with the given input arguments.
%
%      FG_FIND_BENCHMARK('Property','Value',...) creates a new FG_FIND_BENCHMARK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_find_benchmark_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_find_benchmark_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_find_benchmark

% Last Modified by GUIDE v2.5 25-Feb-2019 14:42:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_find_benchmark_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_find_benchmark_OutputFcn, ...
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


% --- Executes just before fg_find_benchmark is made visible.
function fg_find_benchmark_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_find_benchmark (see VARARGIN)

    % Choose default command line output for fg_find_benchmark
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fg_find_benchmark wait for user response (see UIRESUME)

    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    
    project = CssProject(vars.working_folder);
    database = project.LoadDatabase();
    
    hObj = guidata(hObject);
    
    % set the status of the window
    hObj.status = false;
    hObj.modified_benchmarks = [];
    hObj.project = project;
    hObj = UpdateBenchmarkList(database.benchmarks, hObj, handles);
    hObj.working_folder = vars.working_folder;
    
    fp = get(hObject,'position');
    fp = idgetnicedialoglocation(fp, get(hObject,'Units'));
    set(hObject,'position',fp)
    
    if isempty(varargin)
        handles.cmdSelect.String = 'Save end Exit';
        handles.chkHide.Visible = 'off';
    else
        hObj.search_lat = varargin{1};
        hObj.search_lon = varargin{2};
        if and(hObj.search_lat == 0, hObj.search_lon == 0)
            handles.chkHide.Visible = 'off';
        end
        handles.cmdDeleteBenckmark.Visible = 'off';
    end
    
    guidata(hObject, hObj);
        
    uicontrol(handles.lstBenchmarks)
    uiwait(handles.figure1);
        
function hObj = UpdateBenchmarkList(in_benchmarks, hObj, handles, varargin)

    if handles.chkHide.Value
        % user has requested to filter using the search lat lon parameters
        % passed from the parent window
        d = (m_idist([in_benchmarks.lon],[in_benchmarks.lat],hObj.search_lon,hObj.search_lat)./1000)';
    else
        d = 100*ones(length(in_benchmarks),1);
    end
    
    for i = 1:length(in_benchmarks)
        if d(i) <= 80 % km from the search point
            benchmarks(i,1) = {sprintf('<HTML><strong><FONT color="%s">%s (%.0f km)', 'green', in_benchmarks(i).name, d(i))};
        elseif ~isempty(in_benchmarks(i).absolute_g)
            benchmarks(i,1) = {sprintf('<HTML><strong><FONT color="%s">%s (absolute)', 'blue', in_benchmarks(i).name)};
        else
            benchmarks(i,1) = {in_benchmarks(i).name};
        end
    end
    
    if handles.chkHide.Value
        [~,c] = sortrows([{in_benchmarks.name}' num2cell(d)],[2 1]);
        hObj.benchmarks = in_benchmarks(c);
    else
        [~,c] = sortrows({in_benchmarks.name}');
        hObj.benchmarks = in_benchmarks(c);
    end
    
    set(handles.lstBenchmarks, 'string', benchmarks(c))
    
    if ~isempty(varargin)
        % set the selected benchmark
        
        % find the index of this benchmark
        index = find(ismember({in_benchmarks(c).name},{varargin{1}}));

        set(handles.lstBenchmarks,'Value',index)
    end
    
% --- Outputs from this function are returned to the command line.
function varargout = fg_find_benchmark_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    hObj = guidata(hObject);
    
    % determine the exist state
    if hObj.status == true
        selection = get(handles.lstBenchmarks,{'String','Value'});

        varargout{1} = selection{2};
        varargout{2} = hObj.benchmarks(selection{2});
        
        handles.figure1.Pointer = 'watch';
        drawnow;
        if ~isempty(hObj.modified_benchmarks)
            UpdateBenchmarks(hObj.project, hObj.modified_benchmarks);
        end
        handles.figure1.Pointer = 'arrow';
    else
        varargout{1} = [];
        varargout{2} = [];
    end
    
    delete(handles.figure1)


function txtSearch_Callback(hObject, eventdata, handles)
    % hObject    handle to txtSearch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtSearch as text
    %        str2double(get(hObject,'String')) returns contents of txtSearch as a double


% --- Executes during object creation, after setting all properties.
function txtSearch_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to txtSearch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in lstBenchmarks.
function lstBenchmarks_Callback(hObject, eventdata, handles)
    % hObject    handle to lstBenchmarks (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns lstBenchmarks contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from lstBenchmarks
    
    hObj = guidata(hObject);
    
    Click_lstBenchmarks(handles, hObj)
    
    if strcmp(handles.figure1.SelectionType,'open')
        % if double click, trigger selection
        hObj.status = true;
        guidata(hObject, hObj);
        uiresume(handles.figure1);
    end
    
    
function Click_lstBenchmarks(handles, hObj)
    selection = get(handles.lstBenchmarks,{'string','value'});
    
    hObj.selected_benchmark = selection{2};
    
    set(handles.txtLat, 'string', sprintf('%.8f',hObj.benchmarks(selection{2}).lat))
    set(handles.txtLon, 'string', sprintf('%.8f',hObj.benchmarks(selection{2}).lon))
    set(handles.txtHeight, 'string', sprintf('%.3f',hObj.benchmarks(selection{2}).height))
    set(handles.txtX, 'string', sprintf('%.3f',hObj.benchmarks(selection{2}).x))
    set(handles.txtY, 'string', sprintf('%.3f',hObj.benchmarks(selection{2}).y))
    set(handles.txtZ, 'string', sprintf('%.3f',hObj.benchmarks(selection{2}).z))
    set(handles.BenchmarkPanel, 'title', ['Benchmark Metadata: ' hObj.benchmarks(selection{2}).name])


% --- Executes during object creation, after setting all properties.
function lstBenchmarks_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lstBenchmarks (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in cmdSelect.
function cmdSelect_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSelect (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    hObj.status = true;
    
    guidata(hObject, hObj);
    
    uiresume(handles.figure1);

% --- Executes on button press in cmdNew.
function cmdNew_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdNew (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    
    name = inputdlg('Enter the new benchmark name:', 'New Benchmark', [1 50]);
    
    if ~isempty(name)
        name = lower(name{1});
        
        % verify that the benchmark doesn' exist
        if CssBenchmark.exists(hObj.benchmarks, name)
            h = warndlg('The name entered already exists in the database');
            waitfor(h);
        else
            % the name doesn't exist, then request the RINEX data
            % initialize in zeros
            hObj.benchmarks(end+1) = CssBenchmark(name, 0, 0, 0, 0);
            % save the new benchmark to know which ones to process
            hObj.modified_benchmarks = [hObj.modified_benchmarks; hObj.benchmarks(end)];
        end
        
        % update the list and set the selected benchmark
        hObj = UpdateBenchmarkList(hObj.benchmarks, hObj, handles, name);
        % select it
        Click_lstBenchmarks(handles, hObj);
            
        guidata(hObject, hObj);
    end


function txtLat_Callback(hObject, eventdata, handles)
    % hObject    handle to txtLat (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtLat as text
    %        str2double(get(hObject,'String')) returns contents of txtLat as a double
    
    % trigger the conversion button callback
    cmdLLA2XYZ_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function txtLat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLon_Callback(hObject, eventdata, handles)
    % hObject    handle to txtLon (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtLon as text
    %        str2double(get(hObject,'String')) returns contents of txtLon as a double
    
    % trigger the conversion button callback
    cmdLLA2XYZ_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function txtLon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtHeight_Callback(hObject, eventdata, handles)
    % hObject    handle to txtHeight (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtHeight as text
    %        str2double(get(hObject,'String')) returns contents of txtHeight as a double

    % trigger the conversion button callback
    cmdLLA2XYZ_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function txtHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtRINEXPath_Callback(hObject, eventdata, handles)
    % hObject    handle to txtRINEXPath (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtRINEXPath as text
    %        str2double(get(hObject,'String')) returns contents of txtRINEXPath as a double


% --- Executes during object creation, after setting all properties.
function txtRINEXPath_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to txtRINEXPath (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function txtX_Callback(hObject, eventdata, handles)
    % hObject    handle to txtX (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtX as text
    %        str2double(get(hObject,'String')) returns contents of txtX as a double
    
    cmdXYZ2LLA_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function txtX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtY_Callback(hObject, eventdata, handles)
    % hObject    handle to txtY (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtY as text
    %        str2double(get(hObject,'String')) returns contents of txtY as a double

    % trigger the conversion button callback
    cmdXYZ2LLA_Callback(hObject, eventdata, handles);
    
% --- Executes during object creation, after setting all properties.
function txtY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZ_Callback(hObject, eventdata, handles)
    % hObject    handle to txtZ (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtZ as text
    %        str2double(get(hObject,'String')) returns contents of txtZ as a double

    cmdXYZ2LLA_Callback(hObject, eventdata, handles);
    
% --- Executes during object creation, after setting all properties.
function txtZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdLLA2XYZ.
function cmdLLA2XYZ_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdLLA2XYZ (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    lat = str2double(handles.txtLat.String)*pi/180;
    lon = str2double(handles.txtLon.String)*pi/180;
    h   = str2double(handles.txtHeight.String);
    
    xyz = CssBenchmark.e2x([lat; lon; h]);
    
    handles.txtX.String = sprintf('%.3f',xyz(1));
    handles.txtY.String = sprintf('%.3f',xyz(2));
    handles.txtZ.String = sprintf('%.3f',xyz(3));
    
    hObj = guidata(hObject);
    % get the selection
    selection = get(handles.lstBenchmarks,{'String','Value'});
    % apply modification to benchmark
    hObj.benchmarks(selection{2}) = hObj.benchmarks(selection{2}).UpdateCoords(CssBenchmark('dummy', lat*180/pi, lon*180/pi, h, 0));
    % save in the modified list
    hObj = save_modified_benchmark(hObj, selection{2});
    
    guidata(hObject, hObj);

% --- Executes on button press in cmdXYZ2LLA.
function cmdXYZ2LLA_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdXYZ2LLA (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    x = str2double(handles.txtX.String);
    y = str2double(handles.txtY.String);
    z = str2double(handles.txtZ.String);
    
    lla = CssBenchmark.x2e([x; y; z]);
    
    handles.txtLat.String = sprintf('%.8f',lla(1)*180/pi);
    handles.txtLon.String = sprintf('%.8f',lla(2)*180/pi);
    handles.txtHeight.String = sprintf('%.3f',lla(3));
    
    hObj = guidata(hObject);
    % get the selection
    selection = get(handles.lstBenchmarks,{'String','Value'});
    % apply modification to benchmark
    hObj.benchmarks(selection{2}) = hObj.benchmarks(selection{2}).UpdateCoords(CssBenchmark('dummy', x, y, z, 0));
    % save in the modified list
    hObj = save_modified_benchmark(hObj, selection{2});
    
    guidata(hObject, hObj);
    
    
function hObj = save_modified_benchmark(hObj, selection)    

    % save the new benchmark to know which ones to process
    if ~CssBenchmark.exists(hObj.modified_benchmarks, hObj.benchmarks(selection).name)
        hObj.modified_benchmarks = [hObj.modified_benchmarks; hObj.benchmarks(selection)];
    end

    
% --- Executes on button press in cmdRINEXPath.
function cmdRINEXPath_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdRINEXPath (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    if ~strcmp(handles.txtRINEXPath.String,'')
        [path, ~, ~] = fileparts(handles.txtRINEXPath.String);
        
        [filename, pathname] = uigetfile(fullfile(path, '*.??d.Z;*.??o;*.??O'),'Find RINEX file');
    else
        [filename, pathname] = uigetfile('*.??d.Z;*.??o;*.??O','Find RINEX file');
    end
    
    if filename ~= 0
        handles.txtRINEXPath.String = fullfile(pathname, filename);
        % open the rinex and get the antenna height
        [path,rinex] = CssBenchmark.CopyCrinexRinex('temp',handles.txtRINEXPath.String);
        [~, ~, ~, height] = CssBenchmark.ReadRinex(rinex);
        % delete temp dir
        delete(rinex)
        rmdir(path)
        handles.txtAntH.String = sprintf('%.3f', height);
    end

% --- Executes on button press in cmdPPP.
function cmdPPP_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdPPP (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    if exist(handles.txtRINEXPath.String, 'file')
        
        hObj = guidata(hObject);
        
        h = waitbar(0,'Obtaining coordinate using NRCAN PPP...');
    
        selection = get(handles.lstBenchmarks,{'String','Value'});

        hObj.benchmarks(selection{2}).RunPPP(handles.txtRINEXPath.String, str2double(handles.txtAntH.String));
        
        Click_lstBenchmarks(handles, hObj)
        
        close(h)
        
        % save the new benchmark to know which ones to process
        hObj = save_modified_benchmark(hObj, selection{2});
        
        guidata(hObject,hObj);
    else
        warndlg('The provided RINEX file could not be found!')
    end
    


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: delete(hObject) closes the figure
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        uiresume(hObject)
    else
        delete(hObject);
    end



function txtAntH_Callback(hObject, eventdata, handles)
% hObject    handle to txtAntH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAntH as text
%        str2double(get(hObject,'String')) returns contents of txtAntH as a double


% --- Executes during object creation, after setting all properties.
function txtAntH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAntH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkHide.
function chkHide_Callback(hObject, eventdata, handles)
    % hObject    handle to chkHide (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of chkHide
    hObj = guidata(hObject);

    hObj = UpdateBenchmarkList(hObj.benchmarks, hObj, handles);
    
    Click_lstBenchmarks(handles, hObj)
    
    guidata(hObject,hObj);
    


% --- Executes on button press in cmdDeleteBenckmark.
function cmdDeleteBenckmark_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdDeleteBenckmark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    
    project = CssProject(hObj.working_folder);
    database = project.LoadDatabase();
    
    selection = get(handles.lstBenchmarks,{'String','Value'});
    name = hObj.benchmarks(selection{2}).name;
    
    lines = project.LoadLines();
    
    for j = 1:length(lines)
        if CssBenchmark.exists(lines(j).benchmarks, name)
            h = msgbox(['This benchmark is in use in line ' lines(j).line_name '. Remove it and try again!']);
            waitfor(h)
            return
        end
    end
    
    i = ismember({database.benchmarks.name}, name);
    
    database.benchmarks(i) = [];
    
     % set the status of the window
    hObj.status = false;
    hObj.project = project;
    hObj = UpdateBenchmarkList(database.benchmarks, hObj, handles);

    save(fullfile(hObj.working_folder, 'database.mat'), 'database');
    
    guidata(hObject, hObj);
    
    