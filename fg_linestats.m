function varargout = fg_linestats(varargin)
% FG_LINESTATS MATLAB code for fg_linestats.fig
%      FG_LINESTATS, by itself, creates a new FG_LINESTATS or raises the existing
%      singleton*.
%
%      H = FG_LINESTATS returns the handle to a new FG_LINESTATS or the handle to
%      the existing singleton*.
%
%      FG_LINESTATS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_LINESTATS.M with the given input arguments.
%
%      FG_LINESTATS('Property','Value',...) creates a new FG_LINESTATS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_linestats_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_linestats_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_linestats

% Last Modified by GUIDE v2.5 26-Sep-2020 11:59:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_linestats_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_linestats_OutputFcn, ...
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


% --- Executes just before fg_linestats is made visible.
function fg_linestats_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_linestats (see VARARGIN)

    % load the current directory project's database
    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    
    project = CssProject(vars.working_folder);
    
    database = project.LoadDatabase();
    lines    = project.LoadLines();
      
    residual_std = EstimateStd(lines);
    set(handles.lblStd,'String',[num2str(3*residual_std) ' mGal'])
    set(handles.lblLines,'String',num2str(length(lines)))
    set(handles.lblBenchmarks,'String',num2str(length(database.benchmarks)))
    
    % Choose default command line output for fg_absolute
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % save the database to the local handle
    hObj                = guidata(hObject);
    hObj.database       = database;
    hObj.working_folder = vars.working_folder;
    hObj.Lines          = lines;
    hObj.changed_lines  = [];
    
    % visualization object
    hObj.Visualization = CssVisualization(lines, database.benchmarks, handles.axMap, [.5,.5,.5], [1 0 1], [.2,.2,.2], [1 0 0]);
    hObj.Visualization = hObj.Visualization.PlotLines(residual_std);
    hObj.Visualization = hObj.Visualization.PlotNodes();
    
    hObj.std           = residual_std;
    
    guidata(hObject, hObj);
    
    % fill the listbox
    set(handles.lstLines,'String', hObj.Visualization.LineNameList)

    % UIWAIT makes fg_linestats wait for user response (see UIRESUME)
    % uiwait(handles.fg_linestats);

function residual_std = EstimateStd(lines)
    % calculate the residual standard deviation based on the line info

    residual_std = [];
    for i = 1:length(lines)
        res = cell2mat(lines(i).residuals);
        if any(isnan(res))
            disp(['Line ' lines(i).line_name ' has residuals equal to NaN. Check the line to identify the problem.'])
        else
            residual_std = [residual_std; res];
        end
    end
    residual_std = std(residual_std);

% --- Outputs from this function are returned to the command line.
function varargout = fg_linestats_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

% --- Executes on selection change in lstLines.
function lstLines_Callback(hObject, eventdata, handles)
    % hObject    handle to lstLines (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns lstLines contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from lstLines

    vars = guidata(hObject);
    
    selection = get(handles.lstLines,{'string','value'});
    line_index = selection{2};
        
    if strcmp(get(handles.fg_linestats,'SelectionType'),'open')
        
        % double click on a line: open stats
        if handles.optTime.Value ~= 0
            vars.Visualization.PlotLineTime(vars.selected_line, true)
        elseif handles.optDeltas.Value ~= 0
            vars.Visualization.PlotLineInfo(vars.selected_line, true)
        elseif handles.optCompare.Value ~= 0
            vars.Visualization.PlotComparison(vars.selected_line, true)
        end
    else
        % simple clic, show information about instruments
        ShowInstruments(vars, line_index, handles)
    end
    
    hObj               = guidata(hObject);
    hObj.selected_line = line_index;
    hObj.Visualization = hObj.Visualization.SelectLine(line_index);
    
    guidata(hObject, hObj);

function ShowInstruments(vars, line_index, handles)

    % check the state of this instrument (check the state of the
    % benchmarks
    for i = 1:length(vars.Lines(line_index).instruments)
        if all(vars.Lines(line_index).status(i,:) == 0)
            inst_names(i) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', vars.Lines(line_index).instruments{i})};
        else
            inst_names(i) = vars.Lines(line_index).instruments(i);
        end
    end
    set(handles.lstInstruments,'String', inst_names)
    set(handles.lstInstruments,'Value', 1)
    set(handles.lstBenchmarks,'String', {})

% --- Executes during object creation, after setting all properties.
function lstLines_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lstLines (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in lstInstruments.
function lstInstruments_Callback(hObject, eventdata, handles)
    % hObject    handle to lstInstruments (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns lstInstruments contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from lstInstruments
    vars = guidata(hObject);
    
    selection = get(handles.lstInstruments,{'string','value'});
    index = selection{2};
    
    if strcmp(get(handles.fg_linestats,'SelectionType'),'open')
        instrument = vars.Lines(vars.selected_line).instruments{index};
        % double click: disable/enable instrument
        if all(vars.Lines(vars.selected_line).status(index,:) == false)
            % deactivated, activate
            vars.Lines(vars.selected_line) = vars.Lines(vars.selected_line).ActivateObservationPair(instrument, vars.Lines(vars.selected_line).benchmarks);
        else
            % activated, deactivate
            vars.Lines(vars.selected_line) = vars.Lines(vars.selected_line).DeactivateObservationPair(instrument, vars.Lines(vars.selected_line).benchmarks);
        end
        
        % repaint the item to show the current state
        ShowInstruments(vars, vars.selected_line, handles)
        
        residual_std = EstimateStd(vars.Lines);
        set(handles.lblStd,'String',[num2str(3*residual_std) ' mGal'])
        
        % update the standard deviation and lines
        hObj        = guidata(hObject);
        hObj.std    = residual_std;
        hObj.Lines  = vars.Lines;
        hObj.Visualization.lines = vars.Lines;
        hObj.Visualization = hObj.Visualization.UpdateStdMax(residual_std);
        
        hObj.changed_lines = [hObj.changed_lines; vars.selected_line];
        
        set(handles.lstLines,'String', hObj.Visualization.LineNameList)
        
        % if plots for this line are open, refresh
        if handles.optTime.Value ~= 0
            vars.Visualization.PlotLineTime(vars.selected_line, false)
        elseif handles.optDeltas.Value ~= 0
            vars.Visualization.PlotLineInfo(vars.selected_line, false)
        elseif handles.optCompare.Value ~= 0
            vars.Visualization.PlotComparison(vars.selected_line, false)
        end
        
        guidata(hObject, hObj);
    else
        % simple click, show benchmarks
        ShowBenchmarks(vars, index, handles)
    end
    
    hObj                = guidata(hObject);
    hObj.selected_inst  = index;
    
    guidata(hObject, hObj);

function ShowBenchmarks(vars, inst_index, handles)

    for i = 1:size(vars.Lines(vars.selected_line).status(inst_index,:),2)
        if vars.Lines(vars.selected_line).status(inst_index,i) == false
            benchmarks(i) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', vars.Lines(vars.selected_line).benchmarks(i).name)};
        else
            benchmarks(i) = {vars.Lines(vars.selected_line).benchmarks(i).name};
        end
    end
    set(handles.lstBenchmarks,'String', benchmarks)
    set(handles.lstBenchmarks,'Value', 1)
    
% --- Executes during object creation, after setting all properties.
function lstInstruments_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to lstInstruments (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
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
    vars = guidata(hObject);
    
    selection = get(handles.lstBenchmarks,{'string','value'});
    index = selection{2};
    
    if strcmp(get(handles.fg_linestats,'SelectionType'),'open')
        % deactivate/activate this benchmark
        instrument = vars.Lines(vars.selected_line).instruments{vars.selected_inst};
        % double click: disable/enable instrument
        if vars.Lines(vars.selected_line).status(vars.selected_inst,index) == false
            % deactivated, activate
            vars.Lines(vars.selected_line) = vars.Lines(vars.selected_line).ActivateObservationPair(instrument, vars.Lines(vars.selected_line).benchmarks(index));
        else
            % activated, deactivate
            vars.Lines(vars.selected_line) = vars.Lines(vars.selected_line).DeactivateObservationPair(instrument, vars.Lines(vars.selected_line).benchmarks(index));
        end
        
        % repaint the item to show the current state
        ShowInstruments(vars, vars.selected_line, handles)
        ShowBenchmarks (vars, vars.selected_inst, handles)
        
        residual_std = EstimateStd(vars.Lines);
        set(handles.lblStd,'String',[num2str(3*residual_std) ' mGal'])
        
        % update the standard deviation and lines
        hObj        = guidata(hObject);
        hObj.std    = residual_std;
        hObj.Lines  = vars.Lines;
        hObj.Visualization.lines = vars.Lines;
        hObj.Visualization = hObj.Visualization.UpdateStdMax(residual_std);
        
        hObj.changed_lines = [hObj.changed_lines; vars.selected_line];
        
        set(handles.lstLines,'String', hObj.Visualization.LineNameList)
        
        % if plots for this line are open, refresh
        if handles.optTime.Value ~= 0
            hObj.Visualization.PlotLineTime(vars.selected_line, false)
        elseif handles.optDeltas.Value ~= 0
            hObj.Visualization.PlotLineInfo(vars.selected_line, false)
        elseif handles.optCompare.Value ~= 0
            hObj.Visualization.PlotComparison(vars.selected_line, false)
        end
        
        guidata(hObject, hObj);
    end

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


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
    % hObject    handle to uipushtool1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    [FileName,PathName] = uiputfile({'*.png' ;'*.tif';'*.jpg'},'Save figure as...');

    fr = getframe(handles.axMap);
    imwrite(fr.cdata, fullfile(PathName,FileName))
    
    


% --- Executes on button press in cmdZoomIn.
function cmdZoomIn_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdZoomIn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    hObj = guidata(hObject);
    
    hObj.Visualization.ZoomInLine();
    
    
    


% --- Executes on button press in optDeltas.
function optDeltas_Callback(hObject, eventdata, handles)
% hObject    handle to optDeltas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optDeltas


% --- Executes on button press in optTimes.
function optTimes_Callback(hObject, eventdata, handles)
% hObject    handle to optTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optTimes


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % loop through the lines and save them to disk
    vars = guidata(hObject);
    
    h = waitbar(0,'Saving lines to disk...');
    
    vars.changed_lines = unique(vars.changed_lines);
    
    for i = 1:length(vars.changed_lines)
        vars.Lines(vars.changed_lines(i)).SaveGravityLine(fullfile(vars.working_folder,'lines'), true);
        
        waitbar(i/length(vars.Lines))
    end
    
    close(h)
    


% --- Executes on button press in cmdSelectMap.
function cmdSelectMap_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSelectMap (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    vars = guidata(hObject);
    
    % request user's input
    [x, y] = ginput(1);
    
    vars.Visualization = vars.Visualization.PickLine(x, y);
    
    if ~isempty(vars.Visualization.selected_line)
        set(handles.lstLines,'Value', vars.Visualization.selected_line)
    end
    
    % save data
    guidata(hObject, vars);


% --- Executes on button press in cmdNetworkAnalysis.
function cmdNetworkAnalysis_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdNetworkAnalysis (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    vars = guidata(hObject);
    
    figure;
    set(gcf, 'Name', ['Information for network of project ' vars.database.proj_description]);
    vars.Visualization.network.PlotNetwork(vars.database.agravbench)
    
    


% --- Executes on button press in cmdFind.
function cmdFind_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdFind (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    benchmark = inputdlg('Enter benchmark name','Find benchmark');
    
    vars = guidata(hObject);
    
    for i = 1:length(vars.Lines)
        if ismember(benchmark, {vars.Lines(i).benchmarks.name})
            disp(['found at ' num2str(i)])
            % select the line
            set(handles.lstLines, 'Value', i)
            return
        end
    end
    
    warndlg('Could not find specified benchmark name');
    
    
