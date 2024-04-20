function varargout = fg_adjust(varargin)
% FG_ADJUST MATLAB code for fg_adjust.fig
%      FG_ADJUST, by itself, creates a new FG_ADJUST or raises the existing
%      singleton*.
%
%      H = FG_ADJUST returns the handle to a new FG_ADJUST or the handle to
%      the existing singleton*.
%
%      FG_ADJUST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_ADJUST.M with the given input arguments.
%
%      FG_ADJUST('Property','Value',...) creates a new FG_ADJUST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_adjust_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_adjust_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_adjust

% Last Modified by GUIDE v2.5 26-Sep-2020 13:10:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_adjust_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_adjust_OutputFcn, ...
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


% --- Executes just before fg_adjust is made visible.
function fg_adjust_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_adjust (see VARARGIN)

    % Choose default command line output for fg_adjust
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fg_adjust wait for user response (see UIRESUME)
    % uiwait(handles.fg_adjust);
    
    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    
    project = CssProject(vars.working_folder);
    
    database = project.LoadDatabase();
    lines    = project.LoadLines();
    
    if ~isempty(database.agravbench)
        set(handles.lstAbsBenchmarks,'string',sort({database.agravbench.name}));
    end

        % save the database to the local handle
    hObj                = guidata(hObject);
    hObj.database       = database;
    hObj.working_folder = vars.working_folder;
    hObj.LinesInAdjust  = (1:length(lines))';
    hObj.Network        = network;
    % list of activated abs benchmarks
    hObj.UseConstrains  = database.agravbench;
    % initialize the adjustment object
    hObj.Adjustment = CssAdjustment(lines, database.benchmarks);
    hObj.Adjustment = hObj.Adjustment.GetDesign();
    
    % visualization object
    hObj.Visualization = CssVisualization(lines, database.benchmarks, handles.axMap, [.5,.5,.5], [1 0 1], [.2,.2,.2], [1 0 0]);
    hObj.Visualization = hObj.Visualization.PlotLinesAdjustment();
    hObj.Visualization = hObj.Visualization.PlotBenchmarks();
    
    guidata(hObject, hObj);
    
    set(handles.lstLines,'String', hObj.Visualization.LineNameList)

% --- Outputs from this function are returned to the command line.
function varargout = fg_adjust_OutputFcn(hObject, eventdata, handles) 
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
    
    % check the state of this instrument (check the state of the
    % benchmarks
    
    hObj = guidata(hObject);
    
    selection = get(handles.lstAbsBenchmarks,{'string','value'});
    index = selection{2};
    names = sort({hObj.database.agravbench.name});
    name  = names{index};
    
    if strcmp(get(handles.fg_adjust,'SelectionType'),'open')
        
        if ~ismember(name, {hObj.UseConstrains.name})
            hObj.UseConstrains = [hObj.UseConstrains; CssBenchmark.ReturnBenchmark(hObj.database.agravbench, name)];
        else
            hObj.UseConstrains(ismember({hObj.UseConstrains.name}, name)) = [];
        end
        
        for i = 1:length(names)
            if ~ismember(names{i}, {hObj.UseConstrains.name})
                inst_names(i) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', names{i})};
            else
                inst_names(i) = names(i);
            end
        end
        
        set(handles.lstAbsBenchmarks,'String', inst_names)
    end
    
    % show the benchmark selection
    hObj.Visualization = hObj.Visualization.SelectBenchmark(CssBenchmark.ReturnBenchmark(hObj.Visualization.benchmarks, name));
    
    guidata(hObject, hObj);

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
        
    if strcmp(get(handles.fg_adjust,'SelectionType'),'open')
        
        % double click on a line: remove line
        vars.Visualization = vars.Visualization.ToggleLineActivation(line_index);
        
        set(handles.lstLines,'String', vars.Visualization.LineNameList)
    
    end
    
    vars.selected_line = line_index;
    vars.Visualization = vars.Visualization.SelectLineAdjustment(line_index);
    
    guidata(hObject, vars);

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


% --- Executes on button press in cmdAdjust.
function cmdAdjust_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdAdjust (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    
    if sum(hObj.Visualization.lines_status) ~= length(hObj.Adjustment.lines)
        % redefine the Adjustment object (there was a change in line
        % activation)
        hObj.Adjustment = CssAdjustment(hObj.Visualization.lines(hObj.Visualization.lines_status), hObj.database.benchmarks);
        hObj.Adjustment = hObj.Adjustment.GetDesign();
    end
    
    hObj.Adjustment = hObj.Adjustment.Invert(hObj.UseConstrains);

    hObj.Visualization = hObj.Visualization.PlotResiduals(hObj.Adjustment.residuals, hObj.Adjustment.outliers);
    
    % display statistics
    stats = {['       Variance unit weight: ' sprintf('%.3f', hObj.Adjustment.So)]; 
             ['  Res. st.dev (no outliers): ' sprintf('%.3f', nanstd(hObj.Adjustment.residuals.residual(hObj.Adjustment.outliers)))];
             ['Res. st.dev (with outliers): ' sprintf('%.3f', nanstd(hObj.Adjustment.residuals.residual))]};
    
    for i = 1:length(hObj.Adjustment.instruments)
        stats{end+1} = ['Sigmas ' pad(hObj.Adjustment.instruments{i}, 6, 'left') ': A: ' sprintf('%.3f', hObj.Adjustment.apriori_sigmas(i)) ' P: ' sprintf('%.3f', hObj.Adjustment.aposteriori_sigmas(i)) ];
    end
    
    handles.lstStats.String = stats;
    
    % fill the outlier information listbox
    line_names = strcat(hObj.Adjustment.residuals.line_name(~hObj.Adjustment.outliers),': ');
    % get the instrument name and pad it to 6 chars
    instruments = strcat(cellfun(@(x) pad(x, 6, 'left'), hObj.Adjustment.residuals.instrument(~hObj.Adjustment.outliers),'UniformOutput',false), '->');
    % attach instrument to line name
    line_instr = strcat(line_names, instruments);
    % get the residuals
    residuals = hObj.Adjustment.residuals.residual(~hObj.Adjustment.outliers);
    % print the line as it will be shown to the user
    outlier_info = strcat(line_instr, strcat(hObj.Adjustment.residuals.delta_names(~hObj.Adjustment.outliers), sprintfc(': % .3f', residuals)));
    % sort the outlier info by descending residual
    [~,c] = sortrows(abs(residuals),'descend');
    handles.lstOutliers.String = outlier_info(c);
    
    guidata(hObject, hObj);
    
    figure;
    set(gcf, 'Name', 'Adjustment residuals');
    subplot(1,2,1)
    histfit(hObj.Adjustment.residuals.residual(hObj.Adjustment.outliers));
    grid on
    %title('(a) Outliers removed')
    title(['(a) Outliers removed. Total count: ' num2str(length(hObj.Adjustment.residuals.residual(hObj.Adjustment.outliers)))])
    xlabel('Residual [mGal]')
    ylabel('Count') 
    
    subplot(1,2,2)
    histfit(hObj.Adjustment.residuals.residual);
    grid on
    %title('(b) With outliers')
    title(['(b) With outliers. Total count: ' num2str(length(hObj.Adjustment.residuals.residual))])
    xlabel('Residual [mGal]')
    ylabel('Count') 

% --- Executes on button press in cmdOutput.
function cmdOutput_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdOutput (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    hObj = guidata(hObject);
    
    fg_adjustment_save(hObj.Adjustment);
    
    guidata(hObject, hObj);
    

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
    


% --- Executes on selection change in lstStats.
function lstStats_Callback(hObject, eventdata, handles)
% hObject    handle to lstStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstStats contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstStats


% --- Executes during object creation, after setting all properties.
function lstStats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lstOutliers.
function lstOutliers_Callback(hObject, eventdata, handles)
% hObject    handle to lstOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstOutliers contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstOutliers


% --- Executes during object creation, after setting all properties.
function lstOutliers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdView.
function cmdView_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdView (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    fg_view_solution()
    


% --- Executes on button press in cmdToggle.
function cmdToggle_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdToggle (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    vars = guidata(hObject);
        
    for i = 1:max(vars.LinesInAdjust)
        
        vars.Visualization = vars.Visualization.ToggleLineActivation(i);
    
    end
    
    set(handles.lstLines,'String', vars.Visualization.LineNameList)
    
    guidata(hObject, vars);
