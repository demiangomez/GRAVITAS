function varargout = fg_edit_line(varargin)
% FG_EDIT_LINE MATLAB code for fg_edit_line.fig
%      FG_EDIT_LINE, by itself, creates a new FG_EDIT_LINE or raises the existing
%      singleton*.
%
%      H = FG_EDIT_LINE returns the handle to a new FG_EDIT_LINE or the handle to
%      the existing singleton*.
%
%      FG_EDIT_LINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_EDIT_LINE.M with the given input arguments.
%
%      FG_EDIT_LINE('Property','Value',...) creates a new FG_EDIT_LINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_edit_line_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_edit_line_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_edit_line

% Last Modified by GUIDE v2.5 14-Sep-2024 13:36:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_edit_line_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_edit_line_OutputFcn, ...
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


% --- Executes just before fg_edit_line is made visible.
function fg_edit_line_OpeningFcn(hObject, eventdata, handles, varargin)
    
    % Center the figure on the screen
    movegui(hObject, 'center');

    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_edit_line (see VARARGIN)

    % Choose default command line output for fg_edit_line
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fg_edit_line wait for user response (see UIRESUME)
    % uiwait(handles.fg_edit_line);
    
    % load the project
    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    
    project = CssProject(vars.working_folder);
    
    database    = project.LoadDatabase();
    lines       = project.LoadLines();
    instruments = project.LoadInstruments();
    
    if ~isempty(instruments)
        set(handles.dropInstruments,'string',sort({instruments.name}));
    end

        % save the database to the local handle
    hObj                = guidata(hObject);
    hObj.database       = database;
    hObj.working_folder = vars.working_folder;
    hObj.project        = project;
    hObj.lines          = lines;
    hObj.instruments    = instruments;
    hObj.selected_line  = 0;
    hObj.selected_instrument = 0;
    hObj.selected_direction  = 0;
    hObj.selected_row        = 0;
    % create an empty gravity line and have it ready to start
    % editing/adding observations
    hObj.new_line       = CssGravityLine();
    
    hObj = Update_hObj_Sync(hObj, hObj.new_line, handles);
    
    guidata(hObject, hObj);

    set(handles.dropDirections,'String', str([CssDirections.forward; CssDirections.reverse]))
    set(handles.tblObs,'Data', [])


% --- Outputs from this function are returned to the command line.
function varargout = fg_edit_line_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in dropInstruments.
function dropInstruments_Callback(hObject, eventdata, handles)
    % hObject    handle to dropInstruments (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns dropInstruments contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from dropInstruments
    
    % when click, show the information for the current line and gravimeter
    hObj = guidata(hObject);
    
    selection = get(handles.dropInstruments,{'string','value'});
    
    hObj.selected_instrument = selection{2};
    
    % reset the selection
    set(handles.dropDirections,'value', [])
    set(handles.tblObs,'Data', [])
    hObj.selected_direction = 0;
    hObj.selected_row       = 0;
    
    line = GetLine(hObj);
    
    hObj = Update_hObj_Sync(hObj, line, handles);
    
    guidata(hObject, hObj);


function plot_line_info(hObj, handles)

    if hObj.selected_instrument ~=0
        instrument = hObj.instruments(hObj.selected_instrument);

        line = GetLine(hObj);

        if get(handles.optDeltas, 'Value') ~= 0
            if ~isnan(line.GetDeltasResiduals(instrument.name))
                axes(handles.axPlot)
                cla('reset')
                if line.GetDeltasResiduals(instrument.name)
                    % there are deltas and residuals, plot
                    plot(line, instrument.name)
                    set(handles.lblError,'Visible','off')
                else
                    axes(handles.axPlot)
                    set(handles.lblError,'Visible','on')
                end
            else
                axes(handles.axPlot)
                cla('reset')
                set(handles.lblError,'Visible','on')
            end
        else
            axes(handles.axPlot)
            set(handles.lblError,'Visible','off')
            
            if get(handles.optRaw, 'Value') ~= 0
                plotRaw(line, instrument.name)
            else
                if handles.chkBothDirs.Value
                    plotLineTime(line, instrument.name)
                else
                    if hObj.selected_direction ~= 0
                        plotLineTime(line, instrument.name, CssDirections.ToDirection(hObj.selected_direction - 1))
                    else
                        plotLineTime(line, instrument.name)
                    end
                end
            end
        end
    else
        axes(handles.axPlot)
        cla('reset')
        set(handles.lblError,'Visible','on')
    end
    
% --- Executes during object creation, after setting all properties.
function dropInstruments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropInstruments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dropDirections.
function dropDirections_Callback(hObject, eventdata, handles)
    % hObject    handle to dropDirections (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns dropDirections contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from dropDirections
    hObj = guidata(hObject);
    
    if hObj.selected_instrument ~= 0
        selection = get(handles.dropDirections,{'String','Value'});
    
        hObj.selected_direction = selection{2};
        
        line = GetLine(hObj);
        
        hObj = Update_hObj_Sync(hObj, line, handles);
        
        guidata(hObject, hObj);
    end
    

% --- Executes during object creation, after setting all properties.
function dropDirections_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropDirections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    line = GetLine(hObj);
    
    % save the line, and then create a new one
    if hObj.selected_line ~= 0
        % allow overwrite if editing a line
        line.SaveGravityLine(fullfile(hObj.working_folder,'lines'), true);
    else
        % disallow overwrite if editing a line
        line.SaveGravityLine(fullfile(hObj.working_folder,'lines'), false);
    end

    guidata(hObject, hObj);

function [hObj, line, cancel] = SaveLine(hObj, line, handles)
    
    cancel = false;
    
    if hObj.selected_line == 0
        if ~isempty(line.line_name)
            % requesting a new line, ask if user wants to save the current working
            % line
            response = questdlg('Do you want to save the current new line? If you answer No, all changes to the current new line will be lost!','Save changes?','Yes','No','Cancel','Cancel');

            if strcmp(response,'No')
                % create a new line, don't save changes
                line = CssGravityLine();
            elseif strcmp(response,'Yes')
                % save the line, and then create a new one
                if hObj.selected_line ~= 0
                    % allow overwrite if editing a line
                    line.SaveGravityLine(fullfile(hObj.working_folder,'lines'), true);
                else
                    % disallow overwrite if editing a line
                    line.SaveGravityLine(fullfile(hObj.working_folder,'lines'), false);
                end
                line = CssGravityLine();

            elseif strcmp(response,'Cancel')
                cancel = true;
            end
        end
    else
        h = msgbox('The changes performed to the current line were not saved to disk. These changes are stored in memory until you save them or discard them (by closing the edit line window).');
        waitfor(h)
        
        line = CssGravityLine();
    end
    
    if cancel == false
        hObj.selected_line = 0;
        hObj.selected_instrument = 0;
        hObj.selected_direction = 0;
        hObj.selected_row       = 0;
        UpdateCombos(hObj, line, handles);
        
        set(handles.txtLat, 'string', '')
        set(handles.txtLon, 'string', '')
        set(handles.txtHeight, 'string', '')
        set(handles.txtX, 'string', '')
        set(handles.txtY, 'string', '')
        set(handles.txtZ, 'string', '') 
        set(handles.txtComments, 'string', '') 
    end
    
    % if user cancels, then nothing happens...
    hObj = Update_hObj_Sync(hObj, line, handles);


% --- Executes on button press in cmdNew.
function cmdNew_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdNew (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    line = GetLine(hObj);
    
    [hObj, ~] = SaveLine(hObj, line, handles);
    
    guidata(hObject,hObj);
    

% --- Executes on button press in cmdOpenLine.
function cmdOpenLine_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdOpenLine (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % list the lines available:
    hObj = guidata(hObject);
    
    line = GetLine(hObj);
    
    [hObj, ~, cancel] = SaveLine(hObj, line, handles);
    
    if ~cancel
        [line, sel] = listdlg('PromptString', 'Select a line', 'SelectionMode','single','ListString',{hObj.lines.line_filename});

        if sel
            hObj.selected_line = line;
            handles.lblTitle.String = ['Gravity Line - modifying: ' hObj.lines(line).line_name];

            hObj.selected_instrument = 0;
            hObj.selected_direction  = 0;
            hObj.selected_row        = 0;

            line = GetLine(hObj);
            
            hObj = Update_hObj_Sync(hObj, line, handles);

            set(handles.txtComments,'string', line.comments)

            set(handles.txtLat, 'string', '')
            set(handles.txtLon, 'string', '')
            set(handles.txtHeight, 'string', '')
            set(handles.txtX, 'string', '')
            set(handles.txtY, 'string', '')
            set(handles.txtZ, 'string', '')
            set(handles.BenchmarkPanel, 'title', 'Benchmark Metadata')

            guidata(hObject, hObj);
        end
    end
    
function UpdateCombos(hObj, line, handles)

    for i = 1:length(hObj.instruments)
        if ismember(hObj.instruments(i).name, line.instruments)
            instruments(i) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', hObj.instruments(i).name)};
        else
            instruments(i) = {hObj.instruments(i).name};
        end
    end

    if ismember(CssDirections.forward, line.directions)
        directions(1) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', char(str(CssDirections.forward)))};
    else
        directions(1) = str(CssDirections.forward);
    end

    if ismember(CssDirections.reverse, line.directions)
        directions(2) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', char(str(CssDirections.reverse)))};
    else
        directions(2) = str(CssDirections.reverse);
    end
    
    % show deltas
    if ~hObj.selected_instrument == 0
        deltas = line.GetDeltaName(hObj.instruments(hObj.selected_instrument).name, 0);
        set(handles.lstDeltas, 'String', deltas)
        set(handles.lstDeltas, 'Value', [])
    else
        set(handles.lstDeltas, 'String', {})
        set(handles.lstDeltas, 'Value', [])
    end
    
    % reset controls section
    set(handles.dropInstruments,'String', instruments)
    set(handles.dropDirections, 'String', directions)
    
    if ~hObj.selected_instrument == 0
        set(handles.dropInstruments,'Value', hObj.selected_instrument)
    else
        set(handles.dropInstruments,'Value', [])    
    end
    
    if ~hObj.selected_direction == 0
        set(handles.dropDirections,'Value', hObj.selected_direction)
    else
        set(handles.dropDirections,'Value', [])
    end

function txtLat_Callback(hObject, eventdata, handles)
% hObject    handle to txtLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLat as text
%        str2double(get(hObject,'String')) returns contents of txtLat as a double


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
% hObject    handle to lblHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblHeight as text
%        str2double(get(hObject,'String')) returns contents of lblHeight as a double


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

% --- Executes during object creation, after setting all properties.
function lblHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblHeight (see GCBO)
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
% hObject    handle to lblX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblX as text
%        str2double(get(hObject,'String')) returns contents of lblX as a double


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


% --- Executes during object creation, after setting all properties.
function lblX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtY_Callback(hObject, eventdata, handles)
% hObject    handle to lblY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblY as text
%        str2double(get(hObject,'String')) returns contents of lblY as a double


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

% --- Executes during object creation, after setting all properties.
function lblY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZ_Callback(hObject, eventdata, handles)
% hObject    handle to lblZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblZ as text
%        str2double(get(hObject,'String')) returns contents of lblZ as a double

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

% --- Executes during object creation, after setting all properties.
function lblZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblZ (see GCBO)
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


% --- Executes on button press in cmdXYZ2LLA.
function cmdXYZ2LLA_Callback(hObject, eventdata, handles)
% hObject    handle to cmdXYZ2LLA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdOpenRINEX.
function cmdOpenRINEX_Callback(hObject, eventdata, handles)
% hObject    handle to cmdOpenRINEX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdComputeCoord.
function cmdComputeCoord_Callback(hObject, eventdata, handles)
% hObject    handle to cmdComputeCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txtComments_Callback(hObject, eventdata, handles)
    % hObject    handle to txtComments (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of txtComments as text
    %        str2double(get(hObject,'String')) returns contents of txtComments as a double
    
    % save the comments to the structure
    hObj = guidata(hObject);
    line = GetLine(hObj);
    line.comments = handles.txtComments.String;
    hObj = Update_hObj_Sync(hObj, line, handles);
    guidata(hObject,hObj);

% --- Executes during object creation, after setting all properties.
function txtComments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtComments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdDelete.
function cmdDelete_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdDelete (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    if hObj.selected_row ~= 0
        if strcmp(questdlg('Are you sure you want to delete this observation?'),'Yes')
            line = GetLine(hObj);
            
            line = line.DeleteObservation(hObj.observations(hObj.selected_row));
            
            hObj = Update_hObj_Sync(hObj, line, handles);
            guidata(hObject, hObj);
        end
    end

% --- Executes when selected cell(s) is changed in tblObs.
function tblObs_CellSelectionCallback(hObject, eventdata, handles)
    % hObject    handle to tblObs (see GCBO)
    % eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
    %	Indices: row and column indices of the cell(s) currently selecteds
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    index = eventdata.Indices;
    
    if ~isempty(index)
        row = index(1);
        col = index(2);
        hObj.selected_row = row;
    else
        row = [];
        col = [];
        hObj.selected_row = 0;
    end
    
    modifier = handles.fg_edit_line.CurrentModifier;
    if isempty(modifier)
        modifier = '';
    end
    
    if strcmp(modifier,'control')

        if and(hObj.selected_instrument ~= 0, hObj.selected_direction ~= 0)
            
            instrument = hObj.instruments(hObj.selected_instrument);

            line = GetLine(hObj);
            
            if ~isempty(row)

                % perform other operations
                switch col
                    case 1
                        % benchmark
                        if strcmp(questdlg('You are about to modify the benchmark. This operation will replace this benchmark from ALL observations (and all instruments). Do you want to continue?'),'Yes')
                            % bring up the find_benchmark dialog
                            if row > 1
                                lat = hObj.observations(row-1).benchmark.lat;
                                lon = hObj.observations(row-1).benchmark.lon;
                            elseif and(length(hObj.observations) == 1, row == 1)
                                lat = hObj.observations(row).benchmark.lat;
                                lon = hObj.observations(row).benchmark.lon;
                            elseif and(length(hObj.observations) > 1, row == 1)
                                lat = hObj.observations(row+1).benchmark.lat;
                                lon = hObj.observations(row+1).benchmark.lon;
                            else
                                lat = 0;
                                lon = 0;
                            end
                            
                            [~,benchmark] = fg_find_benchmark(lat, lon);
                            
                            if ~isempty(benchmark)
                                if and(CssBenchmark.exists(line.benchmarks, benchmark.name), length(line.directions) > 1)
                                    h = warndlg('You are trying to replace a benchmark with a benchmark that is already a member of this line. This could lead to inconsistencies during delta estimation in lines with forward and reverse directions. The change will not be performed. If you want to edit the coordinate of an existing benchmark go to Lines -> Benchmarks -> Add/Edit/Remove benchmarks.');
                                    waitfor(h);
                                else
                                    if ~isempty(benchmark)
                                        line = line.UpdateBenchmark(hObj.observations(row).benchmark, benchmark, hObj.instruments);
                                    end
                                end
                            end
                        end
                    case 2
                        % time stamp
                        d = datenum(handles.tblObs.Data{row,col},'yyyy-mm-dd HH:MM');
                        newdate = uigetdate(d);
                        if newdate
                            % find the observation and replace the timestamp
                            % this will trigger an update of all the affected
                            % variables
                            newdate = datetime(newdate,'ConvertFrom','datenum');
                            % in this section of the code it seems weird to
                            % manipulate hObj.observations(), since it looks
                            % like a different variable. However, this is a
                            % pointer to the observations in
                            % hObj.lines(selected_line) or hObj.new_line
                            
                            line = line.UpdateObservation(hObj.observations(row), 'timestamp', newdate, instrument.calibration);
                        end
                end
                
                % save the line back to the corresponding variable
                hObj = Update_hObj_Sync(hObj, line, handles);
            end
        end
    else
        if ~isempty(row)
            % when click in an existing cell, get the information of the benchmark
            set(handles.txtLat, 'string', sprintf('%.8f',hObj.observations(row).benchmark.lat))
            set(handles.txtLon, 'string', sprintf('%.8f',hObj.observations(row).benchmark.lon))
            set(handles.txtHeight, 'string', sprintf('%.3f',hObj.observations(row).benchmark.height))
            set(handles.txtX, 'string', sprintf('%.3f',hObj.observations(row).benchmark.x))
            set(handles.txtY, 'string', sprintf('%.3f',hObj.observations(row).benchmark.y))
            set(handles.txtZ, 'string', sprintf('%.3f',hObj.observations(row).benchmark.z))
            set(handles.BenchmarkPanel, 'title', ['Benchmark Metadata: ' hObj.observations(row).benchmark.name])
        end
    end
    
    if or(col == 1, col == 2)
        % clear selection that prevents from re-calling
        % edit on this cell
        jUIScrollPane = findjobj(handles.tblObs);
        jUITable = jUIScrollPane.getViewport.getView;
        jUITable.changeSelection(row-1,col-1, true, false);
    end
    % leave outside if to save the selected_row

    guidata(hObject, hObj);


% --- Executes on button press in optDeltas.
function optDeltas_Callback(hObject, eventdata, handles)
    % hObject    handle to optDeltas (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of optDeltas
    hObj = guidata(hObject);
    
    if hObj.selected_instrument ~= 0
        plot_line_info(hObj, handles)
    end


% --- Executes on button press in optRaw.
function optRaw_Callback(hObject, eventdata, handles)
    % hObject    handle to optRaw (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of optRaw

    hObj = guidata(hObject);
    
    if hObj.selected_instrument ~= 0
        plot_line_info(hObj, handles)
    end


% --- Executes on button press in optLineTimes.
function optLineTimes_Callback(hObject, eventdata, handles)
    % hObject    handle to optLineTimes (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of optLineTimes

    hObj = guidata(hObject);
    
    if hObj.selected_instrument ~= 0
        plot_line_info(hObj, handles)
    end


% --- Executes on button press in cmdInsert.
function cmdInsert_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdInsert (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % insert a blank line in the tblObs
    hObj = guidata(hObject);
    
    % call find benchmark
    
    
    if and(hObj.selected_instrument ~= 0, hObj.selected_direction ~=0)
        
        % pass a search location
        if ~isempty(hObj.observations)
            lat = hObj.observations(end).benchmark.lat;
            lon = hObj.observations(end).benchmark.lon;
        else
            lat = 0;
            lon = 0;
        end
        % open the dialog box
        [~, benchmark] = fg_find_benchmark(lat,lon);
        
        if isempty(benchmark)
            return
        end
        
        instrument = hObj.instruments(hObj.selected_instrument);
        direction  = CssDirections.ToDirection(hObj.selected_direction - 1);
        
        line = GetLine(hObj);
        
        % use the max timestamp to compute
        % the next possible datetime for the observation (to reduce user
        % errors)
        observations = line.GetObservationsByInstrumentDirection(instrument.name, direction);
        if ~isempty(observations)
            max_time = observations(end).timestamp;

            % estimate the average velocity in this line
            bench = [observations.benchmark];
            lat = [bench.lat]';
            lon = [bench.lon]';
            
            if length(observations) >= 2
                % estimate speed
                d = cumsum([0; m_idist(lon(1:end-1),lat(1:end-1),lon(2:end),lat(2:end))])/1000;
                speed = mean(diff(d)./hours(diff([observations.timestamp]'))); % km/h
            else
                % use a default speed
                speed = 25; % km/h
            end

            % distance between last observation and new benchmark
            ts = max_time + hours(m_idist(lon(end),lat(end),benchmark.lon,benchmark.lat)./1000./speed);
        else
            ts = datetime() + hours(1);
        end
    
        % if the benchmark exists in the line, get the reference to it
        % rather than a copy from fg_find_benchmark
        if CssBenchmark.exists(line.benchmarks, benchmark.name)
            benchmark = CssBenchmark.ReturnBenchmark(line.benchmarks, benchmark.name);
        end
                
        % create a new observation object. Will be deep copied by AddObservation
        new_obs = CssObservation(benchmark, direction, year(ts), month(ts), day(ts), hour(ts), minute(ts), instrument, [0,0,0]);
    
        line = line.AddObservation(new_obs);
        
        % save the line back to the corresponding variable
        hObj = Update_hObj_Sync(hObj, line, handles);
        
        guidata(hObject, hObj);
    else
        h = warndlg('Please select an instrument and a direction');
        waitfor(h);
    end

function hObj = Update_hObj_Sync(hObj, line, handles)
    % this function updates the hObj variables and syncs the display with
    % the contents of the line object
    if hObj.selected_line ~=0
        hObj.lines(hObj.selected_line) = line;
    else
        hObj.new_line = line;
        handles.lblTitle.String = ['Gravity Line - New Line: ' line.line_name];
    end
    
    if and(hObj.selected_instrument ~= 0, hObj.selected_direction ~= 0)
        instrument = hObj.instruments(hObj.selected_instrument);
        direction  = CssDirections.ToDirection(hObj.selected_direction - 1);
        
        % now filter all observations for this instrument and direction
        % using the built-in functions. This creates a shallow copy in
        % hObj.observations that points to the actual object in new_line
        observations = line.GetObservationsByInstrumentDirection(instrument.name, direction, true);

    else
        observations = [];
    end

    % sync the display
    list = ToList(line, observations);
        
    % make the combo items blue or gray, depending the contents of the
    % observations list
    UpdateCombos(hObj, line, handles)

    % copy handle to observations
    hObj.observations = observations;

    % update the table
    handles.tblObs.Data = list;
    
    % try to plot, if possible
    plot_line_info(hObj, handles)
    
% --- Executes when entered data in editable cell(s) in tblObs.
function tblObs_CellEditCallback(hObject, eventdata, handles)
    % hObject    handle to tblObs (see GCBO)
    % eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    
    if and(hObj.selected_instrument ~= 0, hObj.selected_direction ~= 0)
        % when click in an existing cell, get the information of the benchmark
        index = eventdata.Indices;
        
        if ~isempty(index)
            
            row   = index(1);
            col   = index(2);
            hObj.selected_row = row;
        
            instrument = hObj.instruments(hObj.selected_instrument);
            
            line = GetLine(hObj);
            
            switch col
                case 3
                    % offset
                    if ~all(ismember(eventdata.EditData, '0123456789+-.eEdD'))
                        h = warndlg('The entered offset value is invalid');
                        waitfor(h);
                        return
                    end
                    
                    line = line.UpdateObservation(hObj.observations(row), 'offset', str2double(eventdata.EditData), instrument.calibration);
                case {4, 5, 6}
                    % reading 1,2,3
                    if ~all(ismember(eventdata.EditData, '0123456789+-.eEdD'))
                        h = warndlg('The entered offset value is invalid');
                        waitfor(h);
                        return
                    end
                    
                    try
                        line = line.UpdateObservation(hObj.observations(row), 'reading', str2double(eventdata.EditData), instrument.calibration, col-3);
                        
                        % check if all other reading are zero. If they are,
                        % help the user by assigning them the same value
                        
                        % get the shift key modifier
                        modifier = handles.fg_edit_line.CurrentModifier;
                        if isempty(modifier)
                            modifier = '';
                        end
                        
                        for i = col-3:3
                            if or(hObj.observations(row).raw_data(i) == 0, strcmp(modifier,'shift'))
                                line = line.UpdateObservation(hObj.observations(row), 'reading', str2double(eventdata.EditData), instrument.calibration, i);
                            end
                        end
                    catch ME
                        switch ME.identifier
                            case 'CssObservation:compute_reduced_g:calibration_error'
                                h = warndlg(ME.message);
                                waitfor(h);
                            otherwise
                                rethrow(ME)
                        end
                    end
                case 9
                    % active
                    line = line.UpdateObservation(hObj.observations(row), 'active', true); % true is a dummy: the status is toggled
            end
        
            % done editing, select next row
            if length(hObj.observations) < row
                jUIScrollPane = findjobj(handles.tblObs);
                jUITable = jUIScrollPane.getViewport.getView;
                jUITable.changeSelection(row,col-1, false, false);
            end
            % save the line back to the corresponding variable
            hObj = Update_hObj_Sync(hObj, line, handles);
            
            guidata(hObject, hObj);
        end
        
    end

function line = GetLine(hObj)

    if hObj.selected_line ~=0
        % work on the selected line
        line = hObj.lines(hObj.selected_line);
    else
        % work on the new_line object
        line = hObj.new_line;
    end

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function tblObs_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to tblObs (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    


% --- Executes on button press in cmdCopy.
function cmdCopy_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdCopy (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    if and(hObj.selected_instrument ~= 0, hObj.selected_direction ~= 0)
        
        line = GetLine(hObj);
        
        direction = CssDirections.ToDirection(hObj.selected_direction - 1);
        
        benchmarks = [hObj.observations.benchmark];
        if isempty(benchmarks)
            h = warndlg('No observations to copy!');
            waitfor(h)
            return
        end
        
        % check that there are no repeated observations (same benchmark) in
        % the line. If there are, prevent from copying
        if length(unique({benchmarks.name})) ~= length(benchmarks)
            h = warndlg('This line contains a direction with more than one observation in one or more benchmarks. This operation is not allowed on this type of lines');
            waitfor(h);
            return
        end
        
        response = questdlg(['You are about to copy and flip all the observations from the ' char(str(direction)) ' line to the ' char(str(Flip(direction))) ' line. Are you sure? You will be requested to enter the datetime of the first observation in the ' char(str(Flip(direction))) ' line.'],'Copy and flip','Yes','No','No');
        
        if strcmp(response,'Yes')
            
            % request a start date
            d = datenum(max([hObj.observations.timestamp]));
            newdate = uigetdate(d);
            if newdate
                % find the observation and replace the timestamp
                % this will trigger an update of all the affected
                % variables
                newdate = datetime(newdate,'ConvertFrom','datenum');
            
                % use the delta t to move the date forward
                % 1) get the delta t with diff
                % 2) flip to get the order right
                % 3) cumsum to get the culumative time difference
                % 4) flip again to match the reverse for loop
                dt = flip(cumsum(flip([diff([hObj.observations.timestamp]'); 0])));
                
                % flip the observations
                for i = length(hObj.observations):-1:1
                    obs = copy(hObj.observations(i));
                    obs = obs.SwitchDirection();

                    obs.timestamp = newdate + dt(i);
                    obs.epoch = cal2jd(year(obs.timestamp),month(obs.timestamp),day(obs.timestamp) + hour(obs.timestamp)/24 + minute(obs.timestamp)/1440);

                    % add observation
                    line = line.AddObservation(obs);
                end

                hObj = Update_hObj_Sync(hObj, line, handles);

                guidata(hObject, hObj);
            end
        end
    end


% --- Executes on button press in cmdClone.
function cmdClone_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdClone (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    % show available instruments
    line = GetLine(hObj);

    [oinstrument, sel] = listdlg('PromptString', 'Select an instrument to clone from', 'SelectionMode','single','ListString',[line.instruments]);

    if sel
        % choose instrument to copy to
        [dinstrument, sel] = listdlg('PromptString', 'Select an instrument to clone to', 'SelectionMode','single','ListString',{hObj.instruments.name});

        if sel
            % copy all observations
            for direction = [CssDirections.forward CssDirections.reverse]
                observations = line.GetObservationsByInstrumentDirection(line.instruments{oinstrument}, direction, true);
                for i = 1:length(observations)
                    obs = copy(observations(i));
                    obs = obs.ChangeInstrument(hObj.instruments(dinstrument).name, hObj.instruments(dinstrument).calibration);
                    % add this observation (manually to avoid problems with
                    % directions in AddObservation).
                    line.observations = [line.observations; copy(obs)];
                end
            end
            % resort the line
            line = line.sort_benchmarks();
        end
        
        hObj = Update_hObj_Sync(hObj, line, handles);
        
        guidata(hObject, hObj);
    end

% --- Executes on button press in cmdPlotEdit.
function cmdPlotEdit_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPlotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in lstDeltas.
function lstDeltas_Callback(hObject, eventdata, handles)
% hObject    handle to lstDeltas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstDeltas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstDeltas


% --- Executes during object creation, after setting all properties.
function lstDeltas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstDeltas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkBothDirs.
function chkBothDirs_Callback(hObject, eventdata, handles)
    % hObject    handle to chkBothDirs (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of chkBothDirs
    hObj = guidata(hObject);
    
    if hObj.selected_instrument ~= 0
        plot_line_info(hObj, handles)
    end


% --- Executes on button press in cmdCompare.
function cmdCompare_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdCompare (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % plot all deltas on the same axes
    hObj = guidata(hObject);
    line = GetLine(hObj);
    
    if ~isempty(line.instruments)
        figure(1)
        clf
        set(gcf, 'Name', ['Instrument comparison for line ' line.line_name]);
        plotComparison(line)
    end
    


% --- Executes on button press in cmdLoad.
function cmdLoad_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdLoad (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % show open file dialog box
    [FileName, PathName] = uigetfile('*.xlsx','Select an XLSX file with the approriate format');
    
    if FileName
        
        [~, sheets] = xlsfinfo(fullfile(PathName, FileName));
        
        if mod(length(sheets),2) ~= 0
            warndlg('The number of sheets should be even: a forward and reverse for each instrument.')
        else
            
            hObj = guidata(hObject);
            
            % number of sheets should be even (forward and reverse)
            for i = 1:length(sheets)
                [obs, text] = xlsread(fullfile(PathName, FileName), sheets{i});

                benchmarks = strtrim(lower(text(:,1)));
                date = text(:,2);
                
                % identify the gravimeter
                ss = strtrim(split(sheets{i},'-'));
                if strcmp(lower(ss{2}), 'f')
                    direction = CssDirections.forward;
                else
                    direction = CssDirections.reverse;
                end
                
                instrument = hObj.instruments(ismember({hObj.instruments.name}, lower(ss{1})));
                
                if isempty(instrument)
                    warndlg(['Could not find instrument name ' lower(ss{1})])
                    return
                end
                
                for j = 1:size(obs,1)
                    % get the current line
                    line = GetLine(hObj);
                    if length(benchmarks{j}) ~= 4
                        warndlg(['Invalid benchmark name ' lower(benchmarks{j})])
                        return
                    end
                    % if the benchmark exists in the line, get the reference to it
                    % rather than a copy from fg_find_benchmark
                    if CssBenchmark.exists(line.benchmarks, benchmarks{j})
                        benchmark = CssBenchmark.ReturnBenchmark(line.benchmarks, benchmarks{j});
                    else
                        % pull the benchmark from the list of all
                        % benchmarks
                        benchmark = CssBenchmark.ReturnBenchmark(hObj.database.benchmarks, benchmarks{j});
                    end
                    
                    if isempty(benchmark)
                        % DDG Jun 7 2018: load benchmark coordinates from
                        % Excel file and add it
                        benchmark = CssBenchmark(benchmarks{j}, obs(j,1), obs(j,2), 0, obs(j,6));
                        UpdateBenchmarks(hObj.project, benchmark)
                    end
                    
                    if ~isempty(regexp(date{j},'\d+-\d+-\d+-\d+:\d+', 'once'))
                       % remove the extra dash
                       s = regexp(date{j},'(\d+)-(\d+)-(\d+)-(\d+):(\d+)','tokens');
                       date{j} = [s{1}{1} '-' s{1}{2} '-' s{1}{3} ' ' s{1}{4} ':' s{1}{5}];
                    end
                    
                    ts = datetime(date{j});

                    new_obs = CssObservation(benchmark, direction, year(ts), month(ts), day(ts), hour(ts), minute(ts), instrument, obs(j,3:5));

                    line = line.AddObservation(new_obs);

                    % save the line back to the corresponding variable
                    hObj = Update_hObj_Sync(hObj, line, handles);
                end
                
                guidata(hObject, hObj);
            end
        end
    end
    


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over cmdOpenLine.
%function cmdOpenLine_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cmdOpenLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in cmdLoadFromTablet.
function cmdLoadFromTablet_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdLoadFromTablet (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Franco S. Sobrero, Ohio State University, Sept 2024
    
    % show open file dialog box
    % [FileName, PathName] = uigetfile('*.xls','Select an XLS file with the approriate format');
    [FileName, PathName] = uigetfile({'*.xls;*.xlsx', 'Excel Files (*.xls, *.xlsx)'}, 'Select an Excel file with the appropriate format');
   
    if isequal(FileName, 0)                                                             % If the user canceled the file selection
        return;
    else
        FullFilePath = fullfile(PathName, FileName);
        [~, sheets] = xlsfinfo(FullFilePath);    
        
        if strcmp(sheets{1}, 'LINE-STRUCT')                                             % If the file comes from the tablet app Log4G, load it
            hObj = guidata(hObject);
            currentLine = GetLine(hObj);                                                % Get the data from any already-loaded lines
            selectedLines = fg_load_from_tablet(FileName, PathName, currentLine);       % Call the function
                    
            for ii = 1:length(selectedLines)
                % Check if this line name is present in an even number of sheets (forward and reverse)
                isPresent = cellfun(@(x) ~isempty(strfind(x, selectedLines{ii})), sheets);

                if mod(sum(isPresent),2) ~= 0
                    warndlg('The number of sheets should be even: a forward and reverse for each instrument.')
                else

                    for jj = 2:length(sheets)                                   % Ignore sheet 1, which is LINE-STRUCT
                        if ~isempty(strfind(sheets{jj}, selectedLines{ii}))     % If the Line Name is present in the sheet name, then load the data in this sheet
                            
                            sheet_data = readtable(FullFilePath, 'Sheet', jj, 'VariableNamingRule', 'preserve');
                            benchmarks = strtrim(lower(sheet_data.Var1));
                            date = cellstr(strrep(extractBefore(sheet_data.Var2, '.'), 'T', ' ')); % Remove the part after the '.' and replace 'T' with a space
                            
                            % if ~isnumeric(sheet_data.Var9)
                            %     sheet_data.Var9 = zeros(size(sheet_data.Var9));
                            % end

                            % Extract gravimeter name, and line Direction from sheet name (text between the second appearance of "-" and ("R" or "F")
                            % !!! IMPORTANT !!! Gravimeters' names cannot include the letters "R" of "F"
                            extract_info = regexp(sheets{jj}, '(?:[^-]*-){2}([^R|F]+)(R|F)', 'tokens');
                            
                            if ~isempty(extract_info)
                                gravime = extract_info{1,1}{1,1};
                                dir = extract_info{1,1}{1,2};
                            else
                                warndlg(['Invalid sheet name: "', sheets{jj}, '". The correct format is: LINE-{line identifier}-g{gravimeter name}{R/F}{line number (e.g., 1, 2, 3,...)}. Note: The letter "R" or "F" must be capitalized.'], 'Invalid Excel Sheet Name');
                                return
                            end
    
                            if strcmpi(dir, 'f')                            % strcmpi compares strings ignoring case
                                direction = CssDirections.forward;
                            else
                                direction = CssDirections.reverse;
                            end
                            
                            instrument = hObj.instruments(ismember({hObj.instruments.name}, lower(gravime)));
                            
                            if isempty(instrument)
                                warndlg(['Could not find instrument name ' lower(gravime)])
                                return
                            end
                            
                            for j = 1:size(sheet_data,1)                    % Loop through the observations
                                line = GetLine(hObj);                       % Get the current line
                                if length(benchmarks{j}) ~= 4
                                    warndlg(['Invalid benchmark name ' lower(benchmarks{j})])
                                    return
                                end
                            
                                % In older versions of Log4G, the app had a bug which exported the Offsets (Var.8) as '????'. In this case, we replace it by 0
                                if isnumeric(sheet_data.Var8)
                                    offs = sheet_data.Var8(j);
                                else 
                                    offs = 0;
                                end

                                % If the benchmark exists in the line, get the reference to it rather than a copy from fg_find_benchmark
                                if CssBenchmark.exists(line.benchmarks, benchmarks{j})
                                    benchmark = CssBenchmark.ReturnBenchmark(line.benchmarks, benchmarks{j});
                                else
                                    % Pull the benchmark from the list of all benchmarks
                                    benchmark = CssBenchmark.ReturnBenchmark(hObj.database.benchmarks, benchmarks{j});
                                end
                               
                                if isempty(benchmark)
                                    % DDG Jun 7 2018: load benchmark coordinates from Excel file and add it
                                    benchmark = CssBenchmark(benchmarks{j}, sheet_data.Var3(j), sheet_data.Var4(j), 0, offs);
                                    UpdateBenchmarks(hObj.project, benchmark)
                                end
                                
                                ts = datetime(date{j});
            
                                new_obs = CssObservation(benchmark, direction, year(ts), month(ts), day(ts), hour(ts), minute(ts), instrument, sheet_data{j, 5:7});
            
                                line = line.AddObservation(new_obs);
            
                                % Save the line back to the corresponding variable
                                hObj = Update_hObj_Sync(hObj, line, handles);
                            end     
                            guidata(hObject, hObj);
                        end
                    end
                end
            end
        else
            % If the first sheet is not LINE-STRUCT, display an error message in a separate window
            errordlg('You must select an XLS file exported by the app Log4G. The first sheet must be "LINE-STRUCT"', 'Invalid File', 'modal');
            % Abort the function
            return;
        end
    end

