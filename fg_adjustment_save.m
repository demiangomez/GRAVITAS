function varargout = fg_adjustment_save(varargin)
% FG_ADJUSTMENT_SAVE MATLAB code for fg_adjustment_save.fig
%      FG_ADJUSTMENT_SAVE, by itself, creates a new FG_ADJUSTMENT_SAVE or raises the existing
%      singleton*.
%
%      H = FG_ADJUSTMENT_SAVE returns the handle to a new FG_ADJUSTMENT_SAVE or the handle to
%      the existing singleton*.
%
%      FG_ADJUSTMENT_SAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_ADJUSTMENT_SAVE.M with the given input arguments.
%
%      FG_ADJUSTMENT_SAVE('Property','Value',...) creates a new FG_ADJUSTMENT_SAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_adjustment_save_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_adjustment_save_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_adjustment_save

% Last Modified by GUIDE v2.5 23-Jan-2018 18:42:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_adjustment_save_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_adjustment_save_OutputFcn, ...
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
    
% --- Executes just before fg_adjustment_save is made visible.
function fg_adjustment_save_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_adjustment_save (see VARARGIN)

    % Choose default command line output for fg_adjustment_save
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    hObj = guidata(hObject);
    
    % get the adjustment object passed by the parent
    hObj.Adjustment = varargin{1};
    
    guidata(hObject, hObj);
    
    % UIWAIT makes fg_adjustment_save wait for user response (see UIRESUME)
    uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fg_adjustment_save_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    % varargout{1} = handles.output;
    delete(handles.figure1)



function txtPath_Callback(hObject, eventdata, handles)
% hObject    handle to txtPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPath as text
%        str2double(get(hObject,'String')) returns contents of txtPath as a double


% --- Executes during object creation, after setting all properties.
function txtPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdSelect.
function cmdSelect_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSelect (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % output the result to a text file
    [file, path] = uiputfile('*.txt','Save Gravity Results');
    
    if ~isempty(file)
        handles.txtPath.String = fullfile(path,file);
    end

% --- Executes on button press in chkZeroHRemove.
function chkZeroHRemove_Callback(hObject, eventdata, handles)
% hObject    handle to chkZeroHRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkZeroHRemove


% --- Executes on button press in chkNoGRemove.
function chkNoGRemove_Callback(hObject, eventdata, handles)
% hObject    handle to chkNoGRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkNoGRemove


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hObj = guidata(hObject);
    
    if ~isempty(handles.txtPath.String)
        
        handles.figure1.Pointer = 'watch';
        
        hObj.Adjustment = hObj.Adjustment.GetAnomalies();
            
        fid = fopen(handles.txtPath.String,'w');
        fprintf(fid,'Number Station Latitude_[deg] Longitude_[deg] Height-WGS84_[m] Gravity_[mGal] Uncertainty_[mGal] Free-Air_Anomaly_[mGal] Bouguer_Anomaly_[mGal] Undulation_[m]\n');
        
        out(:,1) = num2cell((1:length(hObj.Adjustment.benchmarks))');
        out(:,2) = {hObj.Adjustment.benchmarks.name}';
        out(:,3) = {hObj.Adjustment.benchmarks.lat}';
        out(:,4) = {hObj.Adjustment.benchmarks.lon}';
        out(:,5) = {hObj.Adjustment.benchmarks.height}';
        out(:,6) = num2cell(hObj.Adjustment.adjusted_g);
        out(:,7) = num2cell(hObj.Adjustment.adjusted_g_sigma);
        out(:,8) = num2cell(hObj.Adjustment.fa_anomalies);
        out(:,9) = num2cell(hObj.Adjustment.ba_anomalies);
        out(:,10) = num2cell(hObj.Adjustment.ondulations);

        % sort by station name
        if handles.optAlphabetically.Value
            [~,c] = sortrows(out,2);
        else
            [~,c] = sortrows(out,7,'descend');
        end
        
        % sort according to selected output
        out = out(c,:);
        
        % remove any NaN stations (no adjusted gravity)
        if handles.chkNoGRemove.Value
            out = out(~cellfun(@isnan, out(:,6)),:);
        end
        
        % remove stations that have ~ zero height
        if handles.chkZeroHRemove.Value
            out = out(round(cell2mat(out(:,5))) ~= 0,:);
        end
        
        % remove points with uncert > 1 mGal
        out = out(cell2mat(out(:,7)) < 1,:);
        
        % renumber the rows
        out(:,1) = num2cell((1:size(out,1))');
        out = out';
        
        fprintf(fid,'%6d %7s %14.8f %15.8f %16.3f %14.3f %18.3f %23.3f %22.3f %14.3f\n', out{:});
        
        fclose(fid);
        
        handles.figure1.Pointer = 'arrow';
        
        uiresume(handles.figure1);
    else
        h = warndlg('Select a file before saving!');
        waitfor(h);
    end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton3 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    hObj = guidata(hObject);
    
    if ~isempty(handles.txtPath.String)
        
        icons = repmat({'http://maps.google.com/mapfiles/kml/shapes/square.png'}, length(hObj.Adjustment.benchmarks), 1);
        icons(isnan([hObj.Adjustment.adjusted_g])) = repmat({'http://maps.google.com/mapfiles/kml/pal4/icon48.png'}, sum(isnan([hObj.Adjustment.adjusted_g])), 1);
        
        kmlwritepoint([handles.txtPath.String '.kml'], [hObj.Adjustment.benchmarks.lat]', ...
            [hObj.Adjustment.benchmarks.lon]', [hObj.Adjustment.benchmarks.height]', ...
            'name', {hObj.Adjustment.benchmarks.name}',...
            'icon', icons, ...
            'IconScale', 0.4, ...
            'description', cellstr(num2str([[hObj.Adjustment.benchmarks.lat]', [hObj.Adjustment.benchmarks.lon]', [hObj.Adjustment.benchmarks.height]'])));
    end
    