function varargout = fg_instruments(varargin)
% FG_INSTRUMENTS MATLAB code for fg_instruments.fig
%      FG_INSTRUMENTS, by itself, creates a new FG_INSTRUMENTS or raises the existing
%      singleton*.
%
%      H = FG_INSTRUMENTS returns the handle to a new FG_INSTRUMENTS or the handle to
%      the existing singleton*.
%
%      FG_INSTRUMENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_INSTRUMENTS.M with the given input arguments.
%
%      FG_INSTRUMENTS('Property','Value',...) creates a new FG_INSTRUMENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_instruments_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_instruments_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_instruments

% Last Modified by GUIDE v2.5 13-Sep-2017 13:03:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_instruments_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_instruments_OutputFcn, ...
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


% --- Executes just before fg_instruments is made visible.
function fg_instruments_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_instruments (see VARARGIN)

    % Choose default command line output for fg_instruments
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fg_instruments wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    h = findobj('Tag','fg_main');
    
    vars = guidata(h);
    
    project = CssProject(vars.working_folder);
    
    instruments = project.LoadInstruments();
    
    if ~isempty(instruments)
        set(handles.dropInstruments,'String',{instruments.name})
    end
    
    hObj                = guidata(hObject);
    hObj.database       = project.LoadDatabase();
    hObj.working_folder = vars.working_folder;
    
    guidata(hObject, hObj);

% --- Outputs from this function are returned to the command line.
function varargout = fg_instruments_OutputFcn(hObject, eventdata, handles) 
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
    
    i = get(hObject,'Value');
    
    % load the instrument data
    hObj = guidata(hObject);

    project = CssProject(hObj.working_folder);
    instruments = project.LoadInstruments();
    
    % reload instruments
    load(fullfile(hObj.working_folder,['instruments/' instruments(i).name '.mat']));
    
    set(handles.txtCalibration, 'String', sprintf('%8.3f,\n',instrument_struct.calibration))


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


function txtCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to txtCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCalibration as text
%        str2double(get(hObject,'String')) returns contents of txtCalibration as a double


% --- Executes during object creation, after setting all properties.
function txtCalibration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCalibration (see GCBO)
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

    % for the current instrument selected, save the data to the file
    % remove any non-valid chars
    txt = reshape(get(handles.txtCalibration,'String') .',1,[]);
    txt = strrep(txt,',',' ');
    txt = strrep(txt,'...','');
    
    v = char(strsplit(txt));
    
    i = get(handles.dropInstruments,'Value');
    
    % load the instrument data
    hObj = guidata(hObject);

    project = CssProject(hObj.working_folder);
    instruments = project.LoadInstruments();
    
    instrument_struct = struct('name', instruments(i).name, 'calibration', str2num(v)); % do not replace with str2double!
        
    save(fullfile(hObj.working_folder,['instruments/' instruments(i).name '.mat']),'instrument_struct');
    

% --- Executes on button press in cmdNew.
function cmdNew_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdNew (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % open a dialog box to enter the name of the new instrument

    name = inputdlg('Enter the instrument name:', 'New Instrument', [1 50]);
    
    if ~isempty(name)
        hObj = guidata(hObject);
        
        % create the instrument file
        instrument_struct = struct('name', name{1}, 'calibration', 0);
        
        save(fullfile(hObj.working_folder,['instruments/' name{1} '.mat']),'instrument_struct');
        
        % reload instruments
        project = CssProject(hObj.working_folder);
        
        instruments = project.LoadInstruments();
    
        if ~isempty(instruments)
            set(handles.dropInstruments,'String',{instruments.name})
        end
    end
    
