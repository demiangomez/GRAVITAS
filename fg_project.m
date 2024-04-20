function varargout = fg_project(varargin)
% FG_PROJECT MATLAB code for fg_project.fig
%      FG_PROJECT, by itself, creates a new FG_PROJECT or raises the existing
%      singleton*.
%
%      H = FG_PROJECT returns the handle to a new FG_PROJECT or the handle to
%      the existing singleton*.
%
%      FG_PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_PROJECT.M with the given input arguments.
%
%      FG_PROJECT('Property','Value',...) creates a new FG_PROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_project_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_project_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_project

% Last Modified by GUIDE v2.5 24-Aug-2017 13:44:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_project_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_project_OutputFcn, ...
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

% --- Executes just before fg_project is made visible.
function fg_project_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fg_project (see VARARGIN)

% Choose default command line output for fg_project
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fg_project wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fg_project_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    vars = guidata(hObject);
    
    % array to save missing GPS solutions
    missing_soln = [];
    benchmarks   = [];
    
    % save all the CssGravityClasses
    for i = 1:length(vars.cssLines)
        
        % find the points in the benchmarks array
        if isfield(vars,'benchmarks')
            % the solutions var exists
            for j = 1:length(vars.cssLines(i).benchmarks)
                % find the corresponding benchmark in the GPS solutions

                [~,index] = ismember(vars.cssLines(i).benchmarks(j).name, [vars.benchmarks.name]');
                
                % if it wasn't found using the regular name, try to old
                % name
                if index == 0
                    [~,index] = ismember(CssGravityLine.rename_back(vars.cssLines(i).benchmarks(j).name), [vars.benchmarks.name]');
                end
                
                if index ~= 0
                    vars.cssLines(i).benchmarks(j) = vars.cssLines(i).benchmarks(j).UpdateCoords(vars.benchmarks(index));
                else
                    missing_soln = [missing_soln; {vars.cssLines(i).benchmarks(j).name}];
                end
            end
            
        end
        
        % save benchmarks to benchmark names variable
        for j = 1:length(vars.cssLines(i).benchmarks)
            if ~CssBenchmark.exists(benchmarks, vars.cssLines(i).benchmarks(j).name)
                benchmarks = [benchmarks; vars.cssLines(i).benchmarks(j)];
            end
        end
        
        vars.cssLines(i).SaveGravityLine(fullfile(vars.working_folder,'lines'));
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Absolute gravity
    
    for i = 1:length(vars.cssLinesAGrav)
        
        % find the points in the benchmarks array
        if isfield(vars,'benchmarks')
            % the solutions var exists
            for j = 1:length(vars.cssLinesAGrav(i).benchmarks)
                % find the corresponding benchmark in the GPS solutions

                [~,index] = ismember(vars.cssLinesAGrav(i).benchmarks(j).name, [vars.benchmarks.name]');
                
                if index ~= 0
                    vars.cssLinesAGrav(i).benchmarks(j) = vars.cssLinesAGrav(i).benchmarks(j).UpdateCoords(vars.benchmarks(index));
                else
                    missing_soln = [missing_soln; {vars.cssLinesAGrav(i).benchmarks(j).name}];
                end
            end
            
        end
        
        % save benchmarks to benchmark names variable
        for j = 1:length(vars.cssLinesAGrav(i).benchmarks)
            if ~CssBenchmark.exists(benchmarks, vars.cssLinesAGrav(i).benchmarks(j).name)
                benchmarks = [benchmarks; vars.cssLinesAGrav(i).benchmarks(j)];
            end
        end
        
        % after a discussion with Kevin, absolute gravity lines can be
        % treated as regular lines
        vars.cssLinesAGrav(i).SaveGravityLine(fullfile(vars.working_folder,'lines'));
    end
    
    database = struct();
    database.benchmarks = benchmarks;
    database.agravbench = [];
    database.proj_description = get(handles.txtProjectDesc,'String');
    save(fullfile(vars.working_folder, 'database.mat'),'database');
    
    if ~isempty(missing_soln)
        msgbox(['Some benchmarks were not found in the GPS solutions directory. List follows: ' strjoin(unique(missing_soln))])
    end
    

function txtProjectDesc_Callback(hObject, eventdata, handles)
% hObject    handle to txtProjectDesc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtProjectDesc as text
%        str2double(get(hObject,'String')) returns contents of txtProjectDesc as a double


% --- Executes during object creation, after setting all properties.
function txtProjectDesc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtProjectDesc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdDiscard.
function cmdDiscard_Callback(hObject, eventdata, handles)
% hObject    handle to cmdDiscard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdImport.
function cmdImport_Callback(hObject, eventdata, handles)
% hObject    handle to cmdImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdImportRelative.
function cmdImportRelative_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdImportRelative (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % show a dialog box to request the user folder to import data from
    relative_folder = uigetdir('~','Select the folder with the relative gravity lines:');
    
    if relative_folder == 0
        return
    end
        
    % loop through the folder files searching for lines
    files = dir(relative_folder);
    % remove the directories
    files([files.isdir] == 1) = [];

    j = 1;
    h = waitbar(0,'Analyzing lines...');
    for i = 1:length(files)
        % split the fields using the delimiter
        fields = strsplit(files(i).name,'_');

        if length(fields) <= 3
            % to avoid fu$%@ DS_STORE files in the fu$%@ Macs
            continue
        end
        end_points = sortrows([fields(2);fields(3)]);
        
        lines(j,:) = [pad(end_points{1},15) pad(end_points{2},15)];
        j = j + 1;
        
        waitbar(i/length(files))
        
    end
    close(h);
    
    % get the unique pairs
    lines = unique(lines,'rows');
    
    h = waitbar(0,'Importing lines...');
    for i = 1:size(lines,1)
        % load each of the lines
        cssLines(i) = CssGravityLine(relative_folder,strtrim(lines(i,1:15)),strtrim(lines(i,16:30)));
        waitbar(i/length(lines))
    end
    close(h);
    
    % make a list of the cities
    benchmarks = cellstr(unique([lines(:,1:15);lines(:,16:30)],'rows'));
    
    hObj = guidata(hObject);
    hObj.cssLines = cssLines;
    guidata(hObject, hObj);
    
    set(handles.lblTotalCities1,'String',['Total imported lines: ' num2str(length(benchmarks))]);
    
    msgbox('Successful import of relative gravity lines!')
    

% --- Executes on button press in cmdImportSecondOrder.
function cmdImportSecondOrder_Callback(hObject, eventdata, handles)
% hObject    handle to cmdImportSecondOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdImportAbsolute.
function cmdImportAbsolute_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdImportAbsolute (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % show a dialog box to request the user folder to import data from
    absolute_folder = uigetdir('~','Select the folder with the absolute gravity links:');
    
    if absolute_folder == 0
        return
    end
        
    % loop through the folder files searching for lines
    files = dir(absolute_folder);
    % remove the directories
    files([files.isdir] == 1) = [];

    j = 1;
    h = waitbar(0,'Analyzing lines...');
    for i = 1:length(files)
        % split the fields using the delimiter
        fields = strsplit(files(i).name,'_');

        if or(length(fields) <= 3, strcmp(fields{1},'.'))
            % to avoid fu$%@ DS_STORE files in the fu$%@ Macs
            continue
        end
        end_points = sortrows([fields(2);fields(3)]);
        
        lines(j,:) = [pad(end_points{1},15) pad(end_points{2},15)];
        j = j + 1;
        
        waitbar(i/length(files))
        
    end
    close(h);
    
    % get the unique pairs
    lines = unique(lines,'rows');
    
    h = waitbar(0,'Importing lines...');
    for i = 1:size(lines,1)
        % load each of the lines
        cssLinesAGrav(i) = CssGravityLine(absolute_folder,strtrim(lines(i,1:15)),strtrim(lines(i,16:30)));
        waitbar(i/length(lines))
    end
    close(h);
    
    % make a list of the absolute benchmarks
    benchmarks = cellstr(unique([lines(:,1:15);lines(:,16:30)],'rows'));
    
    hObj = guidata(hObject);
    hObj.cssLinesAGrav = cssLinesAGrav;
    guidata(hObject, hObj);
    
    set(handles.lblTotalCities3,'String',['Total imported links: ' num2str(length(benchmarks))]);
    
    msgbox('Successful import of absolute gravity links!')

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdGPSSoln.
function cmdGPSSoln_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdGPSSoln (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % show a dialog box to request the user folder to import data from
    gps_folder = uigetdir('~','Select the folder with the GPS solutions:');
    
    if gps_folder == 0
        return
    end
    
    set(handles.lblGpsPath, 'String', gps_folder);
    
    % loop through the folder files searching for lines
    files = dir(gps_folder);
    % remove the directories
    files([files.isdir] == 1) = [];

    j = 1;
    h = waitbar(0,'Loading GPS solutions...');
    for i = 1:length(files)
        % split the fields using the delimiter
        fields = strsplit(files(i).name,'_');

        if or(strcmp(fields(1), '.DS'), strcmp(fields(1), '.'))
            % to avoid fu$%@ DS_STORE files in the fu$%@ Macs
            continue
        end
        
        % open the file
        loc = load(fullfile(gps_folder,files(i).name));
        
        benchmarks(j) = CssBenchmark(fields(1), loc(1), loc(2), loc(3), 0);
        
        j = j + 1;
        
        waitbar(i/length(files))
        
    end
    close(h);
    
    hObj = guidata(hObject);
    hObj.benchmarks = benchmarks;
    guidata(hObject, hObj);

    msgbox('Successful import of GPS solutions!')
    
