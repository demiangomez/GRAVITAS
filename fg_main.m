function varargout = fg_main(varargin)
% FG_MAIN MATLAB code for fg_main.fig
%      FG_MAIN, by itself, creates a new FG_MAIN or raises the existing
%      singleton*.
%
%      H = FG_MAIN returns the handle to a new FG_MAIN or the handle to
%      the existing singleton*.
%
%      FG_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_MAIN.M with the given input arguments.
%
%      FG_MAIN('Property','Value',...) creates a new FG_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_main

% Last Modified by GUIDE v2.5 29-Sep-2020 18:49:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_main_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_main_OutputFcn, ...
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


% --- Executes just before fg_main is made visible.
function fg_main_OpeningFcn(hObject, eventdata, handles, varargin)

    % Center the figure on the top of the screen
    movegui(hObject, 'center');

    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_main (see VARARGIN)

    % Choose default command line output for fg_main
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes fg_main wait for user response (see UIRESUME)
    % uiwait(handles.fg_main);

    % disable all menus until user sets working path

    set(handles.mnuIntruments,'Enable', 'off');
    set(handles.mnuLines,'Enable', 'off');
    set(handles.mnuGravity,'Enable', 'off');
    set(handles.mnuCreateFolderStruct,'Enable','off');
    set(handles.mnuCheckFolderStruct,'Enable','off');


% --- Outputs from this function are returned to the command line.
function varargout = fg_main_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


% --------------------------------------------------------------------
function mnuIntruments_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuIntruments (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuAddInstrument_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuAddInstrument (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hfg_instruments = fg_instruments();
    

% --------------------------------------------------------------------
function mnuInstStatistics_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuInstStatistics (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % bring up a couple statistic windows for a particular instrument
    hObj = guidata(hObject);
    project = CssProject(hObj.working_folder);
    instruments = project.LoadInstruments();
    
    [inst, sel] = listdlg('PromptString', 'Select an instrument', 'SelectionMode','single','ListString',{instruments.name});
    
    if sel
        % an instrument was selected
        lines = project.LoadLines();
        
        residuals = [];
        % make a histogram of the residuals of this instrument
        for i = 1:length(lines)
            [~,r] = lines(i).GetDeltasResiduals(instruments(inst).name);
            residuals = [residuals; r];
        end
        
        figure('Name', ['Delta residuals for instrument ' instruments(inst).name]);
        histfit(residuals,100);
        grid on
        title(['Standard deviation: ' sprintf('%5.2f', nanstd(residuals)) ' mGal'])
        xlabel('Residual [mGal]')
        ylabel('Count')
    end

% --------------------------------------------------------------------
function mnuFolders_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuFolders (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuLines_Callback(hObject, eventdata, handles)
% hObject    handle to mnuLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuGravity_Callback(hObject, eventdata, handles)
% hObject    handle to mnuGravity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuAbsoluteGravity_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuAbsoluteGravity (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hfg_absolute = fg_absolute();
    
    
% --------------------------------------------------------------------
function mnuAdjustment_Callback(hObject, eventdata, handles)
% hObject    handle to mnuAdjustment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuBenchmarks_Callback(hObject, eventdata, handles)
% hObject    handle to mnuBenchmarks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuGravityLine_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuGravityLine (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    hfg_edit_line = fg_edit_line();

% --------------------------------------------------------------------
function mnuEditGravityLine_Callback(hObject, eventdata, handles)
% hObject    handle to mnuEditGravityLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuRemoveGravityLine_Callback(hObject, eventdata, handles)
% hObject    handle to mnuRemoveGravityLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuSetWorkingFolder_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuSetWorkingFolder (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    hObj.working_folder = uigetdir('~','Select the working folder:');
    hObj.valid_folder_struct = -1;
    
    if hObj.working_folder ~= 0
        
        % trigger the check structure functions
        result = check_folder_struct(hObj.working_folder, handles);

        if result ~= 0
            choice = questdlg('The folder structure of a gravity project was not detected. Would you like to create it?', ...
                'Folder structure check not passed!', ...
                'Yes','No','No');
            
            if choice == 'Yes'
                create_folder_structure(hObj.working_folder)
                hObj.valid_folder_struct = 0;
            end
        else
            hObj.valid_folder_struct = 0;
        end
        
        if hObj.valid_folder_struct == 0
            set(handles.mnuIntruments,'Enable', 'on');
            set(handles.mnuLines,'Enable', 'on');
            set(handles.mnuGravity,'Enable', 'on');
            set(handles.mnuCreateFolderStruct,'Enable','off');
            set(handles.mnuCheckFolderStruct,'Enable','on');
            
            % load the database
            load(fullfile(hObj.working_folder, 'database.mat'))
            
            set(handles.lblFolder, 'String', [hObj.working_folder ' (' database.proj_description ')' ]);
        else
            set(handles.mnuCreateFolderStruct,'Enable','on');
        end
        
        guidata(hObject, hObj);
    
    end
    
    

% --------------------------------------------------------------------
function mnuCreateFolderStruct_Callback(hObject, eventdata, handles)
% hObject    handle to mnuCreateFolderStruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuCheckFolderStruct_Callback(hObject, eventdata, handles)
% hObject    handle to mnuCheckFolderStruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuAddBenchmark_Callback(hObject, eventdata, handles)
% hObject    handle to mnuAddBenchmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuAddEditRemoveBenchmark_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuAddEditRemoveBenchmark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % call fg_find_benchmark
    fg_find_benchmark();

% --------------------------------------------------------------------
function mnuAdjust_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuAdjust (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    fg_adjust();
    

% --------------------------------------------------------------------
function mnuStats_Callback(hObject, eventdata, handles)
% hObject    handle to mnuStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section for user defined functions

function value = check_folder_struct(working_folder, handles)
    % this function check the existance of all the necessary folders and
    % structures.
    
    value = -1;

    if ~exist(fullfile(working_folder, 'lines'),'dir')
        set(handles.lblFolder, 'String', 'checking lines folder...');
        return
    end
    
    if ~exist(fullfile(working_folder, 'instruments'),'dir')
        set(handles.lblFolder, 'String', 'checking instruments...');
        return
    end
    
%     if ~exist(fullfile(working_folder, 'absolute'),'dir')
%         set(handles.lblFolder, 'String', 'checking absolute gravity data...');
%         return
%     end
    
    if ~exist(fullfile(working_folder, 'database.mat'),'file')
        set(handles.lblFolder, 'String', 'checking database mat file...');
        return
    end
    
    value = 0;


function create_folder_structure(working_folder)
    % this function creates all the necessary folders and
    % structures.
    
    if ~exist([working_folder '/lines'],'dir')
        mkdir([working_folder '/lines'])
    end
    
    if ~exist([working_folder '/instruments'],'dir')
        mkdir([working_folder '/instruments'])
    end
    
%     if ~exist([working_folder '/absolute'],'dir')
%         mkdir([working_folder '/absolute'])
%     end
    
    if ~exist([working_folder '/database.mat'],'file')
        % no basic information in the folder! bring up the config window
        hfg_project = fg_project();
        
        hfg_project_data = guidata(hfg_project);
        hfg_project_data.working_folder = working_folder;
        guidata(hfg_project, hfg_project_data);
        waitfor(hfg_project);
    end
        
% end section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------
function cmdLineStats_Callback(hObject, eventdata, handles)
    % hObject    handle to cmdLineStats (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    hfg_linestats = fg_linestats();
    


% --------------------------------------------------------------------
function mnuRejected_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuRejected (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % output a KML file with all 
    vars = guidata(hObject);
    project = CssProject(vars.working_folder);
    lines  = project.LoadLines();
    
    database = project.LoadDatabase();
    benchmarks = database.benchmarks; 
    
    network = CssNetwork(lines, benchmarks);
    
    k = kml('Gravity Lines') ; 
    colors = {'FF0000FF', 'FF00FFFF', 'FF00FF00'};
    icons = {'http://maps.google.com/mapfiles/kml/shapes/donut.png', 'http://maps.google.com/mapfiles/kml/shapes/square.png'};
    
    f = k.createFolder('Rejected benchmarks/lines');
    r = {};
    g = {};
    for i = 1:length(lines)
        if any(sum(lines(i).status) < 2) % the purpose of this if sentence is to check benchmark redundancy
            % if less than 2 positive status (i.e. positive = 1) then no redundancy.
            % this line has to go into the KML
            
            % get the number of observations available
            cc = sum(lines(i).status);
            % if more than 3, then assign 3
            cc(cc > 2) = 2;
            % convert to index from 1-4
            cc = cc + 1;
            
            folder = f.createFolder(lines(i).line_name);
            
            % circle or square depending on node or point
            nodes = ismember({lines(i).benchmarks.name}', {network.nodes.name});
            nodes = nodes + 1;
            
            folder.plot3([lines(i).benchmarks.lon], [lines(i).benchmarks.lat], [lines(i).benchmarks.height], 'lineWidth', 3, 'name', 'gravity line', 'lineColor', colors{min(cc)}, 'altitudeMode', 'clampToGround');
            folder.point([lines(i).benchmarks.lon], [lines(i).benchmarks.lat], [lines(i).benchmarks.height], ...
                'iconScale', 0.5, ...
                'name', {lines(i).benchmarks.name}', ...
                'iconURL', icons(nodes), ...
                'iconColor', colors(cc)', 'altitudeMode', 'clampToGround');
            r{i,1} = {[lines(i).benchmarks.lon]', [lines(i).benchmarks.lat]', [lines(i).benchmarks.height]', {lines(i).benchmarks.name}', colors(cc)'};
        end
    end
    
    f = k.createFolder('Good benchmarks/lines');

    for i = 1:length(lines)
        if all(sum(lines(i).status) >= 2)
            
            % get the number of observations available
            cc = sum(lines(i).status);
            % if more than 3, then assign 3
            cc(cc > 2) = 2;
            % convert to index from 1-4
            cc = cc + 1;
            
            folder = f.createFolder(lines(i).line_name);
            
            % circle or square depending on node or point
            nodes = ismember({lines(i).benchmarks.name}', {network.nodes.name});
            nodes = nodes + 1;
            
            folder.plot3([lines(i).benchmarks.lon], [lines(i).benchmarks.lat], [lines(i).benchmarks.height], 'lineWidth', 3, 'name', 'gravity line', 'lineColor', 'FFFF0000', 'altitudeMode', 'clampToGround');
            folder.point([lines(i).benchmarks.lon], [lines(i).benchmarks.lat], [lines(i).benchmarks.height], ...
                'iconScale', 0.5, ...
                'name', {lines(i).benchmarks.name}', ...
                'iconURL', icons(nodes), ...
                'altitudeMode', 'clampToGround');
        end
        g{i,1} = {[lines(i).benchmarks.lon]', [lines(i).benchmarks.lat]', [lines(i).benchmarks.height]', {lines(i).benchmarks.name}', colors(cc)'};
    end
    g = g(~cellfun('isempty',g));
    r = r(~cellfun('isempty',r));
    k.save('test.kml')
    save('good.mat', 'g')
    save('bad.mat', 'r')

    
    


% --------------------------------------------------------------------
function mnuInstDrifts_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuInstDrifts (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % bring up a couple statistic windows for a particular instrument
    hObj = guidata(hObject);
    project = CssProject(hObj.working_folder);
    instruments = project.LoadInstruments();
    
    [inst, sel] = listdlg('PromptString', 'Select an instrument', 'SelectionMode','single','ListString',{instruments.name});
    
    if sel
        % an instrument was selected
        lines = project.LoadLines();
        
        drifts = [];
        % make a histogram of the residuals of this instrument
        for i = 1:length(lines)
            d = lines(i).GetDrifts(instruments(inst).name);
            drifts = [drifts; d];
        end
        
        figure('Name', ['Drifts for instrument ' instruments(inst).name]);
        histfit(drifts,100);
        grid on
        title(['Standard deviation: ' sprintf('%5.2f', nanstd(drifts)) ' mGal/h'])
        xlabel('[mGal/h]')
        ylabel('Count')
    end


% --------------------------------------------------------------------
function mnuExportDatabase_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuExportDatabase (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % this menu creates a database file (Excel) with all the instruments
    % and points to syncronize the tablet software
    
    hObj = guidata(hObject);
    
    [file, path] = uiputfile('*.xls','Export database for tablet');
    
    if ~isempty(file)
        filename = fullfile(path,file);
        
        project = CssProject(hObj.working_folder);
        
        database = project.LoadDatabase();
        benchmarks = database.benchmarks; 
        
        t = table({benchmarks.name}', [benchmarks.lat]', [benchmarks.lon]',...
                  [benchmarks.height]', [benchmarks.offset]', [benchmarks.x]',...
                  [benchmarks.y]', [benchmarks.z]', ...
                  zeros(size([benchmarks.z]')), zeros(size([benchmarks.z]')), ...
                  'VariableNames', {'name', 'lat', 'lon', 'height', 'offset', ...
                  'x', 'y', 'z', 'absolute_g', 'uncertainty'});
        
        writetable(t, filename, 'Sheet', 'Sheet1')
        
        instruments = project.LoadInstruments();
        
        for i = 1:size(instruments, 1)
            writetable(table(instruments(i).calibration), filename, 'Sheet', instruments(i).name, 'WriteVariableNames', false)
        end
        
        msgbox(['Successfully exported the database to ' filename]);
    end


% --------------------------------------------------------------------
function mnuZeroHeight_Callback(hObject, eventdata, handles)
    % hObject    handle to mnuZeroHeight (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    hObj = guidata(hObject);
    
    [file, path] = uiputfile('*.kmz','Output KML with benchmarks with zero height');
    
    if ~isempty(file)
        k = kml('Benchmarks') ; 
        project = CssProject(hObj.working_folder);
        icons =  {'http://maps.google.com/mapfiles/kml/shapes/donut.png'};
        colors = {'FF0000FF', 'FF00FFFF', 'FF00FF00'};
        
        database = project.LoadDatabase();
        benchmarks = database.benchmarks([database.benchmarks.height] == 0);
        
        c = ones(size(benchmarks,1),1);
        
        for i = 1:size(benchmarks,1)
            lon(i) = benchmarks(i).lon;
            lat(i) = benchmarks(i).lat;
            h(i) = benchmarks(i).height;
            name{i} = benchmarks(i).name;
            if ~isempty(benchmarks(i).absolute_g)
                c(i) = 3;
            end
        end
        folder = k.createFolder('Benchmarks with zero height');

        folder.point(lon, lat, h, ...
            'iconScale', 0.5, ...
            'name', name, ...
            'iconURL', icons(ones(size(benchmarks,1),1)), ...
            'iconColor', colors(c), 'altitudeMode', 'clampToGround');
        
        k.save(fullfile(path,file))
        
        msgbox(['Successfully exported benchmarks without height to ' fullfile(path,file)]);
    end
    
