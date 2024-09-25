function varargout = fg_load_from_tablet(varargin)
% FG_LOAD_FROM_TABLET MATLAB code for fg_load_from_tablet.fig
%      FG_LOAD_FROM_TABLET, by itself, creates a new FG_LOAD_FROM_TABLET or raises the existing
%      singleton*.
%
%      H = FG_LOAD_FROM_TABLET returns the handle to a new FG_LOAD_FROM_TABLET or the handle to
%      the existing singleton*.
%
%      FG_LOAD_FROM_TABLET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_LOAD_FROM_TABLET.M with the given input arguments.
%
%      FG_LOAD_FROM_TABLET('Property','Value',...) creates a new FG_LOAD_FROM_TABLET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_load_from_tablet_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_load_from_tablet_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_load_from_tablet

% Last Modified by GUIDE v2.5 15-Sep-2024 20:26:35

% Franco S. Sobrero, Ohio State University, Sept. 2024

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_load_from_tablet_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_load_from_tablet_OutputFcn, ...
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


% --- Executes just before fg_load_from_tablet is made visible.
function fg_load_from_tablet_OpeningFcn(hObject, eventdata, handles, varargin)
    
    % Center the figure on the screen
    movegui(hObject, 'center');

    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_load_from_tablet (see VARARGIN)
    
    % This assigns handles.output to the GUI's figure handle (hObject), ensuring 
    % the field exists before it is accessed in OutputFcn.
    handles.output = []; %hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Combine PathName and FileName to get the full file path
    PathName = varargin{2};
    FileName = varargin{1};
    FullFilePath = fullfile(PathName, FileName);
    
    % Current line in the fg_edit_line interface
    loadedLine = varargin{3};
     
    if ~isempty(loadedLine.benchmarks)
        loadedBenchmarks = {loadedLine.benchmarks(:, 1).name}';
    else
        loadedBenchmarks = {};
    end

    % Use readtable to load the Excel file
    Text_struc = readtable(FullFilePath, 'Sheet', 1, 'VariableNamingRule', 'preserve');
    
    % Get the unique lines and corresponding data
    line_nr = Text_struc.line;
    lines_tablet = unique(line_nr);
    
    % Initialize variables for station start, end, and date
    station_start = cell(length(lines_tablet), 1);
    station_end = cell(length(lines_tablet), 1);
    date_tablet = cell(length(lines_tablet), 1);
    gravim_used = cell(length(lines_tablet), 1);
    operator = cell(length(lines_tablet), 1);
    line_status = cell(length(lines_tablet), 1);
    
    for k = 1:length(lines_tablet)
        % Get the data corresponding to the current line
        ff = ismember(line_nr, lines_tablet(k));
        benchmarks_tablet{k} = Text_struc.benchmark(ff);
        gravim_tablet = Text_struc.instrument(ff);
        operator_tablet = Text_struc.user(ff);
        line_status_tablet = Text_struc.("line status")(ff);
    
        % Assign the start and end stations
        station_start{k} = char(benchmarks_tablet{1,k}{1,1});               % Convert to char
        station_end{k} = char(benchmarks_tablet{1,k}{end/2,1});             % Convert to char
    
        % Assign the operator of the instrument for this line
        if length(unique(operator_tablet)) > 1
            operator{k} = "MULTPLE";
        else
            operator{k} = char(unique(operator_tablet));                    % Convert to char
        end
    
        % Assign the line_status for this line
        if length(unique(line_status_tablet)) > 1
            line_status{k} = 'ERROR';                                       % In principle this will never occur, but just in case
        elseif char(unique(line_status_tablet)) == "CERRADA"                % Convert to char
            line_status{k} = 'CLOSED';
        else
            line_status{k} = 'INCONSISTENT';
        end
    
        % Assign the gravimeter used in each line
        if length(unique(gravim_tablet)) > 1
            gravim_used{k} = 'ERROR';
        else
            gravim_used{k} = char(unique(gravim_tablet));                   % Convert to char
        end
    
        % Extract and format the date
        date_full_tablet = extractBefore(Text_struc.timestamp(ff), ".");
        date_tablet{k} = char(extractBefore(date_full_tablet{1}, 'T'));     % Convert to char
    end
    
    % Create a logical column for checkboxes (user checks rows of the lines they want to load onto fg_edit_line.m)
    load_column = false(length(lines_tablet), 1);                           % Initially unchecked
    
    if isempty(loadedBenchmarks)                                            % If no benchmarks are loaded, then any line is allowed to be loaded
        benchmark_match = true(length(lines_tablet), 1);                    % No benchmark conflict
        instrument_match = false(length(lines_tablet), 1);                  % No gravimeter conflict
    else
        for k = 1:length(lines_tablet)
            line_benchmarks = unique(benchmarks_tablet{1, k});
            benchmark_match(k) = all(ismember(line_benchmarks, loadedBenchmarks));              % Check if all benchmarks in this new line match the ones already loaded to the fg_edit_line.m interface
            instrument_match(k) = ismember(gravim_used{k}, unique(loadedLine.instruments(:)));  % Check if the gravimeter of this new line is matches any of the gravimeters already loaded to the fg_edit_line.m interface
            
        end
    end

    % Combine the data into a cell array
    data = [lines_tablet, station_start, station_end, gravim_used, operator, date_tablet, line_status, num2cell(load_column)];

    % Set the data and format in the uitable
    set(handles.listOfTabletLines, 'Data', data);
    set(handles.listOfTabletLines, 'ColumnName', {'Line', 'Start Station', 'End Station', 'Gravimeter', 'Operator', 'Date', 'Line Status', 'Load'});
    set(handles.listOfTabletLines, 'ColumnEditable', [false false false false false false false true]);
    set(handles.listOfTabletLines, 'ColumnFormat', {'char', 'char', 'char', 'char', 'char', 'char', 'char', 'logical'});
    
    % Store benchmark matching information in handles
    handles.benchmark_match = benchmark_match;
    handles.instrument_match = instrument_match;
       
    colorRows(hObject, handles);                                            % Set the row colors
    guidata(hObject, handles);                                              % Update handles structure
    uiwait(handles.figure1);                                                % Wait for the user to either Load selected lines, or Cancel the operation

% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
    tableData = get(handles.listOfTabletLines, 'Data');                     % Retrieve the table data
    selectedRows = cell2mat(tableData(:, end));                             % Get the rows where the checkbox is selected, the logical column (checkbox column) is the last one

    if any(selectedRows)
        selectedLines = tableData(selectedRows == true, 1);                 % Get the selected lines
        selectedLinesFeatures = tableData(selectedRows == true, 2:4);       % Export the First station, Last station, and gravimeter of each selected line, for fg_line_edit.m to check if they can all be loaded together
        n = size(selectedLinesFeatures, 1);                                 % Number of selected lines to be loaded
        first_col_equal = all(cellfun(@(x) isequal(x, selectedLinesFeatures{1,1}), selectedLinesFeatures(:,1)));    % Check if the first station is the same at all lines
        second_col_equal = all(cellfun(@(x) isequal(x, selectedLinesFeatures{1,2}), selectedLinesFeatures(:,2)));   % Check if the last station is the same at all lines
        third_col_diff = numel(unique(selectedLinesFeatures(:,3))) == n;                                            % Check if all lines were surveyed with different gravimeters
        
        % First check if the lines I am trying to load have the same benchmarks as the already loaded lines
        if ~all(handles.benchmark_match(selectedRows)) || any(handles.instrument_match(selectedRows))       % If any of the selected lines is in conflict with alreaded loaded lines (differente stations, or repeated gravimeters)
            errordlg('The selected line(s) cannot be loaded due to a conflict with existing line(s). This may be due to mismatched stations or repeated use of gravimeter(s).', 'Invalid Line Selection', 'modal');
            handles.output = [];                                                % No selection made, set output to empty
            guidata(hObject, handles);                                          % Update handles structure
        elseif ~first_col_equal || ~second_col_equal || ~third_col_diff         % If all selected lines they share the same First and Last stations, and if they were measured with different gravimeters
            errordlg('All selected lines must contain the same stations and use different gravimeters', 'Invalid line selection', 'modal');
            handles.output = [];                                                % No selection made, set output to empty
            guidata(hObject, handles);                                          % Update handles structure
        else
            handles.output = selectedLines;                                     % Store selected lines in handles.output
            guidata(hObject, handles);                                          % Update handles structure
            uiresume(handles.figure1);                                          % Resume execution of the parent function
        end
    else
        errordlg('No gravity lines selected', 'Empty selection', 'modal');
        handles.output = [];                                                % No selection made, set output to empty
        guidata(hObject, handles);                                          % Update handles structure
    end

    uiresume(gcf);  % Resume execution of the parent function

% --- Outputs from this function are returned to the command line.
function varargout = fg_load_from_tablet_OutputFcn(hObject, eventdata, handles)
    % Check if handles.output exists and is not empty
    if isfield(handles, 'output') && ~isempty(handles.output)
        % If output is valid, pass it to varargout
        varargout{1} = handles.output;
    else
        % If no output (e.g., Cancel was pressed or no lines selected), return empty
        varargout{1} = [];
    end

    % Try to close the figure if it still exists
    try
        if ishandle(hObject)
            delete(hObject);
        end
    catch ME
        rethrow(ME);                                                        % If an error occurs, show it to me
    end

% --- This function is used to color the table rows based on the line status
% RED: Line with some serious problem, such as: multiple gravimeters in 
% the same line, multiple operators, unknown gravimeter, etc
% YELLOW: Line exported by the tablet app as "INCONSISTENT", which means 
% that it was not closed. Needs inspection before loading it to the 
% adjustment software
function colorRows(hObject, handles)
    % Get the table data
    data = get(handles.listOfTabletLines, 'Data');

    % Get the number of rows and columns
    [numRows, numCols] = size(data);

    % Initialize the color matrix
    colors = 0.97*(ones(numRows, 3));                                           % Gray background for all rows
    
    % Set colors based on line_status (column 7) and benchmark matching
    for row = 1:numRows
        status = data{row, 7};  % Assuming line_status is in the 7th column
        if handles.benchmark_match(row) && ~handles.instrument_match(row)       % If this line's benchamarks match the ones already loaded, and no gravimeter conflict, mark this line in green (as allowed to be loaded)
            colors(row, :) = [0.8, 1, 0.8];                                       % Light red for non-matching benchmarks, or gravimeters conflict
        end
    end

    % Apply the colors to the table
    set(handles.listOfTabletLines, 'BackgroundColor', colors);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    % This function is called when the user tries to close the figure
    handles.output = [];                                                    % Set output to empty
    guidata(hObject, handles);                                              % Update handles structure
    %uiresume(gcbf);                                                        % Resume execution of the parent function using gcbf
    % Hint: delete(hObject) closes the figure
    delete(hObject);
