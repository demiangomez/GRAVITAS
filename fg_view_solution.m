function varargout = fg_view_solution(varargin)
% FG_VIEW_SOLUTION MATLAB code for fg_view_solution.fig
%      FG_VIEW_SOLUTION, by itself, creates a new FG_VIEW_SOLUTION or raises the existing
%      singleton*.
%
%      H = FG_VIEW_SOLUTION returns the handle to a new FG_VIEW_SOLUTION or the handle to
%      the existing singleton*.
%
%      FG_VIEW_SOLUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FG_VIEW_SOLUTION.M with the given input arguments.
%
%      FG_VIEW_SOLUTION('Property','Value',...) creates a new FG_VIEW_SOLUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fg_view_solution_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fg_view_solution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fg_view_solution

% Last Modified by GUIDE v2.5 27-Oct-2017 10:58:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fg_view_solution_OpeningFcn, ...
                   'gui_OutputFcn',  @fg_view_solution_OutputFcn, ...
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


% --- Executes just before fg_view_solution is made visible.
function fg_view_solution_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to fg_view_solution (see VARARGIN)

    % Choose default command line output for fg_view_solution
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    
    h = findobj('Tag','fg_adjust');
    
    hObj = guidata(h);

    showData(hObj, handles)
    
    vars = guidata(hObject);
    vars.Adjustment = hObj.Adjustment;
    
    guidata(hObject, vars);
        
    % UIWAIT makes fg_view_solution wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

function showData(hObj, handles)

    out(:,1) = {hObj.Adjustment.benchmarks.name}';
    out(:,2) = {hObj.Adjustment.benchmarks.lat}';
    out(:,3) = {hObj.Adjustment.benchmarks.lon}';
    out(:,4) = {hObj.Adjustment.benchmarks.height}';
    out(:,5) = num2cell(hObj.Adjustment.adjusted_g);
    out(:,6) = num2cell(hObj.Adjustment.adjusted_g_sigma);

    % sort by station name
    if handles.cboSort.Value == 1
        [~,c] = sortrows(out,1);
    else
        [~,c] = sortrows(out,6,'descend');
    end

    % sort according to selected output
    out = out(c,:);

    % remove any NaN stations (no adjusted gravity)
    if handles.chkNoGRemove.Value
        out = out(~cellfun(@isnan, out(:,5)),:);
    end

    % remove stations that have ~ zero height
    if handles.chkZeroHRemove.Value
        out = out(round(cell2mat(out(:,4))) ~= 0,:);
    end
    
    handles.tblAdj.Data = out;

% --- Outputs from this function are returned to the command line.
function varargout = fg_view_solution_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in cboSort.
function cboSort_Callback(hObject, eventdata, handles)
    % hObject    handle to cboSort (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns cboSort contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from cboSort
    hObj = guidata(hObject);
    showData(hObj, handles)

% --- Executes during object creation, after setting all properties.
function cboSort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cboSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkNoGRemove.
function chkNoGRemove_Callback(hObject, eventdata, handles)
    % hObject    handle to chkNoGRemove (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of chkNoGRemove
    hObj = guidata(hObject);
    showData(hObj, handles)

% --- Executes on button press in chkNoGRemove.
function chkZeroHRemove_Callback(hObject, eventdata, handles)
    % hObject    handle to chkNoGRemove (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of chkNoGRemove
    hObj = guidata(hObject);
    showData(hObj, handles)
    
