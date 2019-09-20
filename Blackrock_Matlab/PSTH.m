function varargout = PSTH(varargin)
% PSTH MATLAB code for PSTH.fig
%      PSTH, by itself, creates a new PSTH or raises the existing
%      singleton*.
%
%      H = PSTH returns the handle to a new PSTH or the handle to
%      the existing singleton*.
%
%      PSTH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSTH.M with the given input arguments.
%
%      PSTH('Property','Value',...) creates a new PSTH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PSTH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PSTH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PSTH

% Last Modified by GUIDE v2.5 01-Jul-2019 20:08:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PSTH_OpeningFcn, ...
                   'gui_OutputFcn',  @PSTH_OutputFcn, ...
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


% --- Executes just before PSTH is made visible.
function PSTH_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PSTH (see VARARGIN)

% Choose default command line output for PSTH
handles.output = hObject;

% --- Begin My Code ---

% Create a timer object that will be used to grab data and refresh
% analysis/plotting
handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly
    'Period', 0.5, ...                      % Initial period is 100 ms
    'TimerFcn', {@updateDisplay,hObject}, ... % callback function.  Pass the figure handle
    'StartFcn', {@startTimer,hObject});     % callback to execute when timer starts

handles.cbmexStatus = 'closed';

% Update handles structure
guidata(hObject, handles);

clc
clear all

% UIWAIT makes PSTH wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PSTH_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----        Figure Objects Create and Callback Functions            --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% --- Executes on button press in cmd_cbmexOpen.
function cmd_cbmexOpen_Callback(hObject, eventdata, handles)
% hObject    handle to cmd_cbmexOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use a TRY-CATCH in case cbmex is already open.  If you try to open it
% when it's already open, Matlab throws a 'MATLAB:unassignedOutputs'
% MException.
try
    cbmex('open');
catch ME
    if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
        % Dont need to do anything because cbmex is already open and it
        % already sends a message stating that
    else
        disp(ME)
    end
end
handles.cbmexStatus = 'open';

cbmex('trialconfig',1,'absolute')
pause(0.1)

% Acquire some data to get channel information.  Determine which channels
% are enabled
[spikeEvents, time, continuousData] = cbmex('trialdata',1);
handles.channelList = [continuousData{:,1}];
% set channel popup meno to hold channels
set(handles.pop_channels,'String',handles.channelList);

guidata(hObject,handles)

% --- Executes on button press in cmd_cbmexClose.
function cmd_cbmexClose_Callback(hObject, eventdata, handles)
% hObject    handle to cmd_cbmexClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cbmex('close')
handles.cbmexStatus = 'closed';
guidata(hObject,handles)

function txt_lfpPreSpike_Callback(hObject, eventdata, handles)
% hObject    handle to txt_lfpPreSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_lfpPreSpike as text
%        str2double(get(hObject,'String')) returns contents of txt_lfpPreSpike as a double
handles.windowSize = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function txt_lfpPreSpike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_lfpPreSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tgl_StartStop.
function tgl_StartStop_Callback(hObject, eventdata, handles)
% hObject    handle to tgl_StartStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tgl_StartStop

handles = guidata(hObject);

% if Start
if get(hObject,'Value') == 0
    
    % Check to make sure cbmex connection is open
    if strcmp(handles.cbmexStatus,'closed')
        errordlg('No cbmex connection.  Open connection before starting','Not Connected')
        return
    end
    
    set(hObject,'String','Stop');

    % This starts the timer and also executes the StartFnc which grabs the
    % data, creates the buffer and plots the first bit of data
    start(handles.timer)
    
% Stop
else
    set(hObject,'String','Start')
    stop(handles.timer)
end

function txt_lfpPostSpike_Callback(hObject, eventdata, handles)
% hObject    handle to txt_lfpPostSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_lfpPostSpike as text
%        str2double(get(hObject,'String')) returns contents of txt_lfpPostSpike as a double


% --- Executes during object creation, after setting all properties.
function txt_lfpPostSpike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_lfpPostSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_channels.
function pop_channels_Callback(hObject, eventdata, handles)
% hObject    handle to pop_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_channels
settingChange(hObject)

% --- Executes during object creation, after setting all properties.
function pop_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cbmex('close')

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard

% --- Executes on selection change in pop_unit.
function pop_unit_Callback(hObject, eventdata, handles)
% hObject    handle to pop_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_unit
settingChange(hObject)

% --- Executes during object creation, after setting all properties.
function pop_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----          Main Loop, Runs every time the timer fires            --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %
function updateDisplay(hObject, eventdata, hfigure)
try
    
    handles = guidata(hfigure);
    if strcmp(handles.cbmexStatus,'closed')
        stop(handles.timer)
    end
    
    [events, time, continuousData] = cbmex('trialdata',1);

    newSpikeTimes = events{handles.channelIndex, handles.unitIndex+1};
    % make sure it's a column vector
    if size(newSpikeTimes,2) ~= 1
        newSpikeTimes = newSpikeTimes';
    end
    newContinuousData = continuousData{handles.channelIndex,3};
    handles.rawDataBuffer = cycleBuffer(handles.rawDataBuffer, newContinuousData);
    handles.lastSampleProcTime = time*30000 + length(newContinuousData) - 1;

    spikeTimes = [handles.unprocessedSpikes; newSpikeTimes];

    guidata(hfigure,handles)
    

    if ~isempty(spikeTimes)
        [lfps, unprocessedEvents] = extractLfp(spikeTimes, hfigure);

        if ~isempty(lfps)
            % update the average LFP based on the previous average, the previous N
            % and the new waveforms and their quantity.  If there 'lfps' is
            % empty, the unprocessedEvents will get added to the list of
            % unprocessedSpikes field below and retried
            handles.lfpAverage = (handles.lfpAverage*handles.numLfp + sum(lfps,2)) / ...
                                  (handles.numLfp + size(lfps,2));
    %             handles.lfpMatrix = [handles.lfpMatrix, lfps];

            handles.numLfp = handles.numLfp + size(lfps,2);

    %         handles.averageLfp = mean(handles.lfpMatrix,2);
    %         handles.stdLfp = std(double(handles.lfpMatrix),[],2);      
    %         [lfpPatchX, lfpPatchY] = createPatchWaveform(handles.averageLfp, handles.stdLfp);
    %         set(handles.h_lfpPatch,'XData', lfpPatchX, 'YData', lfpPatchY);

            set(handles.h_lfps, 'YData', handles.lfpAverage);
            a = axis(handles.ax_rasters);
            set(handles.h_lfpN, 'Position',[(a(2)-a(1))*0.9+a(1), (a(4)-a(3))*0.9+a(3)], ...
                'String', ['N = ' num2str(handles.numLfp)]);

        end

        handles.unprocessedSpikes = unprocessedEvents;

    end

    % update YData of ax_rawData
    set(handles.h_rawDataTrace,'YData',handles.rawDataBuffer)

    guidata(hfigure,handles)

catch ME
    getReport(ME)
end
    
% ----------------------------------------------------------------------- %
% ----                                                                --- %
% ----                      Helper Functions                          --- %
% ----                                                                --- %
% ----------------------------------------------------------------------- %

% Runs once when timer starts
function  startTimer(hObject, eventdata, hfigure)
try
    handles = guidata(hfigure);

    % Create raw data buffer of zeros
    handles.rawDataBuffer = zeros(150000,1);

    lfpTraceLength = (str2double(get(handles.txt_lfpPreSpike,'String')) + ...
                     str2double(get(handles.txt_lfpPostSpike,'String'))) * 30 + 1;
    handles.lfpAverage = zeros(lfpTraceLength,1);
    handles.lfpMatrix = [];
    handles.numLfp = 0;

    % Check which channel is selected and get some data to plot
    handles.channelIndex = get(handles.pop_channels,'Value');
    handles.unitIndex = get(handles.pop_unit,'Value');

    [events, time, continuousData] = cbmex('trialdata',1);
    newSpikeTimes = events{handles.channelIndex,handles.unitIndex+1};
    newContinuousData = continuousData{handles.channelIndex,3};

    % keep track of the proc time of the most recent data point.  This will
    % help pull out chunks based off spike times.  'time' is the time at the
    % first data point of the new chunk of continuous data in seconds.
    handles.lastSampleProcTime = time*30000 + length(newContinuousData) - 1;

    handles.rawDataBuffer = cycleBuffer(handles.rawDataBuffer, newContinuousData);
    handles.h_rawDataTrace = plot(handles.ax_raw,handles.rawDataBuffer);
    handles.h_lfps = plot(handles.ax_rasters, handles.lfpAverage);
    a = axis(handles.ax_rasters);
    handles.h_lfpN = text((a(2)-a(1))*0.9+a(1), (a(4)-a(3))*0.9+a(3), ...
        'N = 0','Parent',handles.ax_rasters);

    % create zeros patch plot for LFP.  Use handles.lfpAverage because it's a
    % zeros vector of the right size.
    % [XPatch, YPatch] = createPatchWaveform(handles.lfpAverage,handles.lfpAverage);
    % handles.h_lfpPatch = patch(XPatch,YPatch,'b');

    % Unless the LFP window is super small, the first chunk of data won't be
    % big enough for extracting LFP around spikes.  Add any captured spikes to
    % the list of unprocessed spikes for later
    handles.unprocessedSpikes = newSpikeTimes;

    guidata(hfigure,handles)
    
catch ME
    getReport(ME)
end

function newBuffer = cycleBuffer(oldBuffer, newData)
N = length(newData);
if N >= length(oldBuffer)
    newBuffer = newData(end-N+1:end);
else
    newBuffer = [oldBuffer(N+1:end); newData];
end


function [lfpMatrix, unprocessedEvents] = extractLfp(spikeTimes, hfigure)

handles = guidata(hfigure);

window = [str2double(get(handles.txt_lfpPreSpike,'String')), ...
          str2double(get(handles.txt_lfpPostSpike,'String'))] * 30;
      
unprocessedEvents = [];
lfpMatrix = [];

rawDataStart = handles.lastSampleProcTime - length(handles.rawDataBuffer) + 1;
rawDataEnd = handles.lastSampleProcTime;

% if using nPlay, when you restart the file, the proctime will be reset and
% rawDataStart becomes negative.  Save up spike times and process them when
% it's back to normal
if rawDataStart < 0
    unprocessedEvents = spikeTimes;
    return
end

for i = 1 : length(spikeTimes)
    
    if spikeTimes(i)-window(1) < rawDataStart || spikeTimes(i)+window(2) > rawDataEnd
        unprocessedEvents = [unprocessedEvents; spikeTimes(i)];
    else
        if length(spikeTimes(i)-window(1):spikeTimes(i)+window(2)) ~= (window(1)+window(2)+1)
            disp(['LFP length wrong: ' ...
                num2str(length(spikeTimes(i)-window(1):spikeTimes(i)+window(2))) ])
            keyboard
        end
        
        lfpMatrix = [lfpMatrix, ...
                     handles.rawDataBuffer( ...
                     (spikeTimes(i)-window(1):spikeTimes(i)+window(2)) - rawDataStart-1 ) ...
                     ];
    end
end

function [patchX, patchY] = createPatchWaveform(meanVector, stdVector)

if size(meanVector,2) == 1 % is column
    patchY = [meanVector+stdVector; meanVector(end:-1:1)-stdVector(end:-1:1)];
    patchX = [(1:length(meanVector))'; (length(meanVector):-1:1)'];
else
    patchX = [1:length(meanVector), length(meanVector):-1:1];
    patchY = [meanVector+stdVector, meanVector(end:-1:1)-stdVector(end:-1:1)];
end

function settingChange(hObject)
handles = guidata(hObject);

% if the timer is running, stop it and restart it (which will use the newly
% selected channel.  If the timer isn't running, don't do anything.
if strcmp(handles.timer.Running,'on')
    stop(handles.timer)
    start(handles.timer)
end
