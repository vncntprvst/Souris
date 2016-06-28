function csvFile=readTrialData(filename,dname)
%% Import trial data from behavior csv file.
% filename = 'C:\Data\Behav\PrV77_32_2016-01-30T02_06_59.csv';

%% Initialize variables.
switch nargin
    case 0
        [filename,dname] = uigetfile({'*.csv','.csv Files';...
            '*.*','All Files' },'Behavior Data');
        cd(dname)
    case 1
        cd('C:\Data\Behav')
    case 2
        cd(dname)
end
delimiter = ',';
startRow = 2;

%% Open .csv file.
fileID = fopen(filename,'r');

%% get file open time from first line
fileStartTime=regexp(fgets(fileID),'\d+','match');
csvFile.fileRecordingDate=datetime([fileStartTime{1} '-' fileStartTime{2} '-' fileStartTime{3}]);
csvFile.fileStartTime=hours(str2double(fileStartTime{4})) + ...
    minutes(str2double(fileStartTime{5})) + ...
    seconds(str2double([fileStartTime{6} '.' fileStartTime{7}]));
csvFile.fileStartTimeSubMilli=str2double(fileStartTime{7});
frewind(fileID);

%% Read data 
formatSpec = '%f%f%f%s%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close file.
fclose(fileID);

%% Allocate imported array to column variable names
csvFile.trialNumberIdx = dataArray{:, 1};
csvFile.trialEventType = dataArray{:, 2};
csvFile.successCount = dataArray{:, 3};
csvFile.eventTime=cellfun(@(x) regexp(x,'\d+','match'),dataArray{:, 4},'UniformOutput',false);
%convert trial times to milliseconds, with exact sub-ms precision
csvFile.eventTime_ms=cellfun(@(x) 1000*(seconds(datetime(str2double(x{1}),str2double(x{2}),str2double(x{3}))...
    -datetime(str2double(x{1}),1,1))+ seconds(duration(str2double(x{4}),str2double(x{5}),str2double([x{6} '.' x{7}])))),...
    csvFile.eventTime,'UniformOutput',false);
csvFile.eventTime_ms=cellfun(@(x) x-csvFile.eventTime_ms{1}, csvFile.eventTime_ms,'UniformOutput',false);
% csvFile.eventTime_ms=cellfun(@(x) 1000*(seconds(duration(str2double(x{4}),str2double(x{5}),str2double([x{6} '.' x{7}])))),...
%     csvFile.trialTime,'UniformOutput',false);
% cellfun(@(x) (floor(seconds(hours(str2double(x{4})) + ...
%     minutes(str2double(x{5})) + ...
%     seconds(str2double(x{6})))) + ...
%     str2double(['0.' x{7}]))*1000, csvFile.eventTime,'UniformOutput',false);

%% Clear temporary variables
clearvars delimiter startRow formatSpec fileID dataArray ans fileStartTime;

%Identify trial starts and ends
eventSequence=diff(mod(csvFile.trialNumberIdx,2));
if eventSequence(end)~=0
lastEvent=find([1;eventSequence]==0,1,'last');
eventSequence=eventSequence(1:lastEvent-1);
csvFile.trialNumberIdx=csvFile.trialNumberIdx(1:lastEvent);
csvFile.trialEventType=csvFile.trialEventType(1:lastEvent);
             csvFile.successCount=csvFile.successCount(1:lastEvent);
                csvFile.eventTime=csvFile.eventTime(1:lastEvent);
             csvFile.eventTime_ms=csvFile.eventTime_ms(1:lastEvent);
end
csvFile.trials=struct('trialNumber',unique(csvFile.trialNumberIdx),...
    'trialRetroIndex',find([1;eventSequence]~=0),...
    'trialStartTime',[csvFile.eventTime_ms{[1;eventSequence]~=0}]',...
    'trialEndTime',[csvFile.eventTime_ms{[eventSequence;1]~=0}]',...
    'trialType',csvFile.trialEventType([1;eventSequence]~=0),...
    'trialOutcomeType',[csvFile.trialEventType([false;[1;diff(mod(csvFile.trialNumberIdx(1:end-1),2))]~=0]),... %first column: First answer
                        csvFile.trialEventType([eventSequence;1]~=0)]);  %2nd column: eventual answer (only different in relaxed constraint conditions

%% failed trials
% old recordings: no failed/timeout trials - bug gave ~ unlimited time to get reward
%  although a spurious 0 was inserted right after trial (without TTLout(2))
% 2/11/16 +> changed that in Training2_Front_then_Side_Ports_PCcontrol_2FrontLEDs
% Timeout after 10 seconds. Timeout code added (9x). x stands for task type (1/2 for left/right)
% 3/16/16 +> Incorrect reward port code added (8x) in TextureDiscriminationTask.
csvFile.trials.errorTrialIdx=round(csvFile.trials.trialOutcomeType/10)==8;
csvFile.trials.timeoutTrialIdx=round(csvFile.trials.trialOutcomeType/10)==9;
csvFile.trials.correctTrialIdx=~(csvFile.trials.errorTrialIdx |csvFile.trials.timeoutTrialIdx);

end