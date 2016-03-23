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
csvFile.fileRecodingDate=datetime([fileStartTime{1} '-' fileStartTime{2} '-' fileStartTime{3}]);
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
csvFile.trialNumber = dataArray{:, 1};
csvFile.trialType = dataArray{:, 2};
csvFile.successCount = dataArray{:, 3};
csvFile.eventTime=cellfun(@(x) regexp(x,'\d+','match'),dataArray{:, 4},'UniformOutput',false);
csvFile.eventTime_ms=cellfun(@(x) 1000*(seconds(datetime(str2double(x{1}),str2double(x{2}),str2double(x{3}))...
    -datetime(str2double(x{1}),1,1))+ seconds(duration(str2double(x{4}),str2double(x{5}),str2double([x{6} '.' x{7}])))),...
    csvFile.trialTime,'UniformOutput',false);
% csvFile.eventTime_ms=cellfun(@(x) 1000*(seconds(duration(str2double(x{4}),str2double(x{5}),str2double([x{6} '.' x{7}])))),...
%     csvFile.trialTime,'UniformOutput',false);
csvFile.eventTime_ms=cellfun(@(x) x-csvFile.eventTime_ms{1}, csvFile.eventTime_ms,'UniformOutput',false);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans fileStartTime;