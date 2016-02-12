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
csvFile.fileStartTime=hours(str2num(fileStartTime{4})) + ...
    minutes(str2num(fileStartTime{5})) + ...
    seconds(str2num([fileStartTime{6} '.' fileStartTime{7}]));
csvFile.fileStartTimeSubMilli=str2num(fileStartTime{7});
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
csvFile.trialTime=cellfun(@(x) regexp(x,'\d+','match'),dataArray{:, 4},'UniformOutput',false);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans fileStartTime;