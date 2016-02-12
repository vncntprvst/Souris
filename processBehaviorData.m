
function Behavior=processBehaviorData(fileName,dirName)
% fileName='C:\Data\Behav\PrV77_32_2016-01-30T02_06_59.csv';
behavDir='C:\Data\Behav';

switch nargin
    case 0
        [fileName,dirName] = uigetfile({'*.csv','.csv Files';...
            '*.*','All Files' },'Behavior Data',behavDir);
        cd(dirName)
    case 1
        cd(behavDir)
    case 2
        cd(dirName)
end

Behavior=readTrialData(fileName,'C:\Data\Behav');

%convert trial times to milliseconds, with exact sub-ms precision
Behavior.trialTime=cellfun(@(x) (floor(seconds(hours(str2double(x{4})) + ...
    minutes(str2double(x{5})) + ...
    seconds(str2double(x{6})))) + ...
    str2double(['0.' x{7}]))*1000, Behavior.trialTime,'UniformOutput',false);

%intervals
Behavior.intervals=[0;diff(vertcat(Behavior.trialTime{:}))];

%index trials
Behavior.trialStartIdx=Behavior.trialType==0;

%% failed trials
% old recordings: no failed/timeout trials - bug gave ~ unlimited time to get reward
%  although a spurious 0 was inserted right after trial (without TTLout(2))
% 2/11/16 +> changed that in Training2_Front_then_Side_Ports_PCcontrol_2FrontLEDs
% Timeout after 10 seconds. Timeout code added (90).
if datetime < Behavior.fileRecodingDate
    if sum(diff(Behavior.trialType)==0)
        Behavior.trialStartIdx(diff(Behavior.trialType)==0)=0;
    end
    Behavior.failedTrials=ceil(find([diff(Behavior.trialType)==0;0]&Behavior.trialStartIdx)/2);
    
    Behavior.trialStartTime=[Behavior.trialTime{Behavior.trialStartIdx}]';
else
    Behavior.timeoutTrial=find(Behavior.trialType(logical([0;Behavior.trialStartIdx(1:end-1)]))==90);
end
end

