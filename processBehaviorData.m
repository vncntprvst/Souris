
function Behavior=processBehaviorData(fileName,dName)
% fileName='C:\Data\Behav\PrV77_32_2016-01-30T02_06_59.csv';
behavDir='C:\Data\Behav';

switch nargin
    case 0
        [fileName,dName] = uigetfile({'*.csv','.csv Files';...
            '*.*','All Files' },'Behavior Data',behavDir);
        cd(dName)
    case 1
        cd(behavDir)
    case 2
        cd(dName)
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
% 2/11/16 - currently, no failed/timeout trials - bug gives ~ unlimited time to get reward
%  but a spurious 0 is inserted right after trial (without TTLout(2)) +> changed that now. 
if sum(diff(Behavior.trialType)==0)
    Behavior.trialStartIdx(diff(Behavior.trialType)==0)=0;
end
Behavior.failedTrials=ceil(find([diff(Behavior.trialType)==0;0]&Behavior.trialStartIdx)/2);

Behavior.trialStartTime=[Behavior.trialTime{Behavior.trialStartIdx}]';

end








