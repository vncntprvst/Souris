function [ephysCommonTrials, behaviorCommonTrials]=match_trials(ephysTrials,Behavior)

% load('PrV77_32__2016_01_30_02_09_00_OEph_CAR'); 
%clocktest_TTLs_OEformat_2016_02_11_15_06_07_OEph_CAR clocktest_1_2016_02_10_22_52_02_OEph_CAR
% Behavior=processBehaviorData;

% no need to downsample - already done
%downSamplingRatio=uint64(Spikes.samplingRate(1,1)/Spikes.samplingRate(1,2));

%msTTL_TrialStart=double(Trials.start)/double(downSamplingRatio);
TTL_TrialStart=ephysTrials.start;
zeroed_TTL_TrialStart=TTL_TrialStart-TTL_TrialStart(1)+1;
zeroed_TTL_TrialStartIdx=false(1,ceil(zeroed_TTL_TrialStart(end)));
zeroed_TTL_TrialStartIdx(int32(round(zeroed_TTL_TrialStart)))=true;

% msTTL_TrialEnd=double(Trials.end)/double(downSamplingRatio);
TTL_TrialEnd=ephysTrials.end;
zeroed_TTL_TrialEnd=TTL_TrialEnd-TTL_TrialStart(1)+1;
zeroed_TTL_TrialEndIdx=false(1,ceil(zeroed_TTL_TrialEnd(end)));
zeroed_TTL_TrialEndIdx(int32(round(zeroed_TTL_TrialEnd)))=true;

Behav_TrialStart=Behavior.trials.trialStartTime;%[Behavior.trialTime{Behavior.trialStartIdx}];
% Behav_TrialStart=Behav_TrialStart(1:end); %(2:end) if missed first trial
% msBehav_TrialStart=[Behavior.trialTime{1:end}]; %plot all timestamps for
% comparison purposes
zeroed_Behav_TrialStart=Behav_TrialStart-Behav_TrialStart(1)+1;
zeroed_Behav_TrialStartIdx=false(1,max([ceil(zeroed_TTL_TrialStart(end)) ceil(zeroed_Behav_TrialStart(end))]));
zeroed_Behav_TrialStartIdx(int32(round(zeroed_Behav_TrialStart)))=true;

% figure;hold on
% plot(zeroed_TTL_TrialStartIdx)
% plot(zeroed_TTL_TrialEndIdx)
% plot(zeroed_Behav_TrialStartIdx*0.5)

%% find common trials
% bin to 30ms bins
binSize=30;
numBin=ceil(size(zeroed_TTL_TrialStartIdx,2)/binSize);
[zeroed_TTL_TrialStartIdx_binned] = histcounts(double(find(zeroed_TTL_TrialStartIdx)), linspace(0,size(zeroed_TTL_TrialStartIdx,2),numBin));
[zeroed_Behav_TrialStartIdx_binned] = histcounts(double(find(zeroed_Behav_TrialStartIdx)), linspace(0,size(zeroed_Behav_TrialStartIdx,2),numBin));

figure;hold on
plot(zeroed_TTL_TrialStartIdx_binned)
plot(zeroed_Behav_TrialStartIdx_binned*0.5)

find(zeroed_TTL_TrialStartIdx_binned,10)
find(zeroed_Behav_TrialStartIdx_binned,10)

behaviorCommonTrials=ismember(find(zeroed_Behav_TrialStartIdx_binned),find(zeroed_TTL_TrialStartIdx_binned)) |...
ismember(find(zeroed_Behav_TrialStartIdx_binned),find(zeroed_TTL_TrialStartIdx_binned)+1) |...
ismember(find(zeroed_Behav_TrialStartIdx_binned),find(zeroed_TTL_TrialStartIdx_binned)-1);

ephysCommonTrials=ismember(find(zeroed_TTL_TrialStartIdx_binned),find(zeroed_Behav_TrialStartIdx_binned)) |...
ismember(find(zeroed_TTL_TrialStartIdx_binned),find(zeroed_Behav_TrialStartIdx_binned)+1) |...
ismember(find(zeroed_TTL_TrialStartIdx_binned),find(zeroed_Behav_TrialStartIdx_binned)-1);
end
