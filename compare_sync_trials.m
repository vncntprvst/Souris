load('PrV77_32__2016_01_30_02_09_00_OEph_CAR'); %clocktest_1_2016_02_10_22_52_02_OEph_CAR
Behavior=processBehaviorData;

downSamplingRatio=uint64(Spikes.samplingRate(1,1)/Spikes.samplingRate(1,2));


msTTL_TrialStart=double(Trials.start)/double(downSamplingRatio);
zeroed_msTTL_TrialStart=msTTL_TrialStart-msTTL_TrialStart(1)+1;
zeroed_msTTL_TrialStartIdx=false(1,size(zeroed_msTTL_TrialStart,1));
zeroed_msTTL_TrialStartIdx(int32(round(zeroed_msTTL_TrialStart)))=true;

msTTL_TrialEnd=double(Trials.end)/double(downSamplingRatio);
zeroed_msTTL_TrialEnd=msTTL_TrialEnd-msTTL_TrialStart(1)+1;
zeroed_msTTL_TrialEndIdx=false(1,size(zeroed_msTTL_TrialEnd,1));
zeroed_msTTL_TrialEndIdx(int32(round(zeroed_msTTL_TrialEnd)))=true;

msBehav_TrialStart=Behavior.trialStartTime(1:end);
% msBehav_TrialStart=[Behavior.trialTime{1:end}]; %plot all timestamps for
% comparison purposes
zeroed_msBehav_TrialStart=msBehav_TrialStart-msBehav_TrialStart(1)+1;
zeroed_msBehav_TrialStartIdx=false(1,size(zeroed_msBehav_TrialStart,1));
zeroed_msBehav_TrialStartIdx(int32(round(zeroed_msBehav_TrialStart)))=true;
zeroed_msBehav_TrialStartIdx=zeroed_msBehav_TrialStartIdx*0.5;


figure;hold on
plot(zeroed_msTTL_TrialStartIdx)
plot(zeroed_msTTL_TrialEndIdx)
plot(zeroed_msBehav_TrialStartIdx)



