
msTTL_TrialStart=double(Trials.start)/double(downSamplingRatio);
zeroed_msTTL_TrialStart=msTTL_TrialStart-msTTL_TrialStart(1)+1;
zeroed_msTTL_TrialStartIdx=false(1,size(zeroed_msTTL_TrialStart,1));
zeroed_msTTL_TrialStartIdx(int32(round(zeroed_msTTL_TrialStart)))=true;

msBehav_TrialStart=Behavior.trialStartTime(2:end); % missed the first trial on that one
zeroed_msBehav_TrialStart=msBehav_TrialStart-msBehav_TrialStart(1)+1;
zeroed_msBehav_TrialStartIdx=false(1,size(zeroed_msBehav_TrialStart,1));
zeroed_msBehav_TrialStartIdx(int32(round(zeroed_msBehav_TrialStart)))=true;
zeroed_msBehav_TrialStartIdx=zeroed_msBehav_TrialStartIdx*0.5;

figure;hold on
plot(zeroed_msTTL_TrialStartIdx)
plot(zeroed_msBehav_TrialStartIdx)


