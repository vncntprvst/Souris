function [ephys,behav]=GetData_From_UnitExplorer(ue)
%% get whisker data
wAngle = cell2mat(ue.dataTable.position);
wVelocity = cell2mat(ue.dataTable.velocity);
wPhase= cell2mat(ue.dataTable.phase);

% fill missing values
wAngle=fillmissing(wAngle,'spline','EndValues','nearest');
wVelocity=fillmissing(wVelocity,'spline','EndValues','nearest');
wPhase=fillmissing(wPhase,'spline','EndValues','nearest');
    
%resample to ms precision
videoSamplingRate=ue.videoSamplingRate;
ephysSamplingRate=1000;

% timeVector=(0:ephysSamplingRate/videoSamplingRate:numel(wAngle)*...
%     ephysSamplingRate/videoSamplingRate-1)/1000;
% wAngle=resample(wAngle,ephysSamplingRate,videoSamplingRate);
% wVelocity=resample(wVelocity,ephysSamplingRate,videoSamplingRate);
% foo=resample(wPhase,ephysSamplingRate,videoSamplingRate,2,40);
% wPhase=resample(wPhase',timeVector,ephysSamplingRate,'spline');
% wPhase = interp(wPhase,2);
% foo=WhiskingFun.ComputePhase(wAngle');
% iwPhase=nan(numel(wAngle)*ephysSamplingRate/videoSamplingRate,1);
% iwPhase(1:ephysSamplingRate/videoSamplingRate:end)=wPhase;

%% get spike times and compute rasters
spikeTimes = cell(length(ue.dataTable.spikeTimes),1);
for epochNum = 1: length(ue.dataTable.spikeTimes)
    epochSpikeTimes = ue.dataTable.spikeTimes{epochNum,1};
    spikeTimes{epochNum} = (epochSpikeTimes + 5.0 * (epochNum - 1));
end
spikeTimes=int32(vertcat(spikeTimes{:})*ephysSamplingRate); %from s. to ms.
unitID=ones(numel(spikeTimes),1);
spikeRasters=EphysFun.MakeRasters(spikeTimes,unitID,ephysSamplingRate,...
    numel(wAngle)*(ephysSamplingRate/videoSamplingRate)); %if no resampling:  numel(wAngle)*(ephysSamplingRate/videoSamplingRate)

ephys=struct('spikeTimes',spikeTimes,'spikeRasters',spikeRasters,'samplingRate',ephysSamplingRate);
behav=struct('angle',wAngle,'velocity',wVelocity,'phase',wPhase,'samplingRate',videoSamplingRate);









