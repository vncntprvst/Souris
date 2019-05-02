recDir='D:\Data\vIRt32\vIRt32_0424\vIRt32_2019-04-24_16-48-53_5185\experiment1\recording1\continuous\Rhythm_FPGA-100.0';
recName='continuous.dat';
[rec,data,spikes,TTLs] = LoadEphysData(recName,recDir);

spikeIds=unique(spikes.clusters);
spikesTimes=spikes.spikeTimes(spikes.clusters==1)-rec.recordingStartTime;
bestTrace=mode(spikes.electrodes(spikes.clusters==1));

figure; hold on
% plot(data(bestTrace,:))
% plot(spikesTimes,ones(numel(spikesTimes),1),'db')
plot(spikesTimes-vidTimes(1),zeros(numel(spikesTimes),1),'dk')

sum(TTLs(2,:)==2)
vidTimes=int64(TTLs(1,TTLs(2,:)<0))-rec.recordingStartTime;
% vidTimes=324958;
load('whiskerTrackingData.mat')


angleTrace=whiskerTrackingData(1,:); %1:numel(vidTimes) %size(whiskerTrackingData,2)-numel(vidTimes)+1:end
vidTimeArray=linspace(1,double(vidTimes(end)-vidTimes(1)+1),size(angleTrace,2));
plot(vidTimeArray,angleTrace-mean(angleTrace)); %vidTimes-vidTimes(1)
plot(spikesTimes-vidTimes(1),zeros(numel(spikesTimes),1),'dk')
angleTrace=whiskerTrackingData(1,1:numel(vidTimes)); % %size(whiskerTrackingData,2)-numel(vidTimes)+1:end
plot(vidTimes-vidTimes(1),angleTrace-mean(angleTrace)); 


% tsbackupID=fopen('vIRt32_0424_5185_pole.mp4');
% tsTTLs=fread(tsbackupID,[2 183377997],'int64');
% fclose(tsbackupID)