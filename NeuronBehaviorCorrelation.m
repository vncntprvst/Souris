%% Correlation between bursts/spike rate and periodical behaviors (whisking, breathing)

%% load periodical behavior data ("thetas" - variable frame rate ~500Hz)
[fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'Perodical Behavior Data','C:\Data\Ephys\Behav');
load([dirName fileName]);

%% get video sync data 
videoFrameTimes=readVideoTTLData;

periodBehavData=[thetas(videoFrameTimes.TTLFrames(1):size(videoFrameTimes.frameTime_ms,1))',... %Trace
    videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1):end)-...
    videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1))]; % Time points
% figure; hold on
% plot(periodBehavData(:,2),periodBehavData(:,1))

% [periodBehavData(:,1),periodBehavData(:,2)] = resample(periodBehavData(:,1),periodBehavData(:,2),'pchip'); 

periodBehavData=interp1(periodBehavData(:,2),periodBehavData(:,1),periodBehavData(1,2):periodBehavData(end,2));
% plot(interp)

%% get recording traces (loads rawData) (sampling rate usually at 30kHz)
[fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'Recording Data','C:\Data\Ephys\export');
load([dirName fileName]);

electrodeNum=10;


%% get spike and TTL times (loads spikeData and TTLs variables)
[fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'Get Spikes and TTL','C:\Data\Ephys\export');
cd(dirName); 

fileName=cell2mat(regexp(fileName,'.+(?=_)','match')); %'039v_0925_2Hz20ms_20mW_28Ch_nopp'; % '039v_0927_2Hz20ms_20mW_28Ch_nopp'; % 'SpVi12_133_2Hz2ms_10mW_nopp';
spikeData=load([fileName '_Ch' num2str(electrodeNum) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
load([fileName '_Ch' num2str(electrodeNum) '.mat'],'TTLs');
unique(spikeData.unitsIdx)
figure; hold on 
for cluster=1:3
    plot(mean(spikeData.waveForms(:,spikeData.unitsIdx==cluster),2))
end

figure; hold on 
plot(rawData(electrodeNum,:))
for cluster=2
    plot(spikeData.spikeTimes(spikeData.unitsIdx==cluster),...
    ones(size(spikeData.spikeTimes(spikeData.unitsIdx==cluster),1))*-1500,'r*')
end

%% downsample trace to 1ms
recordingTrace=rawData(electrodeNum,:);
recordingTrace=recordingTrace(TTLs.TTLtimes(1)*double(spikeData.samplingRate)/TTLs.samplingRate:end);
recordingTrace=decimate(double(recordingTrace),double(spikeData.samplingRate)/TTLs.samplingRate,'fir');

%% bin spike counts in 1ms bins
spikeTimes=double(spikeData.spikeTimes(spikeData.unitsIdx==2)); %figure; plot(spikeTimes)
spikeTimes=spikeTimes-(TTLs.TTLtimes(1)*double(spikeData.samplingRate)/TTLs.samplingRate);
spikeTimes=spikeTimes(spikeTimes>0);
dsSpikeTimes=(spikeTimes/(double(spikeData.samplingRate)/TTLs.samplingRate));

binSize=1;
binSpikeTime=DownSampleToMilliseconds(spikeTimes,binSize,spikeData.samplingRate);
% figure; hold on 
% bar(binSpikeTime)
% plot(dsSpikeTimes,ones(length(dsSpikeTimes),1),'r*')

figure; hold on
plot(recordingTrace)
plot(dsSpikeTimes,ones(length(dsSpikeTimes),1)*-250,'r*')
% plot(dsSpikeTimes,ones(length(dsSpikeTimes),1)*-220,'*')

%% compute sdf
sdf=GaussConv(binSpikeTime,5)*1000;
figure; hold on
plot(sdf)
plot(dsSpikeTimes,zeros(length(dsSpikeTimes),1),'r*')

%% make high-pass and low pass periodic behavior traces 
lp_periodBehavData=FilterTrace(periodBehavData,0.3,1000,'low');
figure; hold on
plot(periodBehavData); plot(lp_periodBehavData,'LineWidth',2)

bp_periodBehavData=FilterTrace(periodBehavData,[1 20],1000,'bandpass');
figure; hold on
plot(periodBehavData-mean(periodBehavData)); plot(bp_periodBehavData,'LineWidth',1)

% lp_periodBehavData=lp_periodBehavData(1:length(sdf));
% bp_periodBehavData=bp_periodBehavData(1:length(sdf));

%% plot traces together
figure; hold on
trace=periodBehavData-mean(periodBehavData); plot(trace); shift=ceil(max(trace));
trace=lp_periodBehavData-min(lp_periodBehavData); plot(trace+shift); shift=shift+ceil(max(trace));
trace=sdf-min(sdf); plot(trace+shift); shift=shift+ceil(max(trace));
trace=bp_periodBehavData*3-min(bp_periodBehavData); plot(trace+shift);

%% plot cross-correlation
[acor,lag] = xcorr(sdf,bp_periodBehavData,2500);
figure
plot(lag,acor)

[acor,lag] = xcorr(sdf,lp_periodBehavData,2500);
figure
plot(lag,acor)


