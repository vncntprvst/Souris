%% Correlation between bursts/spike rate and periodical behaviors (whisking, breathing)

%% load periodical behavior data ("thetas" - variable frame rate ~500Hz)
[fileName,dirName] = uigetfile({'*.mat; *.csv; *.avi','.mat Files; .csv Files; video Files';...
    '*.*','All Files'},'Whisker angle / Perodical Behavior Data / Video file','C:\Data\Ephys\Behav');
if contains(fileName,'.csv') %not yet converted to mat file
    % import data
    thetas=ImportCSVasVector(fullfile(dirName,fileName));
    % smooth data to remove outliers - also takes care of missing data points
    thetas=smoothdata(thetas,'rloess',9);
     %fill missing / NaNs values (if any)
    thetas = fillmissing(thetas,'spline'); %linear %spline %movmedian

    [b,a] = butter(3,20/250,'low');
    foo = filtfilt(b,a,thetas);
    
elseif contains(fileName,'.avi') %video file to extract whisker angle
    thetas=ExtractMultiWhiskerAngle_FFTonContours(fullfile(dirName,fileName));
    thetas=smoothdata(thetas,'rloess',20);
else
    load([dirName fileName]);
end
%% get video sync data 
videoFrameTimes=readVideoTTLData;

%% create array with angle values and time points
periodBehavData=[thetas(videoFrameTimes.TTLFrames(1):size(videoFrameTimes.frameTime_ms,1)),... %Trace
    videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1):end)-...
    videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1))+1]; % Time points
% figure; hold on
% plot(periodBehavData(:,2),periodBehavData(:,1))

%resample to 1ms precision
% [periodBehavData(:,1),periodBehavData(:,2)] = resample(periodBehavData(:,1),periodBehavData(:,2),'pchip'); 
periodBehavData=interp1(periodBehavData(:,2),periodBehavData(:,1),periodBehavData(1,2):periodBehavData(end,2));


%% plot behavior data and find a period with whisking bouts
% no need to keep periods with no whisking
figure; plot(periodBehavData); % select data point and export cursor info
whiskingPeriod=1:cursor_info.Position(1); %in ms

%% get recording traces (loads rawData) (sampling rate usually at 30kHz)
[fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'Recording Data',cd);
load([dirName fileName]);

channelNum=19;

% original traces:
fileName='100_CH19.continuous';
channels=[18  19  20  23  26  28  30  31] -1;
options.refCh=1;
options.filter={[300 6000],'bandpass'};
options.CAR=true;

[or_traces,or_timestamps,or_info]=Load_OE_Traces(fileName,channels,options);

%% get spike and TTL times (loads spikeData and TTLs variables)
[fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'Get Spikes and TTL',cd);
cd(dirName); 
if contains(fileName,'.mat')
    fileName=cell2mat(regexp(fileName,'.+(?=_)','match')); %'039v_0925_2Hz20ms_20mW_28Ch_nopp'; % '039v_0927_2Hz20ms_20mW_28Ch_nopp'; % 'SpVi12_133_2Hz2ms_10mW_nopp';
    spikeData=load([fileName '_Ch' num2str(channelNum) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
    load([fileName '_Ch' num2str(channelNum) '.mat'],'TTLs');
elseif  contains(fileName,'.spikes')
    resorted_spikeData=Load_OE_SpikesData(fileName);
end
% TTL might be empty - then load _trials
if ~exist('TTLs','var') || ~isfield(TTLs,'TTLtimes')
    [fileName,dirName] = uigetfile({'*.mat','.mat Files';...
    '*.*','All Files' },'load trials file',cd);
    load(fullfile(dirName,fileName));
    TTLs.TTLtimes=trials.TTL_times-trials.startClockTime;
    TTLs.samplingRate=trials.samplingRate{1};
end

% also reset spike times if needed
if spikeData.spikeTimes(end) > size(rawData,2)
    spikeData.spikeTimes=spikeData.spikeTimes-uint32(trials.startClockTime);
end

unique(spikeData.unitsIdx)
figure; hold on 
for cluster=1:3
    plot(mean(spikeData.waveForms(:,spikeData.unitsIdx==cluster),2))
end

%% Keep one recording trace and spike times, set to initial TTL 
channelNum=2; 
recordingTrace=rawData(channelNum,:); %select the trace to keep
recordingTrace=recordingTrace(TTLs.TTLtimes(1)*double(spikeData.samplingRate)/...
    double(TTLs.samplingRate):end); % cut out trace that occurs before first TTL

% same for spikes from selected unit
clusterNum=2;
spikeTimes=double(spikeData.spikeTimes(spikeData.unitsIdx==clusterNum)); %figure; plot(spikeTimes)
spikeTimes=spikeTimes-(double(TTLs.TTLtimes(1))*double(spikeData.samplingRate)/double(TTLs.samplingRate));
spikeTimes=spikeTimes(spikeTimes>0);

rs_spikeTimes=double(resorted_spikeData.spikeTimes(resorted_spikeData.unitsIdx==2)); %figure; plot(spikeTimes)

figure; hold on 
plot(recordingTrace)
for cluster=clusterNum
    plot(spikeTimes,ones(size(spikeTimes,1),1)*-500,'*')
    plot(rs_spikeTimes,ones(size(rs_spikeTimes,1),1)*-300,'d')
end

foo=resample(periodBehavData,30,1); foo=foo(1:length(rawData(channelNum,:)));
plot(foo*1000-200)

% params.Fs=30000;params.fpass=[0 25];params.tapers=[2 3];params.pad=1;params.err=[2 0.05];params.trialave=0;
% [C,phi,S12,S1,S2,f]=coherencyc(foo',double(rawData(channelNum,:))',params);
% figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

%% downsample trace to 1ms (pointless, as it will loose the spikes)
Fs=1000;
% foo=resample(double(recordingTrace),Fs,double(spikeData.samplingRate));
dsrecordingTrace=decimate(double(recordingTrace),double(spikeData.samplingRate)/Fs,'fir');

%% bin spike counts in 1ms bins
% with Chronux' binning function
foo=binspikes(spikeTimes/double(spikeData.samplingRate),Fs);
foo=[zeros(round(spikeTimes(1)/double(spikeData.samplingRate)*Fs)-1,1);foo]; %need to padd with zeroes
% with home-made function. Same result, but takes care of the padding
binSize=1;
binSpikeTime=DownSampleToMilliseconds(spikeTimes,binSize,spikeData.samplingRate);

figure; hold on
plot(dsrecordingTrace)
plot(find(binSpikeTime),ones(length(find(binSpikeTime)),1)*-250,'r*')
plot(find(foo),ones(length(find(foo)),1)*-200,'g*')

%% compute sdf
sdf=GaussConv(binSpikeTime,5)*1000;
figure; hold on
plot(sdf)
plot(find(binSpikeTime),ones(length(find(binSpikeTime)),1)*-10,'r*')

%% compute rasters
[indy, indx] = ind2sub(size(binSpikeTime),find(binSpikeTime)); %find row and column coordinates of spikes
% rasters=[indx indy;indx indy+1];

%% filter periodic behavior traces into low-pass and bandpassed versions
LP_periodBehavData=FilterTrace(periodBehavData,0.3,1000,'low'); %set-point
figure; hold on
plot(periodBehavData); plot(LP_periodBehavData,'LineWidth',2)

BP_periodBehavData=FilterTrace(periodBehavData,[0.3 25],1000,'bandpass'); %whisking
figure; hold on
plot(periodBehavData-mean(periodBehavData)); plot(BP_periodBehavData,'LineWidth',1)

%% cut down behavior and spike traces to same length
LP_periodBehavData=LP_periodBehavData(1:length(sdf));
BP_periodBehavData=BP_periodBehavData(1:length(sdf));

if exist('whiskingPeriod','var')
%     sdf=sdf(whiskingPeriod);
%     LP_periodBehavData=LP_periodBehavData(whiskingPeriod);
%     BP_periodBehavData=BP_periodBehavData(whiskingPeriod);
end

figure; hold on
plot(sdf); plot(BP_periodBehavData/min(BP_periodBehavData)*max(sdf) + max(sdf))

%% Hilbert transform
HTBP_periodBehavData=hilbert(BP_periodBehavData);
% figure; plot(imag(HTBP_periodBehavData));
whiskingPhase=angle(HTBP_periodBehavData);
figure; hold on
plot(sdf); plot(whiskingPhase*10 + max(sdf))

%% plot traces together
figure; hold on
trace=periodBehavData-mean(periodBehavData); plot(trace); shift(1)=ceil(max(trace));
trace=LP_periodBehavData-min(LP_periodBehavData); plot(trace+shift(1)); shift(2)=shift(1)+ceil(max(trace));
trace=BP_periodBehavData*3-min(BP_periodBehavData); plot(trace+shift(2)); shift(3)=shift(2)+ceil(max(trace));
trace=sdf-min(sdf); plot(trace+shift(3)); shift(4)=shift(3)+ceil(max(trace));
plot([indx;indx],[indy+shift(4);indy+shift(4)+50],'color','k'); shift(1)=shift(1)+51;% plot rasters
set(gca,'ytick',[0 80 150 300 450],'yticklabels',...
    {'Raw whisking', 'Whisking set point', 'Whisking cycles', 'Spike density function', 'Rasters'});
set(gca,'xtick',[0 50000 100000 150000 200000],'xticklabels',round([0 50000 100000 150000 200000]/60000,1));
xlabel('Time (mn)')
% set(gca,'xlim',);
axis(gca, 'off'); 

%% plot cross-correlation
[acor,lag] = xcorr(sdf,BP_periodBehavData,100,'coeff');
figure
plot(lag,acor); xlabel('Lag (ms)');set(gca,'ylim',[-1 1])
title({'Cross correlation';'Spike density function vs. Whisking angle'})

[acor,lag] = xcorr(sdf,LP_periodBehavData,1000,'coeff');    
figure
plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
title({'Cross correlation';'Spike density function vs. Set point'})

[acor,lag] = xcorr(sdf,whiskingPhase,100,'coeff');
figure
plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
title({'Cross correlation';'Spike density function vs. Whisking phase'})
% hold on; [acor,lag] = xcorr(sdf,BP_periodBehavData,100,'unbiased'); plot(lag,acor);

[acor,lag] = xcorr(sdf,periodBehavData,2500);
figure
plot(lag,acor); xlabel('Lag (ms)')

% shuffle-subtracted correlogram
% shuffle "trials"
%% Phase / coherence analysis
% Make polar plot of the coherence between spiking activity and
% vibrissa motion at the peak frequency of whisking.

% The relationship between unit activities and the whisking or 
% breathing rhythms was estimated by computing the spectral coherence
% (Mitra and Pesaran, 1999). Single whisks or breaths were isolated by
% band-pass filtering the position traces between 1 Hz and 25 Hz with
% a three-pole Butterworth filter run in forward and backward directions,
% and applying the Hilbert transform (Hill et al., 2008).
% The Chronux toolbox (http://chronux.org) was used to compute coherence
% between spike times or intracellular events averaged over 2 s segments
% with a time-bandwidth product of two. Basal respiration was defined as
% the instantaneous respiratory frequency below 4 Hz and sniffing as the
% instantaneous respiratory frequency greater than 4 Hz. We report the
% magnitude and phase of the coherence at the peak frequency of whisking
% or breathing. A phase of zero corresponds to the peak of protraction
% or inspiration.

params.Fs=1000; % sampling frequency
params.fpass=[0 25]; % band of frequencies to be kept
params.tapers=[2 3]; % taper parameters
params.pad=1; % pad factor for fft
params.err=[2 0.05];
params.trialave=0;

[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpt(BP_periodBehavData',spikeTimes/double(spikeData.samplingRate),params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

C(S1==max(S1))

[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(BP_periodBehavData',sdf',params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

C(S1==max(S1))

[C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(whiskingPhase',sdf',params);
figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

C(S1==max(S1))


