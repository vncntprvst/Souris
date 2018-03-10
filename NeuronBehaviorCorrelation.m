%% Correlation between bursts/spike rate and periodic behaviors (whisking, breathing)

% This script will need data from multiple source files:
% - recording trace (...raw.mat)
% - spikes
% - TTL times
% - Behavioral data
% - [optional] frame times (if irregular behavioral sampling rate) 
% Try to put them all in the same folder for easier processing

dirListing=dir;
dirName=cd; 
samplingRate=30000;
%%%%%%%%%%%%%%%%%%
%% neuronal data %
%%%%%%%%%%%%%%%%%%

%% get recording traces (loads allTraces) (sampling rate usually at 30kHz)
try
    fileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'traces.mat'),...
        {dirListing.name},'UniformOutput',false))).name;
catch
    [fileName,dirName] = uigetfile({'*.mat','.mat Files';...
        '*.*','All Files' },'Recording Data',cd);
end
load(fullfile(dirName,fileName));
allTraces=rawData;
channelNum=[26 30 31];
% 
% %% Load "original" traces (.continuous files from Open Ephys)
% fileName='100_CH19.continuous';
% channels=[18  19  20  23  26  28  30  31] -1;
% options.refCh=1;
% options.filter={[300 6000],'bandpass'};
% options.CAR=true;
% 
% [source_traces,source_timestamps,source_info]=Load_OE_Traces(fileName,channels,options); % .continuous files

%% get spike and TTL times (loads spikeData and TTLs variables)
try
    fileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'Ch.csv'),...
        {dirListing.name},'UniformOutput',false))).name;
catch
    [fileName,dirName] = uigetfile({'*.mat; *.csv','.mat and .csv Files';...
        '*.*','All Files' },'Get Spikes and TTL',cd);
    cd(dirName);
end
keepUnits=[1 2 3]; %[8 10 11 12 13 15];
titularChannels=[10 10 10]; %[10 14 14 15 15 15];
if contains (fileName,'.csv') % from JRClust
   allSpikeData=load(fileName);
   spikeData.samplingRate=samplingRate;
   spikeData.spikeTimes=allSpikeData(:,1)*spikeData.samplingRate;
   spikeData.unitsIdx=allSpikeData(:,2);
   spikeData.channelNum=allSpikeData(:,3);
elseif contains(fileName,'.mat')
    fileName=cell2mat(regexp(fileName,'.+(?=_)','match')); %'039v_0925_2Hz20ms_20mW_28Ch_nopp'; % '039v_0927_2Hz20ms_20mW_28Ch_nopp'; % 'SpVi12_133_2Hz2ms_10mW_nopp';
    spikeData=load([fileName '_Ch' num2str(channelNum(end)) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
    load([fileName '_Ch' num2str(channelNum(end)) '.mat'],'TTLs');
elseif  contains(fileName,'.spikes')
    source_spikeData=Load_OE_SpikesData(fileName); % .spikes file
end
% TTL might be empty - then load _trials
if ~exist('TTLs','var') || ~isfield(TTLs,'TTLtimes')
    try
        fileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'trials.mat'),...
            {dirListing.name},'UniformOutput',false))).name;
    catch
        [fileName,dirName] = uigetfile({'*.mat','.mat Files';...
            '*.*','All Files' },'load trials file',cd);
    end
    load(fullfile(dirName,fileName));
    TTLs.TTLtimes=trials.TTL_times-trials.startClockTime;
    TTLs.samplingRate=trials.samplingRate{1};
end

% also reset spike times if needed
if spikeData.spikeTimes(end) > size(allTraces,2)
    spikeData.spikeTimes=spikeData.spikeTimes-uint32(trials.startClockTime);
end

% unique(spikeData.unitsIdx)
% figure; hold on 
% for cluster=1:3
%     plot(mean(spikeData.waveForms(:,spikeData.unitsIdx==cluster),2))
% end

%% Keep selected recording trace and spike times, set to initial TTL 
%first, traces
keepTraces=10; %14; %[10 14 15]; 
recordingTrace=cell(length(keepTraces),1);
for traceNum=1:length(keepTraces)
    recordingTrace{traceNum}=allTraces(keepTraces(traceNum),:); %select the trace to keep
    recordingTrace{traceNum}=recordingTrace{traceNum}(TTLs.TTLtimes(1)*double(spikeData.samplingRate)/...
        double(TTLs.samplingRate):end); % cut out trace that occurs before first TTL
end

% same for spikes from selected units
spikeTimes=cell(length(keepUnits),1);
for clusterNum=1:length(keepUnits)
    spikeTimes{clusterNum}=double(spikeData.spikeTimes(spikeData.unitsIdx==keepUnits(clusterNum))); %figure; plot(spikeTimes)
    spikeTimes{clusterNum}=spikeTimes{clusterNum}-(double(TTLs.TTLtimes(1))*double(spikeData.samplingRate)/double(TTLs.samplingRate));
    spikeTimes{clusterNum}=spikeTimes{clusterNum}(spikeTimes{clusterNum}>0);
end
% rs_spikeTimes=double(source_spikeData.spikeTimes(source_spikeData.unitsIdx==2)); %figure; plot(spikeTimes)

figure('Color','white'); 
for traceNum=1:length(keepTraces)
    subplot(length(keepTraces),1,traceNum); hold on;
    plot(recordingTrace{traceNum});
    correspondingUnits=find(ismember(titularChannels,keepTraces(traceNum)));
    for clusterNum=1:length(correspondingUnits)
        plot(spikeTimes{correspondingUnits(clusterNum)},...
            ones(size(spikeTimes{correspondingUnits(clusterNum)},1),1)*-300,'*');
        %     plot(rs_spikeTimes,ones(size(rs_spikeTimes,1),1)*-300,'d')
    end
    midRec=round(length(recordingTrace{traceNum})/2);
    set(gca,'ylim',[-500 500],'xlim',[midRec-samplingRate midRec+samplingRate]);
end


% params.Fs=30000;params.fpass=[0 25];params.tapers=[2 3];params.pad=1;params.err=[2 0.05];params.trialave=0;
% [C,phi,S12,S1,S2,f]=coherencyc(foo',double(allTraces(channelNum,:))',params);
% figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

%% downsample trace to 1ms (pointless, as it will loose the spikes)
% Fs=1000;
% % foo=resample(double(recordingTrace),Fs,double(spikeData.samplingRate));
% dsrecordingTrace=decimate(double(recordingTrace),double(spikeData.samplingRate)/Fs,'fir');

%% bin spike counts in 1ms bins
% with Chronux' binning function
foo=binspikes(spikeTimes/double(spikeData.samplingRate),Fs);
foo=[zeros(round(spikeTimes(1)/double(spikeData.samplingRate)*Fs)-1,1);foo]; %need to padd with zeroes
% with home-made function. Same result, but takes care of the padding
binSize=1;
binSpikeTimes=cell(length(keepUnits),1);
for clusterNum=1:length(keepUnits)
    binSpikeTimes{clusterNum,1}=DownSampleToMilliseconds(spikeTimes{clusterNum},binSize,spikeData.samplingRate);
end

figure; hold on
% plot(dsrecordingTrace)
plot(find(binSpikeTimes),ones(length(find(binSpikeTimes)),1)*-250,'r*')
plot(find(foo),ones(length(find(foo)),1)*-200,'g*')

%% compute sdfs
SDFs=cell(length(keepUnits),1);
for clusterNum=1:length(keepUnits)
    SDFs{clusterNum}=GaussConv(binSpikeTimes{clusterNum},5)*1000;
end
figure; hold on
plot(SDFs{1})
plot(find(binSpikeTimes{1}),ones(length(find(binSpikeTimes{1})),1)*-10,'r*')

%% compute rasters
[indy, indx]=deal(cell(length(keepUnits),1));
for clusterNum=1:length(keepUnits)
    [indy{clusterNum}, indx{clusterNum}] = ind2sub(size(binSpikeTimes{clusterNum}),find(binSpikeTimes{clusterNum})); %find row and column coordinates of spikes
end
% rasters=[indx indy;indx indy+1];

%%%%%%%%%%%%%%%%%%%%
%% behavioral data %
%%%%%%%%%%%%%%%%%%%%

%% load periodical behavior data ("thetas" - variable frame rate ~500Hz)
% e.g. WhiskerAngle.csv
try
    fileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_WhiskerAngle.csv'),...
        {dirListing.name},'UniformOutput',false))).name;
catch
[fileName,dirName] = uigetfile({'*.mat; *.csv; *.avi','.mat Files; .csv Files; video Files';...
    '*.*','All Files'},'Whisker angle / Perodical Behavior Data / Video file','C:\Data\Ephys\Behav');
end
if contains(fileName,'.csv') %not yet converted to mat file
    % import data
    thetas=ImportCSVasVector(fullfile(dirName,fileName));
    thetas=thetas*180/pi; %convert to degrees.
    % smooth data to remove outliers - also takes care of missing data points
    thetas=smoothdata(thetas,'rloess',9);
     %fill missing / NaNs values (if any)
    thetas = fillmissing(thetas,'spline'); %linear %spline %movmedian

    [b,a] = butter(3,40/250,'low');
    foo = filtfilt(b,a,thetas);
    
elseif contains(fileName,'.avi') %video file to extract whisker angle
    thetas=ExtractMultiWhiskerAngle_FFTonContours(fullfile(dirName,fileName));
    thetas=smoothdata(thetas,'rloess',20);
else
    load([dirName fileName]);
end
%% get video sync data 
videoFrameTimes=ReadVideoFrameTimes;

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

%% filter periodic behavior traces into low-pass and bandpassed versions
LP_periodBehavData=FilterTrace(periodBehavData,0.3,1000,'low'); %set-point
figure; hold on
plot(periodBehavData); plot(LP_periodBehavData,'LineWidth',2)

BP_periodBehavData=FilterTrace(periodBehavData,[0.3 30],1000,'bandpass'); %whisking
figure; hold on
plot(periodBehavData-mean(periodBehavData)); plot(BP_periodBehavData,'LineWidth',1)

HP_periodBehavData=FilterTrace(periodBehavData,0.3,1000,'high')'; %whisking
plot(HP_periodBehavData,'LineWidth',1)

%% cut down behavior and spike traces to same length
shortestDuration=min(cellfun(@(x) length(x), SDFs));
SDFs=cellfun(@(x) x(1:shortestDuration), SDFs,'UniformOutput',false);
LP_periodBehavData=LP_periodBehavData(1:shortestDuration);
BP_periodBehavData=BP_periodBehavData(1:shortestDuration);
HP_periodBehavData=HP_periodBehavData(1:shortestDuration);
periodBehavData=periodBehavData(1:shortestDuration);

if exist('whiskingPeriod','var')
%     sdf=sdf(whiskingPeriod);
%     LP_periodBehavData=LP_periodBehavData(whiskingPeriod);
%     BP_periodBehavData=BP_periodBehavData(whiskingPeriod);
end

% figure; hold on
% plot(SDFs{1}); plot(BP_periodBehavData/min(BP_periodBehavData)*...
%     max(SDFs{1}) + max(SDFs{1}))
 
%% Hilbert transform
HTBP_periodBehavData=hilbert(BP_periodBehavData);
% figure; plot(imag(HTBP_periodBehavData));
whiskingPhase=angle(HTBP_periodBehavData);
% figure; hold on
% plot(SDFs); plot(whiskingPhase*10 + max(SDFs))

%% find instantaneous frequency
Nfft = 1024;
% [Pxx,f] = pwelch(BP_periodBehavData,gausswin(Nfft),Nfft/2,Nfft,1000);
% figure; plot(f,Pxx); ylabel('PSD'); xlabel('Frequency (Hz)'); grid on;
[~,sgFreq,sgTime,sgPower] = spectrogram(BP_periodBehavData,gausswin(Nfft),Nfft/2,Nfft,1);
instantFreq = medfreq(sgPower,sgFreq);
% figure; plot(sgTime,round(instantFreq*1000),'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%
%% plot traces together
%%%%%%%%%%%%%%%%%%%%%%%
displayUnits=1; %4; %1:length(keepUnits); %default: 1:length(keepUnits) [1 4 5]
% plot 4 seconds around a designated time point
minFreqIdx=sgTime(instantFreq<rms(instantFreq)); % low oscillation frequency
whiskingAmp=abs(BP_periodBehavData-LP_periodBehavData);% whisking amplitude
thAmpIdx=minFreqIdx(whiskingAmp(minFreqIdx)==max(whiskingAmp(minFreqIdx))); 
dispWin=[thAmpIdx-2000 thAmpIdx+1999]; timeAxis=dispWin(1):dispWin(2);
cmap=colormap(lines);
figure('Color','white'); hold on;
% axes('Position',[0 1-(1/(4)*(4)) 1 1/(4)]); 
% subplot(2,1,2); hold on; axis(gca,'off');
for clusterNum=1:length(displayUnits) %:-1:1
    % subplot('Position',[0 1-(1/(4)*(clusterNum+3)) 1 1/(4)]);
    winIndx=indx{displayUnits(clusterNum)}...
        (indx{displayUnits(clusterNum)}>=dispWin(1) & indx{displayUnits(clusterNum)}<=dispWin(end))-timeAxis(1);
    winIndy=indy{displayUnits(clusterNum)}...
        (indx{displayUnits(clusterNum)}>=dispWin(1) & indx{displayUnits(clusterNum)}<=dispWin(end));
    plot([winIndx;winIndx],[winIndy+clusterNum-1;winIndy+clusterNum],...
        'color',cmap(clusterNum,:),'LineStyle','-','LineWidth',1.8); % plot rasters% plot rasters
    plot(BP_periodBehavData(timeAxis)/max(abs(BP_periodBehavData(timeAxis)))+...
        clusterNum+0.5,'color','k','LineWidth',0.8);
    plot(whiskingPhase(timeAxis)/max(abs(whiskingPhase(timeAxis)))+...
        clusterNum+0.5,'color','k','LineWidth',0.8);
end
% plot(diff(BP_periodBehavData(timeAxis))/max(abs(diff(BP_periodBehavData(timeAxis))))+length(keepUnits)+3,'color','k','LineWidth',1.8);axis(gca, 'off');
% axis(gca, 'off');
set(gca,'xlim',[1 length(timeAxis)]);
% set(gca,'ytick',[0 80 150 300 450],'yticklabels',...
%     {'Raw whisking', 'Whisking set point', 'Whisking cycles', 'Spike density function', 'Rasters'});
% set(gca,'xtick',[0 50000 100000 150000 200000],'xticklabels',round([0 50000 100000 150000 200000]/60000,1));
xlabel('Time (ms)')
% set(gca,'xlim',);

%% plot cross-correlation
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs{displayUnits(clusterNum)},BP_periodBehavData,500,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor,'color','k','LineWidth',2); xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
    title({['Cross correlation for vIRt unit ' num2str(keepUnits(displayUnits(clusterNum)))];'Spike density function vs. Whisking angle'})
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs{displayUnits(clusterNum)},LP_periodBehavData,1000,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Set point'})
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs{displayUnits(clusterNum)},whiskingPhase,100,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Whisking phase'})
    % hold on; [acor,lag] = xcorr(sdf,BP_periodBehavData,100,'unbiased'); plot(lag,acor);
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs{displayUnits(clusterNum)},periodBehavData,2500,'coeff');
    figure('position',[602   537   560   420]); 
    plot(lag,acor); xlabel('Lag (ms)');set(gca,'ylim',[0 1]);
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Whisking angle'})
end
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

for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpt(BP_periodBehavData',...
        spikeTimes{displayUnits(clusterNum)}/double(spikeData.samplingRate),params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(BP_periodBehavData',...
        SDFs{displayUnits(clusterNum)}',params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskingPhase',...
        SDFs{displayUnits(clusterNum)}',params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
C(S1==max(S1))
end

% for smoother result, re-sample each with 16(?) ms bins
    % 
% test noise in CC by bootstrap
% cut signal randomly into trials 
% test CC for each trial
% take mean CC
% then randomize trials -> flat CC 
% if mean CC is > 2SD of noise (flat CC) 


% for manual coherence
% fft each signal 
% product(fft)/ sqr(auto-corr)


% for phase 


%% R01 figure
traceNum=1; %3;
vIRtTrace=recordingTrace{traceNum};
midRec=round(length(vIRtTrace)/2)/30000;
vibrissaAngle=resample(HP_periodBehavData,30,1); %foo=foo(1:length(recordingTrace{traceNum}));

figure('Color','white');
subplot(2,1,1); hold on
timeAxis=(1:length(vIRtTrace))/30000;
plot(timeAxis,vIRtTrace,'color','k','LineWidth',0.8);
correspondingUnits=find(ismember(titularChannels,keepTraces(traceNum)));
for clusterNum=1:length(correspondingUnits)
    plot(spikeTimes{correspondingUnits(clusterNum)}/30000,...
        ones(size(spikeTimes{correspondingUnits(clusterNum)},1),1)*-1500+clusterNum*100,'d');
end

set(gca,'ylim',[-1500 1000],'xlim',[661 665],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out');
xlabel('Time (s)'); ylabel('Firing rate (Hz)');

subplot(2,1,2); hold on
timeAxis=(1:length(vibrissaAngle))/30000;
plot(timeAxis,vibrissaAngle,'color','k','LineWidth',0.8);
for clusterNum=1:length(correspondingUnits)
    plot(spikeTimes{correspondingUnits(clusterNum)}/30000,...
        ones(size(spikeTimes{correspondingUnits(clusterNum)},1),1)+clusterNum,'d');
end
set(gca,'xlim',[661 665],'ylim',[-15 15],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out');% axis 'tight'
xlabel('Time (s)');  ylabel('Angle (\circ)');

% cross correlation
bestUnit=2; %4;
[acor,lag] = xcorr(SDFs{bestUnit},HP_periodBehavData,500,'coeff');
figure('position',[602   537   560   420]);
plot(lag,acor,'color','k','LineWidth',2); xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(bestUnit))];'Spike density function vs. Whisking angle'})

%% polar plot
%
% unitSDF=SDFs{bestUnit};
% unitSDFExcerpt=unitSDF(int32([midRec-100:midRec+100]*1000));
% get atda
unitSpikes=binSpikeTimes{bestUnit};
vibrissaAngle=BP_periodBehavData; 
unitSpikes=unitSpikes(1:length(vibrissaAngle));
numSegments=round(length(vibrissaAngle)/60000);

for pplotNum=1:numSegments-1
    %excerpt
    timeIndex=(pplotNum-1)*60+1;
    timeWindowIdx=int32(linspace(timeIndex*1000,(timeIndex+60)*1000-1,60*1000)); %int32([midRec-300:midRec+100]*1000);
    vibrissaAngleExcerpt=vibrissaAngle(timeWindowIdx);
    unitSpikesExcerpt=unitSpikes(timeWindowIdx);
    
    % coordConversion=90; %adjust depending on camera position
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt+coordConversion; % *180/pi;
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt/180*pi; %convert back to radians
    
    % Hilbert transform NEEDS TO BE ZERO CENTRED !!!
    HTvibrissaAngle=hilbert(vibrissaAngleExcerpt);
    % figure; plot(whiskingPhase);
    whiskingPhase=angle(HTvibrissaAngle);
    
    spikeOnWPhase=whiskingPhase(logical(unitSpikesExcerpt));
    figure; polarhistogram(spikeOnWPhase,72,'Displaystyle','stairs',...
        'Normalization','count')
end
% [whiskingPhase,sortIdx]=sort(whiskingPhase);
% unitSDFExcerpt=unitSDFExcerpt(sortIdx);
% 
% whiskingPhase=whiskingPhase(unitSDFExcerpt>5);
% unitSDFExcerpt=unitSDFExcerpt(unitSDFExcerpt>5);

% figure;polarplot(whiskingPhase,unitSDFExcerpt)


% timeAxis=(1:length(unitSpikes))*30000;
% figure; hold on 
% plot(timeAxis,vibrissaAngle)
% plot(find(unitSpikes)*30000,zeros(length(find(unitSpikes)),1),'d')
% plot(timeAxis,whiskingPhase)




