% NeuronBehaviorCorrelation_Plots

%% First, get data
[recordingTraces,spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
    SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
    BP_periodBehavData_ms,HP_periodBehavData_ms,LP_periodBehavData_ms,...
    HTBP_periodBehavData_ms, peakWhisking_ms,periodBehavData_ms,...
    whiskingPhase_ms,instantFreq_ms,sgFreq,sgTime,sgPower] = NeuronBehaviorCorrelation_GatherData;

keepTraces=unique(preferredElectrode(ismember(unitID,keepUnits)));
peakWhiskingIdx=find(peakWhisking_ms==max(peakWhisking_ms));
whiskingPeriod=peakWhiskingIdx-5000:peakWhiskingIdx+4999; %in ms
whiskingPeriodIdx=[0 abs(diff(BP_periodBehavData_ms))>=mad(diff(BP_periodBehavData_ms))] &...
                    abs(BP_periodBehavData_ms)>=mad(BP_periodBehavData_ms);
whiskingPeriodIdx=movsum(whiskingPeriodIdx,500)>0;

%% Polar plot by time chunks
% unitSDF=SDFs{bestUnit};
% unitSDFExcerpt=unitSDF(int32([midRec-100:midRec+100]*1000));
% get atda
vibrissaAngle=BP_periodBehavData_ms;
for unitNum=1:size(spikeRasters_ms,1)
unitSpikes=spikeRasters_ms(unitNum,:);
% unitSpikes=unitSpikes(1:length(vibrissaAngle));
numSegments=round(length(vibrissaAngle)/60000);
figure;
for pplotNum=1:numSegments-2
    %excerpt
    timeIndex=(pplotNum-1)*60+1;
    timeWindowIdx=int32(linspace(timeIndex*1000,(timeIndex+60)*1000-1,60*1000)); %int32([midRec-300:midRec+100]*1000);
%     vibrissaAngleExcerpt=vibrissaAngle(timeWindowIdx);
    unitSpikesExcerpt=unitSpikes(timeWindowIdx);
    
    % coordConversion=90; %adjust depending on camera position
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt+coordConversion; % *180/pi;
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt/180*pi; %convert back to radians
    
    % Hilbert transform NEEDS TO BE ZERO CENTRED !!!
%     HTvibrissaAngle=hilbert(vibrissaAngleExcerpt);
    HTvibrissaAngleExcerpt=HTBP_periodBehavData_ms(timeWindowIdx);
    % figure; plot(whiskingPhase);
    whiskingPhaseExcerpt_ms=angle(HTvibrissaAngleExcerpt);
    
    spikeOnWPhase=whiskingPhaseExcerpt_ms(logical(unitSpikesExcerpt));
    subplot(3,2,pplotNum);
    polarhistogram(spikeOnWPhase,72,'Displaystyle','stairs',...
        'Normalization','count')
end
end

%% Polar plot on single time window
bestWhiskingPeriod=5;
timeIndex=(bestWhiskingPeriod-1)*60+1;
% whiskingPeriod=int32(linspace(timeIndex*1000,(timeIndex+60)*1000-1,60*1000));
whiskingPeriod=whiskingPeriodIdx;
for unitNum=5 %1:size(spikeRasters_ms,1)
    unitSpikes=spikeRasters_ms(unitNum,:);
    unitSpikesExcerpt=unitSpikes(whiskingPeriod);
%     vibrissaAngle=BP_periodBehavData_ms;   
%     vibrissaAngleExcerpt=vibrissaAngle(whiskingPeriod);
    % coordConversion=90; %adjust depending on camera position
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt+coordConversion; % *180/pi;
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt/180*pi; %convert back to radians
    
    % Hilbert transform NEEDS ANGLE TO BE ZERO CENTERED !!!
%     HTvibrissaAngleExcerpt=HTBP_periodBehavData_ms(whiskingPeriod);
    % figure; plot(whiskingPhase);
%     whiskingPhaseExcerpt_ms=angle(HTvibrissaAngleExcerpt);
    whiskingPhaseExcerpt_ms=whiskingPhase_ms(whiskingPeriod);
    spikeOnWPhase=whiskingPhaseExcerpt_ms(logical(unitSpikesExcerpt));
    
    figure;
    polarhistogram(spikeOnWPhase,72,'Displaystyle','stairs',...
        'Normalization','count')
end
unitsOfInterestIdx=[5,1,6,10];unitsOfInterest=keepUnits(unitsOfInterestIdx);

%% plot spike times with whisking angle
% angle and phase
timeVector=1:numel(periodBehavData_ms); timeVector(~whiskingPeriodIdx)=NaN;
figure; hold on 
% plot(timeVector,periodBehavData_ms-median(periodBehavData_ms))
plot(timeVector,BP_periodBehavData_ms); %-median(BP_periodBehavData_ms));
plot(timeVector,whiskingPhase_ms)

for unitNum=5 %1:size(spikeRasters_ms,1)
    unitSpikes=spikeRasters_ms(unitNum,:);
%     unitSpikesExcerpt=unitSpikes(whiskingPeriod);
    unitSpikes(isnan(timeVector) | unitSpikes==0)=nan;
    unitSpikes(unitSpikes==1)=unitNum;
    plot(timeVector,unitSpikes,'d')
end


%% Plot traces together
traceId=mode(preferredElectrode(unitID==unitsOfInterest(1)));
traceNum=find(keepTraces==traceId); %3;
vIRtTrace=recordingTraces(traceNum,:);
% midRec=round(length(vIRtTrace)/2)/30000;
vibrissaAngle=resample(HP_periodBehavData_ms,30,1); %upsampling back to 30kHz

figure('Color','white');
subplot(2,1,1); hold on
timeAxis=double(whiskingPeriod)*samplingRate/1000; %(1:length(vIRtTrace))/samplingRate;
plot(timeAxis-timeAxis(1)+1,vIRtTrace(timeAxis),'color','k','LineWidth',0.8);
%     correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
%     for clusterNum=1:numel(correspondingUnits)
%         plot(spikeTimes(unitID==correspondingUnits(clusterNum)),...
% correspondingUnits=find(ismember(titularChannels,keepTraces(traceNum)));
% correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
for clusterNum=1%:length(unitsOfInterest)
    unitSpikeTimes=spikeTimes(unitID==unitsOfInterest(clusterNum)); %/samplingRate;
    unitTWSpikeTimes=unitSpikeTimes(unitSpikeTimes>=timeAxis(1) &...
                                    unitSpikeTimes<=timeAxis(end))-timeAxis(1)+1;
    plot(unitTWSpikeTimes,ones(numel(unitTWSpikeTimes),1)*-500,'d'); %+clusterNum*100
end

set(gca,'ylim',[-550 500],'xlim',[1480000 1580000],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out'); %
xlabel('Time (s)'); ylabel('Firing rate (Hz)');

subplot(2,1,2); hold on
% timeAxis=(1:length(vibrissaAngle))/30000;
plot(timeAxis-timeAxis(1)+1,vibrissaAngle(timeAxis),'color','k','LineWidth',0.8);
for clusterNum=1:length(correspondingUnits)
    plot(spikeTimes(correspondingUnits(clusterNum),:)/samplingRate,...
        ones(size(spikeTimes(correspondingUnits(clusterNum),:),1),1)+clusterNum,'d');
end
for clusterNum=4%:length(unitsOfInterest)
    unitSpikeTimes=spikeTimes(unitID==unitsOfInterest(clusterNum)); %/samplingRate;
    unitTWSpikeTimes=unitSpikeTimes(unitSpikeTimes>=timeAxis(1) &...
                                    unitSpikeTimes<=timeAxis(end))-timeAxis(1)+1;
    plot(unitTWSpikeTimes,zeros(numel(unitTWSpikeTimes),1),'d'); %+clusterNum*100
end
set(gca,'xlim',[1480000 1580000],'ylim',[-10 25],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out');% axis 'tight'
xlabel('Time (s)');  ylabel('Angle (\circ)');

%% cross correlation
bestUnit=2; %4;
[acor,lag] = xcorr(SDFs_ms{bestUnit},HP_periodBehavData_ms,500,'coeff');
figure('position',[602   537   560   420]);
plot(lag,acor,'color','k','LineWidth',2); xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(bestUnit))];'Spike density function vs. Whisking angle'})

%% Plot traces together
%%%%%%%%%%%%%%%%%%%%%%%
displayUnits=1:size(spikeRasters_ms,1);  %4; %1:length(keepUnits); %default: 1:length(keepUnits) [1 4 5]
% plot 4 seconds around a designated time point
minFreqIdx=sgTime(instantFreq_ms<rms(instantFreq_ms)); % low oscillation frequency
whiskingAmp=abs(BP_periodBehavData_ms-LP_periodBehavData_ms);% whisking amplitude
thAmpIdx=minFreqIdx(whiskingAmp(minFreqIdx)==max(whiskingAmp(minFreqIdx)));
dispWin=[thAmpIdx-2000 thAmpIdx+1999]; timeAxis=dispWin(1):dispWin(2);
cmap=colormap(lines);
figure('Color','white'); hold on;
% axes('Position',[0 1-(1/(4)*(4)) 1 1/(4)]);
% subplot(2,1,2); hold on; axis(gca,'off');
for clusterNum=1:length(displayUnits) %:-1:1
    subplot('Position',[0 1-(1/(4)*(clusterNum+3)) 1 1/(4)]);
    winIndx=rasterXInd_ms{displayUnits(clusterNum)}...
        (rasterXInd_ms{displayUnits(clusterNum)}>=dispWin(1) & rasterXInd_ms{displayUnits(clusterNum)}<=dispWin(end))-timeAxis(1);
    winIndy=rasterYInd_ms{displayUnits(clusterNum)}...
        (rasterXInd_ms{displayUnits(clusterNum)}>=dispWin(1) & rasterXInd_ms{displayUnits(clusterNum)}<=dispWin(end));
    plot([winIndx;winIndx],[winIndy+clusterNum-1;winIndy+clusterNum],...
        'color',cmap(clusterNum,:),'LineStyle','-','LineWidth',1.8); % plot rasters% plot rasters
    plot(BP_periodBehavData_ms(timeAxis)/max(abs(BP_periodBehavData_ms(timeAxis)))+...
        clusterNum+0.5,'color','k','LineWidth',0.8);
    plot(whiskingPhaseExcerpt_ms(timeAxis)/max(abs(whiskingPhaseExcerpt_ms(timeAxis)))+...
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

%% Plot cross-correlation
figure('position',[602   537   560   420]); hold on
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms(displayUnits(clusterNum),:),BP_periodBehavData_ms,500,'coeff');
    %     figure('position',[602   537   560   420]);
    plot(lag,acor,'LineWidth',2); %xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
end
legend();
xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(displayUnits(clusterNum)))];...
    'Spike density function vs. Whisking angle'})

for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms(displayUnits(clusterNum),:),LP_periodBehavData_ms,1000,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Set point'})
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms(displayUnits(clusterNum),:),whiskingPhaseExcerpt_ms,100,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Whisking phase'})
    % hold on; [acor,lag] = xcorr(sdf,BP_periodBehavData,100,'unbiased'); plot(lag,acor);
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms(displayUnits(clusterNum),:),periodBehavData_ms,2500,'coeff');
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
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpt(BP_periodBehavData_ms',...
        spikeTimes(unitID==displayUnits(clusterNum))/double(samplingRate),params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(BP_periodBehavData_ms',...
        SDFs_ms(displayUnits(clusterNum),:)',params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskingPhaseExcerpt_ms',...
        SDFs_ms(displayUnits(clusterNum),:)',params);
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
traceNum=2; %3;
vIRtTrace=recordingTraces(traceNum,:);
midRec=round(length(vIRtTrace)/2)/30000;
vibrissaAngle=resample(HP_periodBehavData_ms,30,1); %foo=foo(1:length(recordingTrace{traceNum}));

figure('Color','white');
subplot(2,1,1); hold on
timeAxis=(1:length(vIRtTrace))/30000;
plot(timeAxis,vIRtTrace,'color','k','LineWidth',0.8);
%     correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
%     for clusterNum=1:numel(correspondingUnits)
%         plot(spikeTimes(unitID==correspondingUnits(clusterNum)),...
% correspondingUnits=find(ismember(titularChannels,keepTraces(traceNum)));
correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
for clusterNum=1:length(correspondingUnits)
    plot(spikeTimes(correspondingUnits(clusterNum),:)/samplingRate,...
        ones(size(spikeTimes(correspondingUnits(clusterNum),:),1),1)*-1500+clusterNum*100,'d');
end

set(gca,'ylim',[-1500 1000],'xlim',[661 665],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out');
xlabel('Time (s)'); ylabel('Firing rate (Hz)');

subplot(2,1,2); hold on
timeAxis=(1:length(vibrissaAngle))/30000;
plot(timeAxis,vibrissaAngle,'color','k','LineWidth',0.8);
for clusterNum=1:length(correspondingUnits)
    plot(spikeTimes(correspondingUnits(clusterNum),:)/samplingRate,...
        ones(size(spikeTimes(correspondingUnits(clusterNum),:),1),1)+clusterNum,'d');
end
set(gca,'xlim',[661 665],'ylim',[-15 15],'Box','off','Color','white','FontSize',10,...
    'FontName','Helvetica','TickDir','out');% axis 'tight'
xlabel('Time (s)');  ylabel('Angle (\circ)');

% cross correlation
bestUnit=2; %4;
[acor,lag] = xcorr(SDFs_ms{bestUnit},HP_periodBehavData_ms,500,'coeff');
figure('position',[602   537   560   420]);
plot(lag,acor,'color','k','LineWidth',2); xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(bestUnit))];'Spike density function vs. Whisking angle'})
