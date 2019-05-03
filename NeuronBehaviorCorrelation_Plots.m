% NeuronBehaviorCorrelation_Plots
clearvars;
%% First, get data
dirFiles=dir;
processedData=cellfun(@(x) contains(x,'processedData'), {dirFiles.name});
if sum(processedData)
    load(fullfile(dirFiles(processedData).folder,dirFiles(processedData).name));
else
    [ephys,behav]=NeuronBehaviorCorrelation_GatherData;
    save('processedData','ephys','behav','-v7.3');
end

% figure; hold on
% plot(behav.vidTimes,behav.whiskerTrackingData(1,:)-mean(behav.whiskerTrackingData(1,:)));
% bestUnitTimes=ephys.spikes.times(ephys.spikes.unitID==4);
% plot(bestUnitTimes,zeros(numel(bestUnitTimes),1),'dk')

bestUnits=EphysFunctions.FindBestUnits(ephys.spikes.unitID);
[spikeRasters_ms,unitList]=EphysFunctions.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
    ephys.samplingRate,int32(size(ephys.traces,2)/ephys.samplingRate*1000));
% decide which units to keep
keepUnits=ismember(unitList,bestUnits); %keepUnits=unitList
keepTraces=unique(ephys.spikes.preferredElectrode(ismember(ephys.spikes.unitID,unitList(keepUnits))));
% compute spike density functions
SDFs_ms=EphysFunctions.MakeSDF(spikeRasters_ms);
spikeRasters_ms=spikeRasters_ms(keepUnits,:);
SDFs_ms=SDFs_ms(keepUnits,:);
% for whisking data resampling, use ephys sampling rate, if videao times are
% based on TTLs
periodBehavData_ms=WhiskingAnalysisFunctions.ResampleBehavData...
    (behav.whiskerTrackingData,behav.vidTimes,ephys.samplingRate);
% make sure behavior and spike traces have same length
if size(SDFs_ms,2)~=size(periodBehavData_ms,2)
    % check what's up
    if size(SDFs_ms,2)<size(periodBehavData_ms,2)
        periodBehavData_ms=periodBehavData_ms(:,1:size(SDFs_ms,2));
    else
        SDFs_ms=SDFs_ms(:,1:size(periodBehavData_ms,2));
        spikeRasters_ms=spikeRasters_ms(:,1:size(periodBehavData_ms,2));
    end
end
%% filtered whisking traces
whiskerTraces_whiskFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(periodBehavData_ms,1000,[4 20]);
whiskerTraces_breathFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(periodBehavData_ms,1000,[1 4]);
% HP_periodBehavData_ms=WhiskingAnalysisFunctions.HighPassBehavData(periodBehavData_ms);
% LP_periodBehavData_ms=WhiskingAnalysisFunctions.LowPassBehavData(periodBehavData_ms);

%% find phase
whiskingPhase_ms=WhiskingAnalysisFunctions.ComputePhase(whiskerTraces_whiskFreq_ms); %,whiskingEpochs);

% figure; hold on
% % plot(periodBehavData_ms(1,:))
% % plot(LP_periodBehavData_ms(1,:))
% plot(whiskerTraces_whiskFreq_ms(1,:))
% % plot(HP_periodBehavData_ms(1,:))
% plot(whiskingPhase_ms(1,:))
% plot(whiskerTraces_breathFreq_ms(1,:))

% parameters for coherence computation
params.Fs=1000; % sampling frequency
params.fpass=[4 20]; %peakWhiskFreq; % band of frequencies to be kept
params.tapers=[2 3]; % taper parameters (corresponds to  % time-bandwidth product of 2: NW=2 [NW 2*NW-1]
params.pad=1; % pad factor for fft
params.err=[2 0.05];
params.trialave=0;
global CHRONUXGPU
CHRONUXGPU = 1;
        
for unitNum=1:size(spikeRasters_ms,1)
    for whiskNum=1 %:size(whiskerTraces_whiskFreq_ms,1)
        %% Find whisking periods of at least 500ms
        % peakWhisking_ms=WhiskingAnalysisFunctions.FindWhiskingBouts(periodBehavData_ms(1,:));
        % peakWhiskingIdx=find(peakWhisking_ms(1,:)==max(peakWhisking_ms(1,:)));
        [whiskingEpochsIdx,wAmplitude,setPoint]=WhiskingAnalysisFunctions.FindWhiskingEpochs(...
            whiskerTraces_whiskFreq_ms(whiskNum,:),whiskingPhase_ms(whiskNum,:),500);
        % plot(whiskingEpochsIdx*50)
        if sum(whiskingEpochsIdx)==0
            continue
        else
            whiskingEpochs=whiskingEpochsIdx;
        %restrict to first 90s
        % whiskingEpochs(90000:end)=0;
        end
        %% find peak whisking frequency
        whiskFreqSpectrum=abs(fft(whiskerTraces_whiskFreq_ms(whiskNum,whiskingEpochs)));
        numVals=numel(whiskerTraces_whiskFreq_ms(whiskNum,whiskingEpochs));
        whiskFreqSpectrum = smoothdata(whiskFreqSpectrum(1:floor(numVals/2)+1),...
            'movmean',round(numel(whiskFreqSpectrum)/2/200));
        % peakWhiskFreq(2:end-1) = 2*peakWhiskFreq(2:end-1);
        freqArray = 1000*(0:(numVals/2))/numVals;
%       figure;  plot(freqArray,whiskFreqSpectrum)
        % title('Single-Sided Amplitude Spectrum of X(t)')
        peakWhiskFreq=freqArray(whiskFreqSpectrum==max(whiskFreqSpectrum));

        %% agregate whisker motion components
        whiskerMotion=[whiskingPhase_ms(whiskNum,:)',wAmplitude',setPoint'];
        %% Coherence
        % we want the coherence at the peak frequency

        try
        [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(...
            whiskerTraces_whiskFreq_ms(whiskNum,whiskingEpochs)',spikeRasters_ms(unitNum,whiskingEpochs)',params);
        
%         smoothedCoherence=smoothdata(Coherence(:,1),'movmean',round(length(Coherence(:,1))/diff(params.fpass)));
figure; hold on 
plot(sgFreq,Coherence)
% plot(sgFreq,smoothedCoherence)
%        smoothedCoherence(find(sgFreq>=peakWhiskFreq,1),:)
       
%         allCoherence{unitNum}(:,whiskNum)=Coherence;
        catch
            continue
        end
    end
%     allCoherence{unitNum}
%     if logical(sum(allCoherence{unitNum}(1,:)>0.8))
%         bestWhisker=allCoherence{unitNum}(1,:)==max(allCoherence{unitNum}(1,:));
%         NBC_Plots_PhaseTuning(whiskingPhase_ms(bestWhisker,:),whiskingEpochs,...
%             spikeRasters_ms(unitNum,:),false,['Unit ' str2double(unitNum) ephys.recName(1:end-19)])
%         NBC_Plots_SpikingWhiskAngleCC(periodBehavData_ms(bestWhisker,:),whiskingPhase_ms(bestWhisker,:),...
%             whiskingEpochs,spikeRasters_ms(unitNum,:),SDFs_ms(unitNum,:),false,ephys.recName)
%     end
end

figure;
subplot(3,1,1); plot(sgFreq,Coherence);
subplot(3,1,2); plot(sgFreq,10*log10(S1));
subplot(3,1,3); plot(sgFreq,10*log10(S2))
Coherence(S1==max(S1))

%% Angle tuning
NBC_Plots_SpikingWhiskAngleTuning(whiskerTraces_whiskFreq_ms,whiskingPhase_ms,...
    whiskingEpochs,spikeRasters_ms,false,ephys.recName) %periodBehavData_ms 

%% Is it tuned to Phase? - Polar plot summary
NBC_Plots_PhaseTuning(whiskingPhase_ms,whiskingEpochs,spikeRasters_ms,false,ephys.recName)

%% whisking vs spikes plot
NBC_Plots_SpikesWhisking(whiskerTraces_whiskFreq_ms,whiskingPhase_ms,spikeRasters_ms);

%% whisking snapshot figures
NBC_Plots_WhiskingPhaseVideoFrame(periodBehavData_ms,vidTimes_ms,...
    whiskingPhase_ms,spikeRasters_ms,SDFs_ms);

%% overlay whisking and spikes on video
NBC_Plots_SpikesWhiskingOverlayOnVideo(periodBehavData_ms,vidTimes_ms,...
    spikeRasters_ms,SDFs_ms); %cursor_1,cursor_2



%% Is WR bursting reliable? - Plot spike times with whisking angle
% NBC_Plots_SpikingWhiskAngleCC(periodBehavData_ms,whiskingPhase_ms,...
%     whiskingEpochs,spikeRasters_ms,SDFs_ms,false,ephys.recName)




%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[602   537   560   420]); hold on
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs{displayUnits(clusterNum)},BP_periodBehavData,500,'coeff');
    %     figure('position',[602   537   560   420]);
    plot(lag,acor,'LineWidth',2); %xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
end
legend();
xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(displayUnits(clusterNum)))];...
    'Spike density function vs. Whisking angle'})

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
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpt(BP_periodBehavData',...
        spikeTimes{displayUnits(clusterNum)}/double(spikeData.samplingRate),params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(BP_periodBehavData',...
        SDFs{displayUnits(clusterNum)}',params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskingPhase',...
        SDFs{displayUnits(clusterNum)}',params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end



%% mock up units to test code
% mock up protraction unit
% fakeUnitSpikes=zeros(1,numel(periodBehavData_ms));
% fakeUnitSpikes([0 diff(BP_periodBehavData_ms)]>0 &...
% [0,0,diff(diff(BP_periodBehavData_ms))]<0 &...
% whiskingPhase_ms<0)=1;
% % mock up retraction unit
% fakeUnitSpikes=zeros(1,numel(periodBehavData_ms));
% fakeUnitSpikes([0 diff(BP_periodBehavData_ms)]<0 &...
% [0,0,diff(diff(BP_periodBehavData_ms))]>0 &...
% whiskingPhase_ms>0)=1;
% unitSpikes=fakeUnitSpikes; unitSDF=GaussConv(unitSpikes,5)*1000;

