% NeuronBehaviorCorrelation_Plots

%% First, get data
dirFiles=dir;
processedData=cellfun(@(x) contains(x,'processedData'), {dirFiles.name});
if sum(processedData)
    load(fullfile(dirFiles(processedData).folder,dirFiles(processedData).name));
else
    [ephys,behav]=NeuronBehaviorCorrelation_GatherData;
    save('processedData','ephys','behav','-v7.3');
end

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
BP_periodBehavData_ms=WhiskingAnalysisFunctions.BandPassBehavData(periodBehavData_ms);
% HP_periodBehavData_ms=WhiskingAnalysisFunctions.HighPassBehavData(periodBehavData_ms);
% LP_periodBehavData_ms=WhiskingAnalysisFunctions.LowPassBehavData(periodBehavData_ms);

%% Find whisking periods
% peakWhisking_ms=WhiskingAnalysisFunctions.FindWhiskingBouts(periodBehavData_ms(1,:));
% peakWhiskingIdx=find(peakWhisking_ms(1,:)==max(peakWhisking_ms(1,:)));
whiskingPeriodIdx=WhiskingAnalysisFunctions.FindWhiskingPeriods(BP_periodBehavData_ms);
% whiskingPeriodIdxInfo=bwconncomp(whiskingPeriodIdx);
% whiskingPeriodDuration=cellfun(@(x) numel(x),whiskingPeriodIdxInfo.PixelIdxList);
% longestWhiskingPeriodIdx=whiskingPeriodIdxInfo.PixelIdxList{...
%     whiskingPeriodDuration==max(whiskingPeriodDuration)};
whiskingEpochs=whiskingPeriodIdx;

whiskingPhase_ms=WhiskingAnalysisFunctions.ComputePhase(BP_periodBehavData_ms,whiskingEpochs);
% figure; hold on 
% plot(periodBehavData_ms(1,:))
% % plot(LP_periodBehavData_ms(1,:))
% plot(BP_periodBehavData_ms(1,:))
% % plot(HP_periodBehavData_ms(1,:))
% plot(whiskingPhase_ms(1,:))
%% whisking vs spikes plot
NBC_Plots_SpikesWhisking(BP_periodBehavData_ms,whiskingPhase_ms,spikeRasters_ms);

%% whisking snapshot figures
NBC_Plots_WhiskingPhaseVideoFrame(periodBehavData_ms,vidTimes_ms,...
    whiskingPhase_ms,spikeRasters_ms,SDFs_ms);

%% overlay whisking and spikes on video 
NBC_Plots_SpikesWhiskingOverlayOnVideo(periodBehavData_ms,vidTimes_ms,...
    spikeRasters_ms,SDFs_ms); %cursor_1,cursor_2

%% Is it tuned to Phase? - Polar plot summary 
NBC_Plots_PhaseTuning(whiskingPhase_ms,whiskingEpochs,spikeRasters_ms,false,ephys.recName)

%% Is WR bursting reliable? - Plot spike times with whisking angle
NBC_Plots_SpikingWhiskAngleCC(periodBehavData_ms,whiskingPhase_ms,...
    whiskingPeriodIdx,spikeRasters_ms,false)


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

