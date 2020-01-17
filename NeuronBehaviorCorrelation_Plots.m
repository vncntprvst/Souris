%% NeuronBehaviorCorrelation_Plots
% At this point, loaded data should be processed, whisking traces filtered,
% phase computed, spike sorted, etc

clearvars;
%% First, get data
dirFiles=dir;
processedDataFiles=cellfun(@(x) contains(x,'processedData') ||...
    contains(x,'UnitExplorer'), {dirFiles.name});
if sum(processedDataFiles) 
    load(fullfile(dirFiles(processedDataFiles).folder,dirFiles(processedDataFiles).name));
end
if exist('ue','var')
    % extract data from ue object 
    [ephys,behav]=GetData_From_UnitExplorer(ue);
    ephys.rasters=ephys.spikeRasters; ephys=rmfield(ephys,'spikeRasters');
    ephys.spikeRate=EphysFun.MakeSDF(ephys.rasters);
    ephys.rasters=logical(sum([ephys.rasters(1:2:end);ephys.rasters(2:2:end)]));
    ephys.spikeRate=decimate(ephys.spikeRate,2);
    ephys.selectedUnits=1;
    ephys.recInfo.sessionName=ue.sessionName;
    whiskerAngle=behav.angle;
    whiskerVelocity=behav.velocity; 
    whiskerPhase=-behav.phase; 
    % Convention here: 
    %   Max protraction: 0 Max retraction pi/-pi
    %   (countercycle and -pi as compared to Kleinfeld et al.'s convention, where 
    %   Max protraction: pi Max retraction 0/2pi
else
    if ~exist('ephys','var')
        [ephys,behav,pulses]=NeuronBehaviorCorrelation_LoadData;
        cd(fullfile('../../Analysis',ephys.recInfo.sessionName))
%         save([ephys.recInfo.sessionName '_processedData'],'ephys','behav','pulses','-v7.3');
    end
    %% whisking data
    whiskerAngle=behav.whiskerTrackingData.Angle_f;
    whiskerVelocity=behav.whiskerTrackingData.Velocity;
    whiskerPhase=behav.whiskerTrackingData.Phase;
    whiskerAmplitude=behav.whiskerTrackingData.Amplitude;
    whiskerFrequency=behav.whiskerTrackingData.Freq;
    
    %If conforming to Kyle's convention
%     if max(whiskerPhase)<=pi
%         whiskerPhase=-whiskerPhase;
%     end
    %% compute rasters
    % aim for same length for ephys traces and behavior data
    [ephys.rasters,unitList]=EphysFun.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
        ephys.spikes.samplingRate,size(whiskerAngle,1)); %int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));
    %% compute spike density functions
    ephys.spikeRate=EphysFun.MakeSDF(ephys.rasters);
    % ephys.spikeRate=decimate(ephys.spikeRate,2);
    %% Convert to UnitExplorer
    %     spikesTrialArrayObj = s.get_sorted_spike_times();
    %     % Combine WhiskerTrialLites with SweepArray
    %     whiskerTrialLiteArrayObj = Whisker.WhiskerTrialLiteArray(fullfile(pathMainFolder, 'workspace'));
    %     % Create UnitExplorer
    %     ue = PrV.UnitExplorer(spikesTrialArrayObj, whiskerTrialLiteArrayObj);
    %     % Save UnitExplorer
%     ue.SaveObject(pathMainFolder);
end

% make sure behavior and spike traces have same length
if size(ephys.spikeRate,2)~=numel(whiskerAngle)
    % check what's up
    if size(ephys.spikeRate,2)<size(whiskerAngle,1)
        whiskerAngle=whiskerAngle(1:size(ephys.spikeRate,2),:);
        whiskerVelocity=whiskerVelocity(1:size(ephys.spikeRate,2),:);
        whiskerPhase=whiskerPhase(1:size(ephys.spikeRate,2),:);
    else
        ephys.spikeRate=ephys.spikeRate(:,1:size(whiskerAngle,1));
        ephys.rasters=ephys.rasters(:,1:size(whiskerAngle,1));
    end
end

%% decide which units to keep: 

mostFrqUnits=EphysFun.FindBestUnits(ephys.spikes.unitID,3);%keep ones over x% spikes
keepUnits=ismember(unitList,mostFrqUnits);
keepTraces=unique(ephys.spikes.preferredElectrode(ismember(ephys.spikes.unitID,unitList(keepUnits))));
ephys.selectedUnits=find(keepUnits);
ephys.selectedUnits=unitList; %all units 

% ephys.spikeRate=ephys.spikeRate(keepUnits,:);
% ephys.rasters=ephys.rasters(keepUnits,:);
% keepUnits{1}=ismember(unitList,mostFrqUnits); % most numerous units
% keepUnits{2}=~isnan(phaseTuning(:,1)); % phase tuned units
% % restrict to selected units
% ephys.rasters{1}=ephys.rasters(keepUnits{1},:);
% ephys.rasters{2}=ephys.rasters(keepUnits{2},:);

%% whisking epochs (based on first trace, if multiple whisker tracked)
ampThd=18; %amplitude threshold
freqThld=1; %frequency threshold
minBoutDur=1000; % minimum whisking bout duration: 1s
whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
   whiskerAmplitude(1,:),whiskerFrequency(1,:),...
   ampThd, freqThld, minBoutDur); 
whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
figure; hold on;
plot(whiskerAngle); plot(whiskingEpochs*std(whiskerAngle)+mean(whiskerAngle))
plot(whiskerPhase*std(whiskerAngle)/2+mean(whiskerAngle));

%% Check Phototagging summary 
ephys.selectedUnits=9; %keepUnits;
PhotoTagPlots(ephys,pulses);

%% Is it tuned to Phase? Polar plot summary
phaseTuning=NBC_Plots_PhaseTuning_PolarPlots(whiskerPhase',...%whiskingPhase',...
    whiskingEpochs,ephys.spikeRate(ephys.selectedUnits,:),false,ephys.recInfo.sessionName); %ephys.spikeRate

%% Phase tuning - Individual plots
phaseTuning=NBC_Plots_PhaseTuning(whiskerAngle',whiskerPhase',ephys,whiskingEpochs,false); %ephys.spikeRate

%% whisking vs spikes plot
% whiskerTraces_whiskFreq=whiskerAngle; %WhiskingFun.BandPassBehavData(whiskerAngle,1000,[4 10]);
NBC_Plots_SpikesWhisking(whiskerAngle,whiskerVelocity,whiskerPhase,whiskingEpochs,ephys.rasters(mostFrqUnits,:));

%% Angle tuning
NBC_Plots_SpikingWhiskAngleTuning(whiskerAngle',whiskerPhase',...
    whiskingEpochs,ephys.rasters,false,ephys.recInfo.sessionName) %periodBehavData 

%% whisking snapshot figures
NBC_Plots_WhiskingPhaseVideoFrame(whiskerAngle,vidTimes,...
    whiskerPhase,ephys.rasters,ephys.spikeRate);

%% overlay whisking and spikes on video
NBC_Plots_SpikesWhiskingOverlayOnVideo(whiskerAngle,vidTimes,...
    ephys.rasters,ephys.spikeRate); %cursor_1,cursor_2






% figure; hold on
% % plot(periodBehavData(1,:))
% % plot(LP_periodBehavData(1,:))
% plot(whiskerTraces_whiskFreq(1,:))
% % plot(HP_periodBehavData(1,:))
% plot(whiskingPhase(1,:))
% plot(whiskerTraces_breathFreq(1,:))

% parameters for coherence computation
params.Fs=1000; % sampling frequency
params.fpass=[4 20]; %peakWhiskFreq; % band of frequencies to be kept
params.tapers=[2 3]; % taper parameters (corresponds to  % time-bandwidth product of 2: NW=2 [NW 2*NW-1]
params.pad=1; % pad factor for fft
params.err=[2 0.05];
params.trialave=0;
global CHRONUXGPU
CHRONUXGPU = 1;
        
for unitNum=1:size(ephys.rasters,1)
    for whiskNum=1 %:size(whiskerTraces_whiskFreq,1)
        %% Find whisking periods of at least 500ms
        % peakWhisking=WhiskingFun.FindWhiskingBouts(periodBehavData(1,:));
        % peakWhiskingIdx=find(peakWhisking(1,:)==max(peakWhisking(1,:)));
        [whiskingEpochsIdx,wAmplitude,setPoint]=WhiskingFun.FindWhiskingEpochs(...
            whiskerTraces_whiskFreq(:,whiskNum)',whiskerPhase(:,whiskNum)',500);
        % plot(whiskingEpochsIdx*50)
        if sum(whiskingEpochsIdx)==0
            continue
        else
            whiskingEpochs=whiskingEpochsIdx;
        %restrict to first 90s
        % whiskingEpochs(90000:end)=0;
        end
        %% find peak whisking frequency
        whiskFreqSpectrum=abs(fft(whiskerTraces_whiskFreq(whiskNum,whiskingEpochs)));
        numVals=numel(whiskerTraces_whiskFreq(whiskNum,whiskingEpochs));
        whiskFreqSpectrum = smoothdata(whiskFreqSpectrum(1:floor(numVals/2)+1),...
            'movmean',round(numel(whiskFreqSpectrum)/2/200));
        % peakWhiskFreq(2:end-1) = 2*peakWhiskFreq(2:end-1);
        freqArray = 1000*(0:(numVals/2))/numVals;
%       figure;  plot(freqArray,whiskFreqSpectrum)
        % title('Single-Sided Amplitude Spectrum of X(t)')
        peakWhiskFreq=freqArray(whiskFreqSpectrum==max(whiskFreqSpectrum));

        %% agregate whisker motion components
        whiskerMotion=[whiskerPhase(whiskNum,:)',wAmplitude',setPoint'];
        %% Coherence
        % we want the coherence at the peak frequency

        try
        [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(...
            whiskerTraces_whiskFreq(whiskNum,whiskingEpochs)',ephys.rasters(unitNum,whiskingEpochs)',params);
        
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
%         NBC_Plots_PhaseTuning(whiskingPhase(bestWhisker,:),whiskingEpochs,...
%             ephys.rasters(unitNum,:),false,['Unit ' str2double(unitNum) ephys.recName(1:end-19)])
%         NBC_Plots_SpikingWhiskAngleCC(periodBehavData(bestWhisker,:),whiskingPhase(bestWhisker,:),...
%             whiskingEpochs,ephys.rasters(unitNum,:),SDFs(unitNum,:),false,ephys.recName)
%     end
end

% figure;
% subplot(3,1,1); plot(sgFreq,Coherence);
% subplot(3,1,2); plot(sgFreq,10*log10(S1));
% subplot(3,1,3); plot(sgFreq,10*log10(S2))
% Coherence(S1==max(S1))




 %% get average waveform
    avWaveform=nan(numel(bestUnits), size(ephys.spikes.waveforms,2));
    %     figure; hold on
    for unitNum=1:numel(bestUnits)
        avWaveform(unitNum,:)=nanmean(ephys.spikes.waveforms(...
            ephys.spikes.unitID==bestUnits(unitNum),:));
        %         plot(avWaveform(unitNum,:))
    end
    % neuron with biggest waveform
    biggestUnit=sum(abs(diff(avWaveform,[],2)),2)==max(sum(abs(diff(avWaveform,[],2)),2));
    
    % filters trace
    whiskingEpochsInfo= bwconncomp(whiskingEpochsIdx) ;
    whiskBoutDuration=cellfun(@(whiskBout) numel(whiskBout), whiskingEpochsInfo.PixelIdxList);
    longestWhiskBout=whiskBoutDuration==max(whiskBoutDuration);
    longWhiskBoutInit=whiskingEpochsInfo.PixelIdxList{longestWhiskBout}(1);
    displayWindow=max([longWhiskBoutInit-1000 1]):min([longWhiskBoutInit+58999 numel(whiskingEpochsIdx)-1]) ;
    filtTraces=PreProcData(ephys.traces(:,displayWindow*30),ephys.spikes.samplingRate,{'CAR','all'});
    filtTraces=filtTraces(keepTraces(biggestUnit),:);
    




%% Is WR bursting reliable? - Plot spike times with whisking angle
% NBC_Plots_SpikingWhiskAngleCC(periodBehavData,whiskingPhase,...
%     whiskingEpochs,ephys.rasters,SDFs,false,ephys.recName)




%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[602   537   560   420]); hold on
for clusterNum=1:length(ephys.selectedUnits)
    [acor,lag] = xcorr(ephys.spikeRate(ephys.selectedUnits,whiskingEpochs),...
        whiskerAngle(whiskingEpochs)',500,'coeff');
    %     figure('position',[602   537   560   420]);
    plot(lag,acor,'LineWidth',2); 
    %xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
end
legend();
xlabel('Lag (ms)');set(gca,'ylim',[-0.5 1])
title({['Cross correlation for vIRt unit ' num2str(ephys.selectedUnits)];...
    'Spike density function vs. Whisking angle'})

for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(ephys.spikeRate{displayUnits(clusterNum)},LP_periodBehavData,1000,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Set point'})
end
for clusterNum=1:length(ephys.selectedUnits)
    [acor,lag] = xcorr(ephys.spikeRate(ephys.selectedUnits,whiskingEpochs),...
        whiskerPhase(whiskingEpochs)',100,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(ephys.selectedUnits)];'Spike density function vs. Whisking phase'})
    % hold on; [acor,lag] = xcorr(sdf,BP_periodBehavData,100,'unbiased'); plot(lag,acor);
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(ephys.spikeRate{displayUnits(clusterNum)},periodBehavData,2500,'coeff');
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
displayUnits=find(keepUnits);
spikeTimes=cell(numel(keepUnits),1);
for clusterNum=1:length(keepUnits)
    spikeTimes{clusterNum}=double(ephys.spikes.times(ephys.spikes.unitID==clusterNum));
end

params.Fs=1000; % sampling frequency
params.fpass=[0 25]; % band of frequencies to be kept
params.tapers=[2 3]; % taper parameters
params.pad=1; % pad factor for fft
params.err=[2 0.05];
params.trialave=0;

for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpt(whiskerAngle,...
        spikeTimes{displayUnits(clusterNum)}/double(ephys.spikes.samplingRate),params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskerAngle,...
        ephys.spikeRate(displayUnits(clusterNum),:)',params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskerPhase,...
        ephys.spikeRate(displayUnits(clusterNum),:)',params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end



%% mock up units to test code
% mock up protraction unit
% fakeUnitSpikes=zeros(1,numel(periodBehavData));
% fakeUnitSpikes([0 diff(BP_periodBehavData)]>0 &...
% [0,0,diff(diff(BP_periodBehavData))]<0 &...
% whiskingPhase<0)=1;
% % mock up retraction unit
% fakeUnitSpikes=zeros(1,numel(periodBehavData));
% fakeUnitSpikes([0 diff(BP_periodBehavData)]<0 &...
% [0,0,diff(diff(BP_periodBehavData))]>0 &...
% whiskingPhase>0)=1;
% unitSpikes=fakeUnitSpikes; unitSDF=GaussConv(unitSpikes,5)*1000;

