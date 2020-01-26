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
    whiskerAngle=behav.whiskerTrackingData.Angle; %whiskerAngle=fillmissing(behav.whiskerTrackingData.Angle_BP,'nearest');
    whiskerVelocity=behav.whiskerTrackingData.Velocity;
    whiskerPhase=behav.whiskerTrackingData.Phase;
    whiskerAmplitude=behav.whiskerTrackingData.Amplitude;
%     whiskerFrequency=behav.whiskerTrackingData.Freq;
    whiskerSetPoint=behav.whiskerTrackingData.SetPoint;
    %% compute whisking frequency (different from instantaneous frequency
    whisksIdx = bwconncomp(diff(whiskerPhase)>0);
    peakIdx = zeros(1,length(whiskerVelocity));
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
    whiskerFrequency=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);

    %If conforming to Kyle's convention
%     if max(whiskerPhase)<=pi
%         whiskerPhase=-whiskerPhase;
%     end
    %% compute rasters
    % aim for same length for ephys traces and behavior data
    [ephys.rasters,unitList]=EphysFun.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
        ephys.spikes.samplingRate,size(whiskerAngle,2)); %int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));
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

%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% make sure behavior and spike traces have same length
if size(ephys.spikeRate,2)~=numel(whiskerAngle)
    % check what's up
    if size(ephys.spikeRate,2)<size(whiskerAngle,2)
        whiskerAngle=whiskerAngle(:,1:size(ephys.spikeRate,2));
        whiskerVelocity=whiskerVelocity(:,1:size(ephys.spikeRate,2));
        whiskerPhase=whiskerPhase(:,1:size(ephys.spikeRate,2));
        whiskerAmplitude=whiskerAmplitude(:,1:size(ephys.spikeRate,2));
        whiskerFrequency=whiskerFrequency(:,1:size(ephys.spikeRate,2));
        whiskerSetPoint=whiskerSetPoint(:,1:size(ephys.spikeRate,2));
    else
        ephys.spikeRate=ephys.spikeRate(:,1:length(whiskerAngle));
        ephys.rasters=ephys.rasters(:,1:length(whiskerAngle));
    end
end

%% decide which units to keep: 
%most frequent
mostFrqUnits=EphysFun.FindBestUnits(ephys.spikes.unitID,1);%keep ones over x% spikes
keepUnits=ismember(unitList,mostFrqUnits);
keepTraces=unique(ephys.spikes.preferredElectrode(ismember(ephys.spikes.unitID,unitList(keepUnits))));
ephys.selectedUnits=find(keepUnits);
%all of them
ephys.selectedUnits=unitList; %all units 
% only SU
[unitQuality,RPVIndex]=SSQualityMetrics(ephys.spikes);
unitQuality=[unique(ephys.spikes.unitID),unitQuality];
unitIdx=ismember(ephys.spikes.unitID,unitQuality(unitQuality(:,2)>0.6,1));
unitQuality(unitQuality(:,2)>0.6,3)=hist(ephys.spikes.unitID(unitIdx),unique(ephys.spikes.unitID(unitIdx)))/sum(unitIdx);
qualityUnits=unitQuality(unitQuality(:,2)>0.6 & unitQuality(:,3)>0.01,:);
ephys.selectedUnits=qualityUnits(:,1);
% plot ACG and ISI
figure;
for unitNum=1:numel(ephys.selectedUnits) %ephys.selectedUnits
    acgH=subplot(numel(ephys.selectedUnits),3,unitNum*3-2);
    EphysFun.PlotACG(ephys.spikes.unitID,ephys.spikes.times,ephys.selectedUnits(unitNum),...
        ephys.spikes.samplingRate,acgH,cmap)
    isiH=subplot(numel(ephys.selectedUnits),3,unitNum*3-1:unitNum*3);
    EphysFun.PLotISI(ephys.spikes.unitID,ephys.spikes.times,ephys.selectedUnits(unitNum),...
        ephys.spikes.samplingRate,isiH,cmap)
end
% only + PT (see below for PT)

%% may try to clean up RPVs
% FixRPV(ephys.spikes,RPVIndex);

%% whisking epochs (based on first trace, if multiple whisker tracked)
ampThd=18; %amplitude threshold
freqThld=1; %frequency threshold
minBoutDur=1000; % minimum whisking bout duration: 1s
whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
   whiskerAmplitude(1,:),whiskerFrequency(1,:),...
   ampThd, freqThld, minBoutDur); 
whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
figure; hold on;
plot(whiskerAngle); plot(whiskingEpochs*nanstd(whiskerAngle)+nanmean(whiskerAngle))
plot(whiskerPhase*nanstd(whiskerAngle)/2+nanmean(whiskerAngle));


%% whisking mode 
% Four main modes:
% foveal: high frequency > 10Hz, medium amplitude >25 <35, high setpoint/angular values >70 at start of whisk
% exploratory: lower frequency < 10, high amplitude >35, medium setpoint/angular values ?
% resting: lower frequency < 10Hz, low/medium amplitude <25, low setpoint/angular values <60
% twiches: high frequency > 25Hz, low amplitude <10, low setpoint/angular values <70 at start of whisk
whiskingEpochs=WhiskingFun.FindWhiskingModes(whiskerAngle,whiskerVelocity,whiskerAmplitude,whiskerFrequency,whiskerSetPoint);
figure; hold on
plot(whiskerAngle)
for wModeNum=1:numel(whiskingEpochs)
    whiskModeVals=nan(size(whiskerAngle));
    whiskModeVals(whiskingEpochs(wModeNum).index)=whiskerAngle(whiskingEpochs(wModeNum).index);
    plot(whiskModeVals);
end
legend('angle',whiskingEpochs.type)

%% Overview plot
NBC_Plots_Overview(whiskerAngle,whiskerPhase,whiskerSetPoint,ephys,false);

%% Check Phototagging summary 
PhotoTagPlots(ephys,pulses);
PTunits=[7];

%% Is it tuned to Phase? Polar plot summary
phaseTuning=NBC_Plots_PhaseTuning_PolarPlots(whiskerPhase',...%whiskingPhase',...
    whiskingEpochs,ephys.spikeRate(ephys.selectedUnits,:),false,ephys.recInfo.sessionName); %ephys.spikeRate

%% Phase tuning - Individual plots
phaseTuning=NBC_Plots_PhaseTuning(whiskerAngle,whiskerAngle,ephys,whiskingEpochs,false); %ephys.spikeRate

%% Phase tuning - Individual plots
phaseTuning=NBC_Plots_PhaseTuning(whiskerAngle,ephys,whiskingEpochs,false); %ephys.spikeRate

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

%% Power spectrum and coherence
unitNum=15;
fs = 1000; winL= 2000; overlapL= 1900; % nsc = floor(numel(whiskerAngle_m)/4.5); nov = floor(nsc/2);

whiskerAngle_m=fillmissing(whiskerAngle,'nearest');
wspectrumVals=spectrogram(whiskerAngle_m,winL,overlapL,[],fs,'yaxis');
wspectrumVals=abs(wspectrumVals);
wf_wps=sum(wspectrumVals(1:15,:)); %power spectrum at whisker frequency range
figure; axis2=subplot(2,1,1); hold on
plot(wf_wps)
%plot whisking epochs
misTvalSize=mod(numel(whiskingEpochs),numel(wf_wps));
epochTarray=linspace(1,numel(wf_wps),numel(whiskingEpochs)-misTvalSize);
epochVarray=whiskingEpochs(misTvalSize/2+1:end-misTvalSize/2)*nanstd(wf_wps)+nanmean(wf_wps);
plot(epochTarray,epochVarray,'r','linewidth',1.5)
xlabel('Whisking power spectrum')

for unitNum=1:numel(ephys.selectedUnits)
    unitRasters=ephys.rasters(ephys.selectedUnits(unitNum),:);
    spectrumVals_s=spectrogram(unitRasters,winL,overlapL,[],fs,'yaxis');
    spectrumVals_s=abs(spectrumVals_s);
    sum_spectrumVals_s(unitNum,:)=sum(spectrumVals_s(1:15,:));
end
wf_sps=mean(sum_spectrumVals_s([7,10,12,13,15],:)); %sum(sspectrumVals(1:15,:)); %power spectrum at whisker frequency range
axis1=subplot(2,1,2);  hold on
plot(wf_sps)
epochVarray=whiskingEpochs(misTvalSize/2+1:end-misTvalSize/2)*nanstd(wf_sps)+nanmedian(wf_sps);
plot(epochTarray,epochVarray,'r','linewidth',1.5)
xlabel('Spiking power spectrum')
linkaxes([axis1,axis2],'x')

figure; %imagesc(sspectrumVals)
spectrogram(unitRasters,winL,overlapL,[],fs,'yaxis');
set(gca,'ylim',[1 30]); %'YDir','normal'); 
caxis([-60 -25]);hold on
epochTarray=linspace(1,max(get(gca,'xlim')),numel(whiskingEpochs));
epochVarray=whiskingEpochs*30;
plot(epochTarray,epochVarray,'r','linewidth',1.5)

% coherence within each whisking epoch
whiskerPhase_m=fillmissing(whiskerPhase,'nearest');
whiskBoutList = bwconncomp(whiskingEpochs) ; 
durationThd=cellfun(@(x) length(x),whiskBoutList.PixelIdxList)>4000; %keep only whisking bouts longer than 4s 
whiskBoutList.PixelIdxList=whiskBoutList.PixelIdxList(durationThd);
whiskBoutList.NumObjects=sum(durationThd);
figure('color','w','name','Spike / Angle Coherence','position',[1000 70 400 1200]); hold on
histHand = gobjects(2,whiskBoutList.NumObjects);
cxy=nan(257,whiskBoutList.NumObjects);
binNum=10; shiftF=0.6;
for wEpoch=whiskBoutList.NumObjects:-1:1
    signalVals=whiskerPhase_m(whiskBoutList.PixelIdxList{wEpoch}); %whiskerAngle_m whiskerPhase_m
    signalVals=signalVals-mean(signalVals);
%     signalVals=detrend(unwrap(signalVals))-mean(detrend(unwrap(signalVals)));
    cxy(:,wEpoch)=mscohere(unitRasters(whiskBoutList.PixelIdxList{wEpoch}),...
        signalVals,500,250,[],fs);
    freqVals=linspace(0,fs/2,numel(cxy(:,wEpoch)));
    xVals = [freqVals freqVals(end) 0]; 
    yVals = [cxy(:,wEpoch)',0,0]; 
    fill(gca,xVals,yVals+(wEpoch*shiftF)-1,cmap(wEpoch,:),'FaceAlpha',0.4,'EdgeColor','w','LineWidth',1);
end
% delete(histHand)
set(gca,'xlim',[0 50],'ytick',(shiftF:shiftF:(shiftF*whiskBoutList.NumObjects))-0.5,...
    'ylim',[-0.5 whiskBoutList.NumObjects*shiftF],'yticklabel',1:whiskBoutList.NumObjects,...
    'tickdir','out','box','off','XGrid', 'on');
ylabel('Whisking epoch')
xlabel('Frequency (Hz)')
title('Spike / Phase Magnitude-squared Coherence')
% legend(num2str((1:whiskBoutList.NumObjects)'))

%ridge plot
lgd=mat2cell(num2str((1:whiskBoutList.NumObjects)'),ones(whiskBoutList.NumObjects,1),2);
ridgeplot(mat2cell(cxy',ones(whiskBoutList.NumObjects,1),257),binNum,[0 0.2],lgd,cmap)

% parameters for coherence computation
params.Fs=1000; % sampling frequency
params.fpass=[4 30]; % % band of frequencies to be kept
params.tapers=[2 3]; % taper parameters (corresponds to  % time-bandwidth product of 2: NW=2 [NW 2*NW-1]
params.pad=1; % pad factor for fft
params.err=[2 0.05];
params.trialave=0;
global CHRONUXGPU
CHRONUXGPU = 1;
        
for unitNum=1:size(ephys.rasters,1)
    for whiskNum=1 %:size(whiskerTraces_whiskFreq,1)
        %% Coherence
        try
        [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(...
            whiskerAngle_m(whiskNum,whiskingEpochs)',ephys.rasters(unitNum,whiskingEpochs)',params);
        
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
    avWaveform=nan(numel(ephys.selectedUnits), size(ephys.spikes.waveforms,2));
    figure; hold on
    for unitNum=1:numel(bestUnits)
            avWaveform=nan(7, size(ephys.spikes.waveforms,2));
        for profileNum=1:9
            avWaveform(profileNum,:)=nanmean(ephys.spikes.waveforms(...
                ephys.spikes.unitID==ephys.selectedUnits(unitNum),:,profileNum));
            plot(avWaveform(profileNum,:))
        end
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

