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
    whisker.Angle=behav.angle;
    whisker.Velocity=behav.velocity; 
    whisker.Phase=-behav.phase; 
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
    whisker.Angle=behav.whiskerTrackingData.Angle; %whisker.Angle=fillmissing(behav.whiskerTrackingData.Angle_BP,'nearest');
    whisker.Velocity=behav.whiskerTrackingData.Velocity;
    whisker.Phase=behav.whiskerTrackingData.Phase;
    whisker.Amplitude=behav.whiskerTrackingData.Amplitude;
    whisker.Frequency=behav.whiskerTrackingData.Freq;
    whisker.SetPoint=behav.whiskerTrackingData.SetPoint;
    whisker.Timestamps=behav.whiskerTrackingData.Timestamp;
    %% compute whisking frequency (different from instantaneous frequency
    whisksIdx = bwconncomp(diff(whisker.Phase)>0);
    peakIdx = zeros(1,length(whisker.Velocity));
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
    whisker.Frequency=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);

    %% other behavior data
    if ~isempty(behav.breathing)
        breathing.data=double(behav.breathing');
        breathing.data=breathing.data*(range(whisker.SetPoint)/range(breathing.data));
        breathing.ts=linspace(0,ephys.recInfo.duration_sec,numel(behav.breathing)); 
    end
    
    %If conforming to Kyle's convention
%     if max(whisker.Phase)<=pi
%         whisker.Phase=-whisker.Phase;
%     end
    %% compute rasters
    % aim for same length for ephys traces and behavior data
    [ephys.rasters,unitList]=EphysFun.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
        1,size(whisker.Angle,2)); %int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));
    %% compute spike density functions
    ephys.spikeRate=EphysFun.MakeSDF(ephys.rasters);
    % ephys.spikeRate=decimate(ephys.spikeRate,2);
    
    %create timeline
    ephys.timestamps = 0:0.001:ephys.recInfo.duration_sec;   
    
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

%% make sure behavior and spike traces have same length
cropTraces = false;
if cropTraces
    if size(ephys.spikeRate,2)~=numel(whisker.Angle)
        % check what's up
        if size(ephys.spikeRate,2)<size(whisker.Angle,2)
            whisker.Angle=whisker.Angle(:,1:size(ephys.spikeRate,2));
            whisker.Velocity=whisker.Velocity(:,1:size(ephys.spikeRate,2));
            whisker.Phase=whisker.Phase(:,1:size(ephys.spikeRate,2));
            whisker.Amplitude=whisker.Amplitude(:,1:size(ephys.spikeRate,2));
            whisker.Frequency=whiskerFrequency(:,1:size(ephys.spikeRate,2));
            whisker.SetPoint=whisker.SetPoint(:,1:size(ephys.spikeRate,2));
        else
            ephys.spikeRate=ephys.spikeRate(:,1:length(whisker.Angle));
            ephys.rasters=ephys.rasters(:,1:length(whisker.Angle));
        end
    end
end

%% whisking epochs (based on first trace, if multiple whisker tracked)
ampThd=18; %12; %18 %amplitude threshold 
freqThld=1; %frequency threshold
minBoutDur=1000; %500; % 1000 % minimum whisking bout duration: 1s
whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
   whisker.Amplitude(1,:),whisker.Frequency(1,:),...
   ampThd, freqThld, minBoutDur); 
whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
whiskingEpochsList=bwconncomp(whiskingEpochs);
[~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList.PixelIdxList),'descend');
whiskingEpochsList.PixelIdxListSorted=whiskingEpochsList.PixelIdxList(wBoutDurSort);

if false
    figure; hold on;
    plot(whisker.Angle);
    plot(whiskingEpochs*nanstd(whisker.Angle)+nanmean(whisker.Angle))
    % plot(whisker.Phase*nanstd(whisker.Angle)/2+nanmean(whisker.Angle));
end

%% decide which units to keep: 
%most frequent units
% mostFrqUnits=EphysFun.FindBestUnits(ephys.spikes.unitID,1);%keep ones over x% spikes
%most frequent units during whisking periods
reconstrUnits=ephys.rasters(:,whiskingEpochs).*(1:size(ephys.rasters,1))';
reconstrUnits=reshape(reconstrUnits,[1,size(reconstrUnits,1)*size(reconstrUnits,2)]);
reconstrUnits=reconstrUnits(reconstrUnits>0);
mostFrqUnits=EphysFun.FindBestUnits(reconstrUnits,1);
keepUnits=ismember(unitList,mostFrqUnits);
keepTraces=unique(ephys.spikes.preferredElectrode(ismember(ephys.spikes.unitID,unitList(keepUnits))));
ephys.selectedUnits=find(keepUnits);
%all of them
ephys.selectedUnits=unitList; %all units 
% only SU
[unitQuality,RPVIndex]=SSQualityMetrics(ephys.spikes);
unitQuality=[unique(double(ephys.spikes.unitID)),unitQuality];
unitIdx=ismember(ephys.spikes.unitID,unitQuality(unitQuality(:,2)>0.6,1));
unitQuality(unitQuality(:,2)>0.6,3)=hist(double(ephys.spikes.unitID(unitIdx)),...
    unique(double(ephys.spikes.unitID(unitIdx))))/sum(unitIdx);
qualityUnits=unitQuality(unitQuality(:,2)>0.6 & unitQuality(:,3)>0.01,:);
ephys.selectedUnits=qualityUnits(:,1);
% add manually, e.g., 
% ephys.selectedUnits=[ephys.selectedUnits;54]; 1;2;19];
% set just one
% ephys.selectedUnits=12;

%% organize selectedUnits by depth
if size(ephys.recInfo.probeGeometry,2)>size(ephys.recInfo.probeGeometry,1)
    ephys.recInfo.probeGeometry=ephys.recInfo.probeGeometry';
end
unitLoc=nan(numel(ephys.selectedUnits),2);
for unitNum=1:numel(ephys.selectedUnits)
    unitID=ephys.selectedUnits(unitNum);
%     [foo,faa]=hist(double(ephys.spikes.preferredElectrode(ephys.spikes.unitID==unitID)),...
%         double(unique(ephys.spikes.preferredElectrode(ephys.spikes.unitID==unitID))))
    prefElec=mode(ephys.spikes.preferredElectrode(ephys.spikes.unitID==unitID));
    unitLoc(unitNum,:)=ephys.recInfo.probeGeometry(prefElec,:);
end
[~,unitDepthOrder]=sort(unitLoc(:,2));
ephys.selectedUnits=ephys.selectedUnits(unitDepthOrder);
ephys.unitCoordinates=unitLoc(unitDepthOrder,:);
% rasters=ephys.rasters(ephys.selectedUnits,:);
% save([ephys.recInfo.sessionName 'Traces_Rasters.mat'],'whisker.Angle','whisker.Phase','rasters')

% only + PT (see below for PT)

%% may try to clean up RPVs
% FixRPV(ephys.spikes,RPVIndex);

%% whisking mode 
if false
    % Four main modes:
    % foveal: high frequency > 10Hz, medium amplitude >25 <35, high setpoint/angular values >70 at start of whisk
    % exploratory: lower frequency < 10, high amplitude >35, medium setpoint/angular values ?
    % resting: lower frequency < 10Hz, low/medium amplitude <25, low setpoint/angular values <60
    % twiches: high frequency > 25Hz, low amplitude <10, low setpoint/angular values <70 at start of whisk
    whiskingEpochs=WhiskingFun.FindWhiskingModes(whisker.Angle,whisker.Velocity,whisker.Amplitude,whiskerFrequency,whisker.SetPoint);
    figure; hold on
    plot(whisker.Angle)
    for wModeNum=1:numel(whiskingEpochs)
        whiskModeVals=nan(size(whisker.Angle));
        whiskModeVals(whiskingEpochs(wModeNum).index)=whisker.Angle(whiskingEpochs(wModeNum).index);
        plot(whiskModeVals);
    end
    legend('angle',whiskingEpochs.type)
end

%% Spike sorting integrity over whisking epochs
if false
% plot waveforms, ACG ISI
    NBC_Plots_Continuity(ephys.spikes,ephys.selectedUnits,whisker.Angle,1:numel(whiskingEpochs),false); %whiskingEpochs
end

%% Consolidate epochs for summary figure
%% Make collision test for antidromic stim

%% Overview plot
NBC_Plots_Overview(whisker,whiskingEpochs,breathing,ephys,pulses.TTLTimes,false,false);

%% Check Phototagging summary 
% ephys.selectedUnits=[60 23]; 10; 2; 37; %12;
if ~isfield(pulses,'duration'); pulses.duration=0.01; end
PhotoTagPlots(ephys,pulses); %Implement SALT test
% PTunits=[12,26,37];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any tuning computation has to be done outside of stimulation periods 
pulseMask=false(1,size(whisker.Angle,2));
pulseMask(round(pulses.TTLTimes(1)*1000):round(pulses.TTLTimes(end)*1000))=true;
whiskingEpochs_m=whiskingEpochs;whiskingEpochs_m(pulseMask)=false;

% ephys.selectedUnits=2; 2; 37; %12;
% ephys.selectedUnits=find(keepUnits);
% ephys.selectedUnits=qualityUnits(:,1);

%% Phase tuning - Individual plots
phaseTuning=NBC_Plots_PhaseTuning(whisker.Angle,whisker.Phase,ephys,whiskingEpochs_m,'whisking',false,false); %whiskingEpochs_m %ephys.spikeRate

% Set point phase tuning 
setpointPhase=WhiskingFun.ComputePhase(whisker.SetPoint,1000);
NBC_Plots_PhaseTuning(whisker.SetPoint,setpointPhase,ephys,whiskingEpochs_m,'setpoint oscillation',false,false);

% Breathing phase tuning 
breathingPhase=WhiskingFun.ComputePhase(smooth(breathing.data,50)',1000);
NBC_Plots_PhaseTuning(breathing.data,breathingPhase,ephys,whiskingEpochs_m,'breathing oscillation',false,false);

%% Angle tuning - Individual plots
whiskingEpochs_m=true(size(whisker.Angle));whiskingEpochs_m(pulseMask)=false;
NBC_Plots_AngleTuning(whisker.Angle,ephys,whiskingEpochs_m,'whisker',false,false);

%% Power spectrum and coherence
unitNum=46;
NBC_Plots_Coherence(whisker.Angle,whisker.Phase,ephys,whiskingEpochs,unitNum,cmap);

%% make video of whisking bouts
boutNum=4 ; %19 7
cellNum=54; %37 26 
traceIndex=whiskingEpochsList.PixelIdxList{boutNum}; %350000:352000;%whiskingEpochsList.PixelIdxList{2}
[boutFrames,frameIndex]=WhiskingBoutVideo(ephys.recInfo.likelyVideoFile,ephys.recInfo.dirName,...
    traceIndex,behav.vidTimes,false);
%remove any drift by adding frames when fps > intended fps  (crude)
boutTimeLine=behav.vidTimes(frameIndex);boutTimeLine=boutTimeLine-boutTimeLine(1);
extraFrameIdx=find(diff(round(boutTimeLine))>mode(diff(round(boutTimeLine))));  
for frameIdx=1:numel(extraFrameIdx)
boutFrames=[boutFrames(1:extraFrameIdx(frameIdx)+frameIdx-1),...
    boutFrames(extraFrameIdx(frameIdx)+frameIdx-1),...
    boutFrames(extraFrameIdx(frameIdx)+frameIdx:end)];
end
boutFrames=boutFrames(1:end-numel(extraFrameIdx));
% vidDims=size(boutFrames(1).cdata);
% figure('position',[1500 450  vidDims(2) vidDims(1)],'color','k');
% movie(boutFrames,1,500);

% with trace added
% add audiovidDims=size(boutFrames(1).cdata);
spikeTimes = movmean(ephys.rasters(cellNum,traceIndex),2);
spikeTimes = logical(spikeTimes(1:2:end));
figure('position',[500 450  vidDims(2) vidDims(1)],'color','k');
boutFrames=FrameByFrame_Overlay(boutFrames,[whisker.Phase(traceIndex(1:2:end));spikeTimes]); %whisker.Angle

%add audio for a given cell, given speed
slowFactor=20; %500/25

% Code perso
% wBoutAudio=WhiskingBoutAudio(ephys.rasters(cellNum,traceIndex(1:5000)),1000,20*slowFactor);
% FacePro
wBoutAudio = FacePro.MakeSpikeAudio(ephys.rasters(cellNum,traceIndex(1:5000)),...
    slowFactor*10, [1 2500], 1000);
% 
% samplingRatio=round(numel(traceIndex)/numel(boutFrames));
% wBoutAudio = zeros(200,numel(boutFrames));
% spikeTimes = round(find(ephys.rasters(cellNum,traceIndex))/samplingRatio);
% waveforms = gausswin(20); 
% 
% for spikeNum = 1 : numel(spikeTimes)
% %     wBoutAudio(spikeTimes(spikeNum)-4:spikeTimes(spikeNum)+size(waveforms,1)-5) = waveforms;
%     wBoutAudio(1:20,spikeTimes(spikeNum)) = waveforms;
% end
% wBoutAudio = int16(wBoutAudio / max(abs(wBoutAudio(:))) * double(intmax('int16')));

% figure; plot(wBoutAudio)
% figure; imagesc(wBoutAudio)

%write video
frameRate=500/slowFactor;
videoFWriter = vision.VideoFileWriter(fullfile(cd,...
    [ephys.recInfo.likelyVideoFile(1:end-4) '_Bout' num2str(boutNum)...
    'x' num2str(slowFactor) '_Unit' num2str(cellNum) '_FP_PhaseTuning.avi']));
videoFWriter.FrameRate =frameRate ;
videoFWriter.AudioInputPort = true;
videoFWriter.VideoCompressor = 'None (uncompressed)'; % 'MJPEG Compressor';

for frameNum = 1 : 2500 %numel(boutFrames)
%     videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(:,(frameNum*2)-1));  %if double the number of frames
    videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(:,frameNum)); %if matrix form factor
%     videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(frameNum));
end
release(videoFWriter);

%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Older code - to be removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Is it tuned to Phase? Polar plot summary
phaseTuning=NBC_Plots_PhaseTuning_PolarPlots(whisker.Phase',...%whiskingPhase',...
    whiskingEpochs_m,ephys.spikeRate(ephys.selectedUnits,:),false,ephys.recInfo.sessionName); %ephys.spikeRate

%% Angle tuning
NBC_Plots_SpikingWhiskAngleTuning(whisker.Angle',whisker.Phase',...
    whiskingEpochs_m,ephys.rasters,false,ephys.recInfo.sessionName) %periodBeh

%% whisking vs spikes plot
% whiskerTraces_whiskFreq=whisker.Angle; %WhiskingFun.BandPassBehavData(whisker.Angle,1000,[4 10]);
NBC_Plots_SpikesWhisking(whisker.Angle,whisker.Velocity,whisker.Phase,whiskingEpochs_m,ephys.rasters(mostFrqUnits,:));


%% whisking snapshot figures
NBC_Plots_WhiskingPhaseVideoFrame(whisker.Angle,vidTimes,...
    whisker.Phase,ephys.rasters,ephys.spikeRate);

%% overlay whisking and spikes on video
NBC_Plots_SpikesWhiskingOverlayOnVideo(whisker.Angle,vidTimes,...
    ephys.rasters,ephys.spikeRate); %cursor_1,cursor_2

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
            whiskerAngle_m(whiskNum,whiskingEpochs_m)',ephys.rasters(unitNum,whiskingEpochs_m)',params);
        
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
%         NBC_Plots_PhaseTuning(whiskingPhase(bestWhisker,:),whiskingEpochs_m,...
%             ephys.rasters(unitNum,:),false,['Unit ' str2double(unitNum) ephys.recName(1:end-19)])
%         NBC_Plots_SpikingWhiskAngleCC(periodBehavData(bestWhisker,:),whiskingPhase(bestWhisker,:),...
%             whiskingEpochs_m,ephys.rasters(unitNum,:),SDFs(unitNum,:),false,ephys.recName)
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
    whiskingEpochsInfo= bwconncomp(whiskingEpochs) ;
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
        whisker.Angle(whiskingEpochs)',500,'coeff');
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
        whisker.Phase(whiskingEpochs)',100,'coeff');
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
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpt(whisker.Angle,...
        spikeTimes{displayUnits(clusterNum)}/double(ephys.spikes.samplingRate),params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whisker.Angle,...
        ephys.spikeRate(displayUnits(clusterNum),:)',params);
    figure; subplot(311); plot(sgFreq,Coherence);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    Coherence(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [Coherence,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whisker.Phase,...
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

