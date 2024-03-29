%% Correlation between bursts/spike rate and periodic behaviors (whisking, breathing)
% Simplified version

%% Copy files to analyze to a single folder
% Required file(s):
%     with spike times
%     with whisker position/angle
%     with ephys and video recording times to sync the two
% Optional files:
%     ephys traces
%     video recording

%% Locate data
spikeSortingFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.result.hdf5','*_jrc.mat','*.csv','*_spikesResorted.mat'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% do not include those files:
spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker'}),...
    {spikeSortingFiles.name}));
sessionDir=cd;

dataFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*raw.kwd','*RAW*Ch*.nex','*.ns*'},'UniformOutput', false);
dataFiles=vertcat(dataFiles{~cellfun('isempty',dataFiles)});
% keep those files
TTLFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_TTLs'}),...
    {dataFiles.name}));
dataFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_export'}),...
    {dataFiles.name}));

% try
%     fileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_WhiskerAngle.csv'),...
%         {dirListing.name},'UniformOutput',false))).name;
% catch
%     [fileName,dirName] = uigetfile({'*.mat; *.csv; *.avi','.mat Files; .csv Files; video Files';...
%         '*.*','All Files'},'Whisker angle / Perodical Behavior Data / Video file','C:\Data\Ephys\Behav');
% end
whiskerTrackingFiles=cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.csv'},'UniformOutput', false);
whiskerTrackingFiles=vertcat(whiskerTrackingFiles{~cellfun('isempty',whiskerTrackingFiles)});
whiskerTrackingFiles=whiskerTrackingFiles(cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker'}),...
    {whiskerTrackingFiles.name}));

%% Get video sync data
videoFrameTimeFiles=cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat'},'UniformOutput', false);
videoFrameTimeFiles=vertcat(videoFrameTimeFiles{~cellfun('isempty',videoFrameTimeFiles)});
videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,{'_VideoFrameTimes'}),...
    {videoFrameTimeFiles.name}));

% decide which file to use
recName='vIRt22_1016_5100_50ms1Hz10mW';
spikeFileNum=1; dataFileNum=1; TTLFileNum=1; wTrackNumFile=[1,2]; vFrameFileNum=1;

%also get some info about the recording if possible
% e.g., load rec_info from the spike file:
load(spikeSortingFiles(spikeFileNum).name,'rec_info')
numElectrodes=numel(rec_info.exportedChan); %32 %numel(unique(spikes.preferredElectrode));

%% Load spikes and recording traces
recDir=spikeSortingFiles(spikeFileNum).folder;
recName=spikeSortingFiles(spikeFileNum).name;
cd(recDir)
dataFileIdx=cellfun(@(datF) contains(datF,regexp(recName,'\S+?(?=\.\w+\.\w+$)','match','once')) ,...
    {dataFiles.name});
dataFileName=dataFiles(dataFileIdx).name;
dataFileDir=dataFiles(dataFileIdx).folder;

traces = memmapfile(fullfile(dataFileDir,dataFileName),'Format','int16');
allTraces=double(traces.Data); 
recDuration=int32(length(allTraces)/numElectrodes);
allTraces=reshape(allTraces,[numElectrodes recDuration]);
filterTraces=true;

spikes=LoadSpikeData(recName,traces);

%% load TTLs (are they needed?)
TTLDir=TTLFiles(TTLFileNum).folder;
TTLFileName= TTLFiles(TTLFileNum).name;%[regexp(recName,'\S+?(?=_export)','match','once') '_TTLs.dat'];
fid = fopen(fullfile(TTLDir,TTLFileName), 'r');
TTLTimes = fread(fid,[2,Inf],'int32');
fclose(fid);
TTLs.times=TTLTimes(1,:);
TTLs.samplingRate=1000;

%% Read frame times
videoFrameDir=videoFrameTimeFiles(vFrameFileNum).folder;
videoFrameFileName= videoFrameTimeFiles(vFrameFileNum).name;
if contains(videoFrameFileName,'.dat')
    fid = fopen(fullfile(videoFrameDir,videoFrameFileName), 'r');
    vFrameTimes = fread(fid,[2,Inf],'double');
    fclose(fid);
else
    vFrameTimes=ReadVideoFrameTimes;
    % videoFrameTimes=readVideoTTLData(dirListing);
end

%% Import whisker tracking data (aka "thetas") 
% variable frame rate typically ~500Hz
if numel(wTrackNumFile)==1
    whiskerTrackDir=whiskerTrackingFiles(wTrackNumFile).folder;
    whiskerTrackFileName= whiskerTrackingFiles(wTrackNumFile).name;
    if contains(whiskerTrackFileName,'.csv') % could be npy
        whiskerTrackingData = ImportDLCWhiskerTrackingCSV(fullfile(...
            whiskerTrackDir,whiskerTrackFileName));
    elseif contains(whiskerTrackFileName,'.avi') %video file to extract whisker angle
        whiskerTrackingData=ExtractMultiWhiskerAngle_FFTonContours(fullfile(dirName,fileName));
        whiskerTrackingData=smoothdata(whiskerTrackingData,'rloess',20);
    else
        load([dirName fileName]);
    end
else % contains(fileName,'.csv') % e.g. WhiskerAngle.csv
    whiskerTrackDir=whiskerTrackingFiles(wTrackNumFile(1)).folder;
    whiskerTrackFileName= whiskerTrackingFiles(wTrackNumFile(1)).name;
    mainWhiskerTrackingData=ImportCSVasVector(fullfile(whiskerTrackDir,whiskerTrackFileName));
    whiskerTrackDir=whiskerTrackingFiles(wTrackNumFile(2)).folder;
    whiskerTrackFileName= whiskerTrackingFiles(wTrackNumFile(2)).name;
    multiWhiskerTrackingData=ImportCSVasVector(fullfile(whiskerTrackDir,whiskerTrackFileName));
    % Multiwhisker export from Bonsai has three columns: Orientation Centroid.X Centroid.Y
    % for up to 5 main whiskers (NaN if less)
    whiskerTrackingData=multiWhiskerTrackingData(:,7); %posterior most whisker
    whiskerTrackingData=RadianToDegreesSmoothFill(whiskerTrackingData); 
%     figure; hold on; plot(whiskerTrackingData)
end

%% Recording start time (mostly for OE)
% Processor: Rhythm FPGA Id: 100 subProcessor: 0 start time: 27306000@30000Hz
if exist('rec_info','var') & isfield(rec_info,'recordingStartTime')
    startTime=double(rec_info.recordingStartTime); % 27306000; % !!!
else
    startTime=0;
    %well, make sure time indices are properly aligned
end

%% Add voltage scaling factor and sampling rate
bitResolution=0.195; %for Open Ephys
spikes.waveforms=double(spikes.waveforms.*bitResolution);
if isfield(spikes,'samplingRate')
    samplingRate=unique(spikes.samplingRate);
else
    samplingRate=30000;
end

%% Filter traces if needed
if filterTraces == true
    allTraces=FilterTrace(allTraces,samplingRate);
end

%% Sync ephys and behavior (video) 
% convert video frame times to native recording frame rate (typically 30kHz) if needed
% first remove recording start clock time
vidTimes=vFrameTimes(1,vFrameTimes(2,:)<0)-double(startTime); %when using Paul's OE Basler module
% % also reset spike times if needed
if spikes.times(end) > size(allTraces,2)
    spikes.times=spikes.times-startTime;
end
    TTLs.times=TTLs.times-double(startTime/(samplingRate/1000));
% then sync
if vidTimes(1)>=0
    allTraces=allTraces(:,vidTimes(1):vidTimes(end));
    spikeReIndex=spikes.times>=vidTimes(1) & spikes.times<=vidTimes(end);
    spikes.unitID=spikes.unitID(spikeReIndex);
    spikes.preferredElectrode=spikes.preferredElectrode(spikeReIndex);
    spikes.waveforms=spikes.waveforms(spikeReIndex,:);
    spikes.times=spikes.times(spikeReIndex)-vidTimes(1);
    vidTimes=vidTimes-vidTimes(1);
else % need to cut behavior trace
    %     videoReIndex=...
    %     whiskerTrackingData=...
end

%% Find best units
spikes.unitID=double(spikes.unitID);
spikes.times=double(spikes.times);
% find most frequent units
[unitFreq,uniqueUnitIDs]=hist(spikes.unitID,unique(spikes.unitID));
[unitFreq,freqIdx]=sort(unitFreq','descend');
unitFreq=unitFreq./sum(unitFreq)*100; uniqueUnitIDs=uniqueUnitIDs(freqIdx);
bestUnitsIdx=find(unitFreq>2);
keepUnits=uniqueUnitIDs(bestUnitsIdx); keepUnits=sort(keepUnits(keepUnits~=0));
if isfield(spikes,'preferredElectrode')
    titularChannels = unique(spikes.preferredElectrode(ismember(spikes.unitID,keepUnits)));
end
% keepUnits=[1 2 3];
% titularChannels=[10 10 10];
keepTraces=titularChannels; %14; %[10 14 15];% keepTraces=1:16; %[10 14 15];

%% Keep selected recording trace and spike times,
recordingTraces=allTraces(keepTraces,:); %select the trace to keep
keepUnitsIdx=ismember(spikes.preferredElectrode,keepTraces);
unitID=spikes.unitID(keepUnitsIdx);
preferredElectrode=spikes.preferredElectrode(keepUnitsIdx);
waveForms=spikes.waveforms(keepUnitsIdx,:);
spikeTimes=spikes.times(keepUnitsIdx);
% keepUnits=unique(unitID);
% recordingTrace=cell(length(keepTraces),1);
% for traceNum=1:length(keepTraces)
%     recordingTrace{traceNum}=allTraces(keepTraces(traceNum),:); %select the trace to keep
%     recordingTrace{traceNum}=recordingTrace{traceNum}(TTLs.times(1)*double(samplingRate)/...
%         double(TTLs.samplingRate):end); % cut out trace that occurs before first TTL
% end
% 
% % same for spikes from selected units
% spikeTimes=cell(length(keepUnits),1);
% for clusterNum=1:length(keepUnits)
%     spikeTimes{clusterNum}=double(spikes.times(spikes.unitID==keepUnits(clusterNum))); %figure; plot(spikeTimes)
%     spikeTimes{clusterNum}=spikeTimes{clusterNum}-(double(TTLs.times(1))*double(samplingRate)/double(TTLs.samplingRate));
%     spikeTimes{clusterNum}=spikeTimes{clusterNum}(spikeTimes{clusterNum}>0);
% end
% % rs_spikeTimes=double(source_spikes.times(source_spikes.unitID==2)); %figure; plot(spikeTimes)

figure('Color','white');
for traceNum=1:length(keepTraces)
%     subplot(length(keepTraces),1,traceNum); 
    figure('Color','white'); hold on;
    plot(recordingTraces(traceNum,:));
    correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
    for clusterNum=1:numel(correspondingUnits)
        plot(spikeTimes(unitID==correspondingUnits(clusterNum)),...
            ones(numel(spikeTimes(unitID==correspondingUnits(clusterNum))),1)*-300,'*');
        %     plot(rs_spikeTimes,ones(size(rs_spikeTimes,1),1)*-300,'d')
    end
    midRec=round(size(recordingTraces,2)/2);
    set(gca,'ylim',[-500 500],'xlim',[midRec-samplingRate midRec+samplingRate]);
end

% foo=resample(periodBehavData,30,1); foo=foo(1:length(allTraces(channelNum,:)));
% plot(foo*1000-200)


% params.Fs=30000;params.fpass=[0 25];params.tapers=[2 3];params.pad=1;params.err=[2 0.05];params.trialave=0;
% [C,phi,S12,S1,S2,f]=coherencyc(foo',double(allTraces(channelNum,:))',params);
% figure; subplot(311); plot(f,C);subplot(312); plot(f,10*log10(S1));subplot(313); plot(f,10*log10(S2))

%% Downsample trace to 1ms (pointless, as it will loose the spikes)
% Fs=1000;
% % foo=resample(double(recordingTrace),Fs,double(samplingRate));
% dsrecordingTrace=decimate(double(recordingTrace),double(samplingRate)/Fs,'fir');

%% Bin spike counts in 1ms bins
% with Chronux' binning function
% foo=binspikes(spikeTimes/double(samplingRate),Fs);
% foo=[zeros(round(spikeTimes(1)/double(samplingRate)*Fs)-1,1);foo]; %need to padd with zeroes
% With home-made function. Same result, but takes care of the padding
binSize=1;
spikeRasters_ms=zeros(numel(keepUnits),ceil(size(recordingTraces,2)/samplingRate*1000));
for clusterNum=1:length(keepUnits)
    unitIdx=unitID==keepUnits(clusterNum);
    lengthUnitTimeArray=ceil(spikeTimes(find(unitIdx,1,'last'))/samplingRate*1000);
    spikeRasters_ms(clusterNum,1:lengthUnitTimeArray)=DownSampleToMilliseconds(...
        spikeTimes(unitIdx),binSize,samplingRate);
end

% figure; hold on
% % plot(dsrecordingTrace)
% plot(find(binSpikeTimes),ones(length(find(binSpikeTimes)),1)*-250,'r*')
% plot(find(foo),ones(length(find(foo)),1)*-200,'g*')

%% Compute sdfs
SDFs_ms=nan(length(keepUnits), ceil(size(recordingTraces,2)/samplingRate*1000));
for clusterNum=1:length(keepUnits)
    SDFs_ms(clusterNum,:)=GaussConv(spikeRasters_ms(clusterNum,:),5)*1000;
end
% figure; hold on
% plot(SDFs{1})
% plot(find(binSpikeTimes{1}),ones(length(find(binSpikeTimes{1})),1)*-10,'r*')

%% Compute raster indices
[rasterYInd_ms, rasterXInd_ms]=deal(cell(length(keepUnits),1));
for clusterNum=1:length(keepUnits)
    [rasterYInd_ms{clusterNum}, rasterXInd_ms{clusterNum}] =...
        ind2sub(size(spikeRasters_ms(clusterNum,:)),find(spikeRasters_ms(clusterNum,:))); %find row and column coordinates of spikes
end
% rasters=[indx indy;indx indy+1];

%% Create array with angle values and time points
periodBehavData=[whiskerTrackingData,vidTimes'];
% periodBehavData=[whiskerTrackingData(videoFrameTimes.TTLFrames(1):...
%                   size(videoFrameTimes.frameTime_ms,1)),... %Trace
%     videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1):end)-...
%     videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1))+1]; % Time points
% figure; hold on
% plot(periodBehavData(:,2),periodBehavData(:,1))

%resample to 1ms precision
periodBehavData_ms=periodBehavData;
periodBehavData_ms(:,2)=periodBehavData_ms(:,2)/samplingRate*1000;
% [periodBehavData(:,1),periodBehavData(:,2)] = resample(periodBehavData(:,1),periodBehavData(:,2),'pchip');
periodBehavData_ms=interp1(periodBehavData_ms(:,2),periodBehavData_ms(:,1),...
    periodBehavData_ms(1,2):periodBehavData_ms(end,2));
% figure; plot(periodBehavData);

%% Plot behavior data and find a period with whisking bouts
% no need to keep periods with no whisking
% figure; plot(periodBehavData); % select data point and export cursor info
% whiskingPeriod=1:cursor_info.Position(1); %in ms
peakWhisking_ms=diff(cummax(abs(diff(periodBehavData_ms))));
peakWhiskingIdx=find(peakWhisking_ms==max(peakWhisking_ms));
whiskingPeriod=peakWhiskingIdx-5000:peakWhiskingIdx+4999; %in ms

%% Filter periodic behavior traces into low-pass and bandpassed versions
LP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,0.3,'low'); %set-point
% figure; hold on
% plot(periodBehavData); plot(LP_periodBehavData,'LineWidth',2)

BP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,[0.3 30],'bandpass'); %whisking
% figure; hold on
% plot(periodBehavData-mean(periodBehavData)); plot(BP_periodBehavData,'LineWidth',1)

HP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,0.3,'high')'; %whisking
plot(HP_periodBehavData_ms,'LineWidth',1)

% make sure behavior and spike traces have same length
if size(SDFs_ms,2)~=numel(periodBehavData_ms)
    % check what 
end

if exist('whiskingPeriod','var')
    %     sdf=sdf(whiskingPeriod);
    %     LP_periodBehavData=LP_periodBehavData(whiskingPeriod);
    %     BP_periodBehavData=BP_periodBehavData(whiskingPeriod);
end

% figure; hold on
% plot(SDFs{1}); plot(BP_periodBehavData/min(BP_periodBehavData)*...
%     max(SDFs{1}) + max(SDFs{1}))

%% Hilbert transform
HTBP_periodBehavData_ms=hilbert(BP_periodBehavData_ms);
% figure; plot(imag(HTBP_periodBehavData));
whiskingPhase_ms=angle(HTBP_periodBehavData_ms);
% figure; hold on
% plot(SDFs); plot(whiskingPhase*10 + max(SDFs))

%% Find instantaneous frequency
Nfft = 1024;
% [Pxx,f] = pwelch(BP_periodBehavData,gausswin(Nfft),Nfft/2,Nfft,1000);
% figure; plot(f,Pxx); ylabel('PSD'); xlabel('Frequency (Hz)'); grid on;
[~,sgFreq,sgTime,sgPower] = spectrogram(BP_periodBehavData_ms,gausswin(Nfft),Nfft/2,Nfft,1);
instantFreq_ms = medfreq(sgPower,sgFreq);
% figure; plot(sgTime,round(instantFreq*1000),'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%

%% plot traces together
%%%%%%%%%%%%%%%%%%%%%%%
displayUnits=1; %4; %1:length(keepUnits); %default: 1:length(keepUnits) [1 4 5]
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
    % subplot('Position',[0 1-(1/(4)*(clusterNum+3)) 1 1/(4)]);
    winIndx=rasterXInd_ms{displayUnits(clusterNum)}...
        (rasterXInd_ms{displayUnits(clusterNum)}>=dispWin(1) & rasterXInd_ms{displayUnits(clusterNum)}<=dispWin(end))-timeAxis(1);
    winIndy=rasterYInd_ms{displayUnits(clusterNum)}...
        (rasterXInd_ms{displayUnits(clusterNum)}>=dispWin(1) & rasterXInd_ms{displayUnits(clusterNum)}<=dispWin(end));
    plot([winIndx;winIndx],[winIndy+clusterNum-1;winIndy+clusterNum],...
        'color',cmap(clusterNum,:),'LineStyle','-','LineWidth',1.8); % plot rasters% plot rasters
    plot(BP_periodBehavData_ms(timeAxis)/max(abs(BP_periodBehavData_ms(timeAxis)))+...
        clusterNum+0.5,'color','k','LineWidth',0.8);
    plot(whiskingPhase_ms(timeAxis)/max(abs(whiskingPhase_ms(timeAxis)))+...
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


%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[602   537   560   420]); hold on
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms{displayUnits(clusterNum)},BP_periodBehavData_ms,500,'coeff');
    %     figure('position',[602   537   560   420]);
    plot(lag,acor,'LineWidth',2); %xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
end
legend();
xlabel('Lag (ms)');set(gca,'ylim',[-0.5 0.5])
title({['Cross correlation for vIRt unit ' num2str(keepUnits(displayUnits(clusterNum)))];...
    'Spike density function vs. Whisking angle'})

for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms{displayUnits(clusterNum)},LP_periodBehavData_ms,1000,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Set point'})
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms{displayUnits(clusterNum)},whiskingPhase_ms,100,'coeff');
    figure('position',[602   537   560   420]);
    plot(lag,acor); xlabel('Lag (ms)'); set(gca,'ylim',[-1 1])
    title({['Cross correlation for unit' num2str(displayUnits(clusterNum))];'Spike density function vs. Whisking phase'})
    % hold on; [acor,lag] = xcorr(sdf,BP_periodBehavData,100,'unbiased'); plot(lag,acor);
end
for clusterNum=1:length(displayUnits)
    [acor,lag] = xcorr(SDFs_ms{displayUnits(clusterNum)},periodBehavData_ms,2500,'coeff');
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
        spikeTimes{displayUnits(clusterNum)}/double(samplingRate),params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(BP_periodBehavData_ms',...
        SDFs_ms{displayUnits(clusterNum)}',params);
    figure; subplot(311); plot(sgFreq,C);subplot(312); plot(sgFreq,10*log10(S1));subplot(313); plot(sgFreq,10*log10(S2))
    C(S1==max(S1))
end
for clusterNum=1:length(displayUnits)
    [C,phi,S12,S1,S2,sgFreq,zerosp,confC,phierr,Cerr]=coherencycpb(whiskingPhase_ms',...
        SDFs_ms{displayUnits(clusterNum)}',params);
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

%% Polar plot
%
% unitSDF=SDFs{bestUnit};
% unitSDFExcerpt=unitSDF(int32([midRec-100:midRec+100]*1000));
% get atda
for unitNum=1:size(spikeRasters_ms,1)
unitSpikes=spikeRasters_ms(unitNum,:);
vibrissaAngle=BP_periodBehavData_ms;
unitSpikes=unitSpikes(1:length(vibrissaAngle));
numSegments=round(length(vibrissaAngle)/60000);
figure;
for pplotNum=1:numSegments-2
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
    whiskingPhase_ms=angle(HTvibrissaAngle);
    
    spikeOnWPhase=whiskingPhase_ms(logical(unitSpikesExcerpt));
    subplot(3,2,pplotNum);
    polarhistogram(spikeOnWPhase,72,'Displaystyle','stairs',...
        'Normalization','count')
end
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




