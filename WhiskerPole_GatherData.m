function [spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
    SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
    templates,templateIdx,stimType,stimTimes_ms,videoFile] = WhiskerPole_GatherData

% [recordingTraces,spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
%     SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
%     BP_periodBehavData_ms,HP_periodBehavData_ms,LP_periodBehavData_ms,...
%     HTBP_periodBehavData_ms, peakWhisking_ms,periodBehavData_ms,...
%     whiskingPhase_ms,instantFreq_ms,sgFreq,sgTime,sgPower,recName] = NeuronBehaviorCorrelation_GatherData

% Gather data for analysis of whisker touch behaviors and spike responses

%% Place files to analyze in current folder
% Required file(s):
%     spike times (from spike sorting)
%     video frame times (see BatchExport)
%     whisker stimulation times (from Bonsai)
% Optional files:
%     ephys traces
%     video recording

% rootDir=cd; 
%% Locate data
spikeSortingFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.result.hdf5','*_jrc.mat','*_spikesResorted.mat'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'',''}),...
%     {spikeSortingFiles.name}));

vSyncFiles = dir([cd filesep '**' filesep '*.dat']);
vSyncFiles=vSyncFiles(cellfun(@(flnm) contains(flnm,{'_VideoFrameTimes'}),...
    {vSyncFiles.name}));

whiskerStimFiles= dir([cd filesep '**' filesep '*.csv']);
whiskerStimFiles=whiskerStimFiles(cellfun(@(flnm) contains(flnm,{'_ROIWhiskerContact'}),...
    {whiskerStimFiles.name}));

videoFile=[];
%also get some info about the recording if possible
% e.g., load rec_info from the spike file:
% load(spikeSortingFiles(spikeFileNum).name,'rec_info')
% numElectrodes=numel(rec_info.exportedChan); %32 %numel(unique(spikes.preferredElectrode));

%% Load spikes and recording traces
% recDir=spikeSortingFiles(spikeFileNum).folder;
% cd(recDir)
% dataFileIdx=cellfun(@(datF) contains(datF,regexp(recName,'\S+?(?=\.\w+\.\w+$)','match','once')) ,...
%     {dataFiles.name});
% dataFileName=dataFiles(dataFileIdx).name;
% dataFileDir=dataFiles(dataFileIdx).folder;
% 
% traces = memmapfile(fullfile(dataFileDir,dataFileName),'Format','int16');
% allTraces=double(traces.Data); 
% recDuration=int32(length(allTraces)/numElectrodes);
% allTraces=reshape(allTraces,[numElectrodes recDuration]);
% filterTraces=true;

cd(spikeSortingFiles.folder);
spikes=LoadSpikeData(spikeSortingFiles.name); %,traces);

%% Add voltage scaling factor and sampling rate
% spikes.waveforms=double(spikes.waveforms.*bitResolution);
if isfield(spikes,'samplingRate')
    samplingRate=unique(spikes.samplingRate);
else
    samplingRate=30000;
end

%% load frame times
vsyncFID=fopen(vSyncFiles.name, 'r');
frameTimes = fread(vsyncFID,[2,Inf],'double');
fclose(vsyncFID);

%% load whisker stim times
stimulationTypes={'ManStim','Pole','Texture'};
stimType=stimulationTypes{cellfun(@(stimT) contains(whiskerStimFiles.name,stimT),...
        stimulationTypes)};
    
whiskerStim_FrameIdx=GetFrame_WhiskerTouch(whiskerStimFiles.name,whiskerStimFiles.folder);

%% get stimulus time
stimTimes_ms=cell(size(whiskerStim_FrameIdx,2),1);
for ROInum= 1:size(whiskerStim_FrameIdx,2)
    % remove periods that are too short
    stimPeriods=bwconncomp(whiskerStim_FrameIdx(:,ROInum));
    stimFrameIdx=cellfun(@(x) x(1), {stimPeriods.PixelIdxList{1,cellfun(@numel,stimPeriods.PixelIdxList)>3}});
    % then get corresponding frame time
    stimTimes_ms{ROInum}=frameTimes(1,stimFrameIdx)/samplingRate*1000;
end

%% Filter traces if needed
% if filterTraces == true
%     allTraces=FilterTrace(allTraces,samplingRate);
% end

%% Find best units
spikes.unitID=double(spikes.unitID);
spikes.times=double(spikes.times);
% find most frequent units
[unitFreq,uniqueUnitIDs]=hist(spikes.unitID,unique(spikes.unitID));
[unitFreq,freqIdx]=sort(unitFreq','descend');
unitFreq=unitFreq./sum(unitFreq)*100; uniqueUnitIDs=uniqueUnitIDs(freqIdx);
bestUnitsIdx=find(unitFreq>0.8);
keepUnits=uniqueUnitIDs(bestUnitsIdx); keepUnits=sort(keepUnits(keepUnits~=0));
if isfield(spikes,'preferredElectrode')
    titularChannels = unique(spikes.preferredElectrode(ismember(spikes.unitID,keepUnits)));
end
% keepUnits=[1 2 3];
% titularChannels=[10 10 10];
keepTraces=titularChannels; %14; %[10 14 15];% keepTraces=1:16; %[10 14 15];

%% Keep selected recording trace and spike times,
% recordingTraces=allTraces(keepTraces,:); %select the trace to keep
keepUnitsIdx=ismember(spikes.preferredElectrode,keepTraces);
unitID=spikes.unitID(keepUnitsIdx);
preferredElectrode=spikes.preferredElectrode(keepUnitsIdx);
waveForms=spikes.waveforms; %(:,:,keepUnitsIdx);
templates=spikes.templates(keepUnitsIdx);
templateIdx=spikes.templatesIdx;
spikeTimes=spikes.times(keepUnitsIdx);

%% Bin spike counts in 1ms bins
% with Chronux' binning function
% foo=binspikes(spikeTimes/double(samplingRate),Fs);
% foo=[zeros(round(spikeTimes(1)/double(samplingRate)*Fs)-1,1);foo]; %need to padd with zeroes
% With home-made function. Same result, but takes care of the padding
binSize=1;
spikeRasters_ms=zeros(numel(keepUnits),ceil(max(spikeTimes)/samplingRate*1000)); %size(recordingTraces,2)
for clusterNum=1:length(keepUnits)
    unitIdx=unitID==keepUnits(clusterNum);
    lengthUnitTimeArray=ceil(spikeTimes(find(unitIdx,1,'last'))/samplingRate*1000);
    spikeRasters_ms(clusterNum,1:lengthUnitTimeArray)=DownSampleToMilliseconds(...
        spikeTimes(unitIdx),binSize,samplingRate);
end


%% Compute sdfs
SDFs_ms=nan(length(keepUnits), ceil(max(spikeTimes)/samplingRate*1000));%size(recordingTraces,2)
for clusterNum=1:length(keepUnits)
    SDFs_ms(clusterNum,:)=GaussConv(spikeRasters_ms(clusterNum,:),5)*1000;
end
% figure; hold on
% plot(SDFs_ms(1,:))
% plot(find(binSpikeTimes{1}),ones(length(find(binSpikeTimes{1})),1)*-10,'r*')

%% Compute raster indices
[rasterYInd_ms, rasterXInd_ms]=deal(cell(length(keepUnits),1));
for clusterNum=1:length(keepUnits)
    [rasterYInd_ms{clusterNum}, rasterXInd_ms{clusterNum}] =...
        ind2sub(size(spikeRasters_ms(clusterNum,:)),find(spikeRasters_ms(clusterNum,:))); %find row and column coordinates of spikes
end
% rasters=[indx indy;indx indy+1];





