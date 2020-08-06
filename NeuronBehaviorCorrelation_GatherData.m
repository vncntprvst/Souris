function [allDataFileNames,targetDir]=NeuronBehaviorCorrelation_GatherData(vararg)

% Gather data for analysis of correlation between bursts/spike rate 
% and periodic behaviors (whisking, breathing)

%% Directory structure assumed be
% > Recording session   (may contain "raw" data)
%     > Spike Sorting Folder
%         > Recording #1
%         > Recording #2
%         ...
%     > Whisker Tracking Folder (all session data processed together)
%     > Analysis Folder (will be created of does not exist yet)

%% Place files to analyze in current folder
% If current folder isn't the Analysis folder, files will be copied there
% Assuming to be in a given recording folder (e.g. Recording #1)
% Required file(s):
%     spike times
%     whisker position/angle
%     ephys and video recording times to sync the two
% Optional files:
%     ephys traces
%     video recording

if nargin~=0
    startingDir=vararg{1};
else
    startingDir=cd;
end
directoryHierarchy=regexp(startingDir,['\S+?(?=\' filesep ')'],'match');

%%%%%%%%%%%%%%%%%
%% Locate data %%
%%%%%%%%%%%%%%%%%
% use dir([cd filesep '**' filesep fileFormat]) search strategy in case
% files are in subdirectories

%% Spikes data file
spikeSortingFiles = cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*_spikes.mat','*.result.hdf5','*_rez.mat','*_res.mat','*_jrc.mat','*.csv','*_spikesResorted.mat'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% do not include those files:
spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker','trial'}),...
    {spikeSortingFiles.name}));

%% Spike waveforms
spikeWaveformFiles = cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*_filt.jrc'},'UniformOutput', false);
spikeWaveformFiles=vertcat(spikeWaveformFiles{~cellfun('isempty',spikeWaveformFiles)});

%% Ephys recording data files
ephysTraceFiles = cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*.dat','*.bin','*raw.kwd','*RAW*Ch*.nex','*.ns*'},'UniformOutput', false);
ephysTraceFiles=vertcat(ephysTraceFiles{~cellfun('isempty',ephysTraceFiles)});
ephysTraceFiles=ephysTraceFiles(cellfun(@(flnm) contains(flnm,{'_export';'_traces'}),...
    {ephysTraceFiles.name}));

%% Recording info
infoFiles = cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*info*'},'UniformOutput', false);
infoFiles=vertcat(infoFiles{~cellfun('isempty',infoFiles)});

%% Probe file (may be needed for channel map, etc)
probeFiles = cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*prb','*Probe*'},'UniformOutput', false);
probeFiles =vertcat(probeFiles{~cellfun('isempty',probeFiles )});

%% Whisker tracking files
% Typically exported from ConvertWhiskerData as *_wMeasurements.mat files. If not there, run it.
whiskerFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*_wMeasurements.mat'},'UniformOutput', false); %'*.csv','whiskerTrackingData',
if isempty(whiskerFiles{:})
    % check other format (e.g., from DLC)
    whiskerFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
        {'*.csv','whiskerTrackingData'},'UniformOutput', false);
    whiskerFiles=vertcat(whiskerFiles{~cellfun('isempty',whiskerFiles)});
    whiskerFiles=whiskerFiles(~cellfun(@(flnm) contains(flnm,{'trial';'analysis'}),...
    {whiskerFiles.name}));
    if ~isempty(whiskerFiles)
        ConvertWhiskerData;
        whiskerFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
            {'*_wMeasurements.mat'},'UniformOutput', false); %'*.csv','whiskerTrackingData',
    else
        %% Ask location
        disp('no whisker tracking file')
        [whiskerFiles,whiskerFilesPath] = uigetfile({'*.mat';'*.*'},...
            'Select the whisker tracking file',startingDir,'MultiSelect','on');
        if ~isempty(whiskerFiles)
            if ~iscell(whiskerFiles); whiskerFiles={whiskerFiles}; end
            whiskerFiles=cellfun(@(fName) fullfile(whiskerFilesPath,fName), whiskerFiles);
        else
            return
        end
    end
end
whiskerFiles=vertcat(whiskerFiles{~cellfun('isempty',whiskerFiles)});
% % velocity 
% whiskerVelocityFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
%     {'*whiskervelocity.mat'},'UniformOutput', false);
% whiskerVelocityFiles=vertcat(whiskerVelocityFiles{~cellfun('isempty',whiskerVelocityFiles)});
% % phase
% whiskerPhaseFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
%     {'*whiskerphase.mat'},'UniformOutput', false);
% whiskerPhaseFiles=vertcat(whiskerPhaseFiles{~cellfun('isempty',whiskerPhaseFiles)});

%% TTL files (other than sync to video, e.g., laser)
TTLFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*TTLOnset.csv','whiskerTrackingData','*trialTS.csv','*trial.mat','*_TTLs.dat'},'UniformOutput', false);
TTLFiles=vertcat(TTLFiles{~cellfun('isempty',TTLFiles)});

%% Video sync data
videoFrameTimeFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*.dat','*.csv'},'UniformOutput', false);
videoFrameTimeFiles=vertcat(videoFrameTimeFiles{~cellfun('isempty',videoFrameTimeFiles)});
videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,{'_VideoFrameTimes','vSync'}),...
    {videoFrameTimeFiles.name}));

%% Video sync files info
videoSyncInfoFiles=cellfun(@(fileFormat) dir([startingDir filesep '**' filesep fileFormat]),...
    {'*_VideoSyncFilesLoc*';'*_WhiskerSyncFilesLoc*'},'UniformOutput', false);
videoSyncInfoFiles=vertcat(videoSyncInfoFiles{~cellfun('isempty',videoSyncInfoFiles)});

%% Decide which file to use
% Keep only the most recent data file
allDataFiles=struct('spikeSortingFiles',spikeSortingFiles,...
    'ephysTraceFiles',ephysTraceFiles,...
    'spikeWaveformFiles',spikeWaveformFiles,...
    'infoFiles',infoFiles,...
    'probeFiles',probeFiles,...
    'TTLFiles',TTLFiles,...
    'videoFrameTimeFiles',videoFrameTimeFiles,...
    'videoSyncInfoFiles',videoSyncInfoFiles,...
    'whiskerFiles',whiskerFiles);
adf_fn=fields(allDataFiles);
for dataFileNum=1:numel(adf_fn)
    [~,dateSort]=sort({allDataFiles.(adf_fn{dataFileNum}).date});
    allDataFiles.(adf_fn{dataFileNum})=allDataFiles.(adf_fn{dataFileNum})(dateSort(dateSort==max(dateSort)));
    allDataFiles.(adf_fn{dataFileNum}).exportname=allDataFiles.(adf_fn{dataFileNum}).name;
end
% mark spike and recording data as such
allDataFiles.(adf_fn{1}).exportname=...
    [allDataFiles.(adf_fn{1}).exportname(1:end-4) ...
    '.spikes'  allDataFiles.(adf_fn{1}).exportname(end-3:end)];
allDataFiles.(adf_fn{2}).exportname=...
    [allDataFiles.(adf_fn{2}).exportname(1:end-4) ...
    '.rec'  allDataFiles.(adf_fn{2}).exportname(end-3:end)];

%% copy processed files to Analysis folder
% try creating folder within Analysis folder: find common file part
allDataFileNames=cellfun(@(fName) getfield(allDataFiles,{1},fName,{1},'exportname'),...
    adf_fn,'UniformOutput', false);
commonStr = GetCommonString(allDataFileNames([1:4,6:9])); %exclude probe file name
if ~isempty(commonStr)
    commonStr=regexprep(commonStr,'[^a-zA-Z0-9]+$','');
end
targetDir=fullfile(directoryHierarchy{1:end-1},'Analysis',commonStr);
if ~exist(targetDir,'dir')
    mkdir(targetDir);
end

for dataFileNum=1:numel(adf_fn)
    copyfile(fullfile(allDataFiles.(adf_fn{dataFileNum}).folder,...
        allDataFiles.(adf_fn{dataFileNum}).name),...
        fullfile(targetDir,... %directoryHierarchy{1:end-1},'Analysis',...
        allDataFiles.(adf_fn{dataFileNum}).exportname));
end







