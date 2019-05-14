function [ephys,behav]=NeuronBehaviorCorrelation_GatherData

% Gather data for analysis of correlation between bursts/spike rate and periodic behaviors (whisking, breathing)
% Simplified version
startingDir=cd;
%% Place files to analyze in current folder
% Required file(s):
%     spike times
%     whisker position/angle
%     ephys and video recording times to sync the two
% Optional files:
%     ephys traces
%     video recording

%% Locate data
% Spikes
spikeSortingFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*_spikes.mat','*.result.hdf5','*_rez.mat','*_res.mat','*_jrc.mat','*.csv','*_spikesResorted.mat'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% do not include those files:
spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker','trial'}),...
    {spikeSortingFiles.name}));
sessionDir=cd;

dataFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*.bin','*raw.kwd','*RAW*Ch*.nex','*.ns*'},'UniformOutput', false);
dataFiles=vertcat(dataFiles{~cellfun('isempty',dataFiles)});
% keep those files
% TTLFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_trialTTLs'}),...
%     {dataFiles.name})); %not used here. Used for Phototag plots
dataFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_export';'_traces'}),...
    {dataFiles.name}));

% Whisker tracking files
whiskerTrackingFiles=cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.csv','whiskerTrackingData'},'UniformOutput', false);
whiskerTrackingFiles=vertcat(whiskerTrackingFiles{~cellfun('isempty',whiskerTrackingFiles)});
whiskerTrackingFiles=whiskerTrackingFiles(cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker'}),...
    {whiskerTrackingFiles.name}));

TTLFiles=cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*TTLOnset.csv','whiskerTrackingData'},'UniformOutput', false);
TTLFiles=vertcat(TTLFiles{~cellfun('isempty',TTLFiles)});

%% Get video sync data
videoFrameTimeFiles=cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*.csv'},'UniformOutput', false);
videoFrameTimeFiles=vertcat(videoFrameTimeFiles{~cellfun('isempty',videoFrameTimeFiles)});
videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,{'_VideoFrameTimes','vSync'}),...
    {videoFrameTimeFiles.name}));

% decide which file to use
spikeFileNum=1; wTrackNumFile=1;
if sum(cellfun(@(flnm) contains(flnm,'vSync'),{videoFrameTimeFiles.name}))
    videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,'vSync'),...
        {videoFrameTimeFiles.name}));
else
    videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,'_VideoFrameTimes'),...
        {videoFrameTimeFiles.name})); 
end

%% Load spikes and recording traces
spikeSortFromExportIdx=cellfun(@(flnm) contains(flnm,'_export'),{spikeSortingFiles.name});
if sum(spikeSortFromExportIdx)
    % prefer offline spike sorting results
    spikeSortingFiles=spikeSortingFiles(spikeSortFromExportIdx);
end
recDir=spikeSortingFiles(spikeFileNum).folder;
recName=spikeSortingFiles(spikeFileNum).name;
cd(recDir)
%also get some info about the recording if possible
% e.g., load rec_info from the spike file:
if logical(sum(cellfun(@(x) contains(x,'rec_info'),...
        who(matfile(spikeSortingFiles(spikeFileNum).name)))))
    load(spikeSortingFiles(spikeFileNum).name,'rec_info')
elseif ~exist('rec_info','var') && exist([dataFiles(spikeFileNum).name(1:end-10) 'recInfo.mat'],'file')
    load([dataFiles(spikeFileNum).name(1:end-10) 'recInfo']);
else
    recInfo=rec_info; clear rec_info; %load('vIRt12-2018-01-18-15-53-03-Acute-KAionto-5100-cutW_nopp_recInfo.mat', 'rec_info')
end
if isfield(recInfo,'exportedChan')
    numElectrodes=numel(recInfo.exportedChan);
elseif isfield(recInfo,'numRecChan')
    numElectrodes=numel(recInfo.numRecChan); %32 %numel(unique(spikes.preferredElectrode));
else
    %Assuming 32 %Could also check numel(unique(spikes.preferredElectrode));
    numElectrodes=32;
end
if numElectrodes==35
    numElectrodes=32; %AUX channels removed
end

dataFileIdx=cellfun(@(datF) contains(datF,regexp(recName,'\S+?(?=\.\w+\.\w+$)','match','once')) ,...
    {dataFiles.name});
dataFileName=dataFiles(dataFileIdx).name;
dataFileDir=dataFiles(dataFileIdx).folder;

traces = memmapfile(fullfile(dataFileDir,dataFileName),'Format','int16');
allTraces=double(traces.Data);
recDuration=int64(length(allTraces)/numElectrodes);
try
    allTraces=reshape(allTraces,[numElectrodes recDuration]);
catch
    allTraces=reshape(allTraces',[recDuration numElectrodes]);
end
filterTraces=false;

spikes=LoadSpikeData(recName,traces);
if isfield(spikes,'samplingRate') && isempty(spikes.samplingRate)
    spikes.samplingRate=recInfo.samplingRate;
    spikes.bitResolution=recInfo.bitResolution;
end

%% Read frame times
videoFrameDir=videoFrameTimeFiles(1).folder;
videoFrameFileName= videoFrameTimeFiles(1).name;
if contains(videoFrameFileName,'vSyncTTLs.dat')
    fid = fopen(fullfile(videoFrameDir,videoFrameFileName), 'r');
    vFrameTimes = fread(fid,[2,Inf],'int32');
    fclose(fid);
elseif contains(videoFrameFileName,'VideoFrameTimes.dat')
    fid = fopen(fullfile(videoFrameDir,videoFrameFileName), 'r');
    vFrameTimes = fread(fid,[2,Inf],'double');
    fclose(fid);
else % csv file from Bonsai
    vFrameTimes=ReadVideoFrameTimes;
end

%% Import whisker tracking data 
startDirFiles=dir(startingDir);
isWTData=cellfun(@(x) contains(x,'whiskerTrackingData'), {startDirFiles.name});
if sum(isWTData)
    load(fullfile(startDirFiles(isWTData).folder,startDirFiles(isWTData).name));
else
    % variable frame rate typically ~500Hz
    whiskerTrackDir=whiskerTrackingFiles(1).folder;
    whiskerTrackFileName= whiskerTrackingFiles(1).name;
    if contains(whiskerTrackFileName,'.csv') % e.g. WhiskerAngle.csv
        if contains(whiskerTrackFileName,'.npy')
        elseif contains(whiskerTrackFileName,'.csv')
            if contains(whiskerTrackFileName,'DeepCut') || contains(whiskerTrackFileName,'DLC')
                whiskerTrackingData = ImportDLCWhiskerTrackingCSV(fullfile(...
                    whiskerTrackDir,whiskerTrackFileName));
            else %assuming from Bonsai
                %             depending on version, export from Bonsai has either
                %               one column: Orientation
                %               three columns:  Centroid.X Centroid.Y Orientation
                %               6 times three columns: Base, Centroid.X and Centroid.Y for each whisker
                if numel(wTrackNumFile)==1
                    if contains(whiskerTrackFileName,'BaseCentroid')
                        whiskerTrackingData=readtable(fullfile(whiskerTrackDir,whiskerTrackFileName));
                        whiskerTrackingData=ContinuityWhiskerID(whiskerTrackingData);
                    else
                        delimiter=' ';hasHeader=false;
                        whiskerTrackingData=ImportCSVasVector(...
                            fullfile(whiskerTrackDir,whiskerTrackFileName),delimiter,hasHeader);
                        if size(whiskerTrackingData,2)>1
                            whiskerTrackingData=whiskerTrackingData(:,1:2);
                        end
                    end
                else
                    whiskerTrackDir=whiskerTrackingFiles(wTrackNumFile(2)).folder;
                    whiskerTrackFileName= whiskerTrackingFiles(wTrackNumFile(2)).name;
                    multiWhiskerTrackingData=ImportCSVasVector(fullfile(whiskerTrackDir,whiskerTrackFileName));
                    % Multiwhiskerfor up to 5 main whiskers (NaN if less)
                    whiskerTrackingData=multiWhiskerTrackingData(:,7); %posterior most whisker
                end
                whiskerTrackingData=WhiskerAngleSmoothFill(whiskerTrackingData); %(:,1),whiskerTrackingData(:,2));
            end
        elseif contains(whiskerTrackFileName,'.avi') %video file to extract whisker angle
            whiskerTrackingData=ExtractMultiWhiskerAngle_FFTonContours(fullfile(dirName,fileName));
            whiskerTrackingData=smoothdata(whiskerTrackingData,'rloess',20);
        else
            load([dirName fileName]);
        end
    end
    cd(startingDir)
    save('whiskerTrackingData','whiskerTrackingData');
end
%% Recording start time (mostly for OE)
if exist('recInfo','var') & isfield(recInfo,'recordingStartTime')
    startTime=double(recInfo.recordingStartTime); 
else
    startTime=0;
    %well, make sure time indices are properly aligned
end

%% Add voltage scaling factor and sampling rate
if ~isfield(spikes,'bitResolution') | (isfield(spikes,'bitResolution') && isempty(spikes.bitResolution))
    spikes.bitResolution=0.195; %for Open Ephys
end
spikes.waveforms=double(spikes.waveforms.*spikes.bitResolution);
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
% convert video frame times to native ephys recording frame rate (typically 30kHz) if needed
% first remove recording start clock time
vidTimes=vFrameTimes(1,vFrameTimes(2,:)<0)-double(startTime); %when using Paul's OE Basler module

% % also reset spike times if needed
if spikes.times(end) > size(allTraces,2)
    spikes.times=spikes.times-startTime;
end
% then sync
if vidTimes(1)>=0
    allTraces=allTraces(:,vidTimes(1):vidTimes(end));
    spikeReIndex=spikes.times>=vidTimes(1) & spikes.times<=vidTimes(end);
    spikes.unitID=spikes.unitID(spikeReIndex);
    try
        spikes.preferredElectrode=spikes.preferredElectrode(spikeReIndex);
    catch %might referenced by cluster or electrodes
    end
    try
        spikes.waveforms=spikes.waveforms(spikeReIndex,:);
    catch
        %?
    end
    spikes.times=spikes.times(spikeReIndex)-vidTimes(1);
    vidTimes=vidTimes-vidTimes(1);
else % need to cut behavior trace
    %     videoReIndex=...
    %     whiskerTrackingData=...
end

% Adjust behavior array if necessary
if numel(whiskerTrackingData(1,:))~=numel(vidTimes)
    vidTimes=linspace(1,double(vidTimes(end)+1),size(whiskerTrackingData,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% group data in structure
ephys=struct('traces',allTraces,'spikes',spikes,...
    'samplingRate',samplingRate,'recName',recName);

behav=struct('whiskerTrackingData',whiskerTrackingData,'vidTimes',vidTimes,...
    'vFrameTimes',vFrameTimes,'startTime',startTime);
 
cd(startingDir);




