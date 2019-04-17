function [recordingTraces,spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
    SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
    BP_periodBehavData_ms,HP_periodBehavData_ms,LP_periodBehavData_ms,...
    HTBP_periodBehavData_ms, peakWhisking_ms,periodBehavData_ms,...
    whiskingPhase_ms,instantFreq_ms,sgFreq,sgTime,sgPower,recName,vidTimes_ms] = NeuronBehaviorCorrelation_GatherData

% Gather data for analysis of correlation between bursts/spike rate and periodic behaviors (whisking, breathing)
% Simplified version

%% Place files to analyze in current folder
% Required file(s):
%     spike times
%     whisker position/angle
%     ephys and video recording times to sync the two
% Optional files:
%     ephys traces
%     video recording

%% Locate data
spikeSortingFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*_spikes.mat','*.result.hdf5','*_rez.mat','*_jrc.mat','*.csv','*_spikesResorted.mat'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% do not include those files:
spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'DeepCut','Whisker','trial'}),...
    {spikeSortingFiles.name}));
sessionDir=cd;

dataFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*.bin','*raw.kwd','*RAW*Ch*.nex','*.ns*'},'UniformOutput', false);
dataFiles=vertcat(dataFiles{~cellfun('isempty',dataFiles)});
% keep those files
TTLFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_trialTTLs'}),...
    {dataFiles.name})); %not used here. Used for Phototag plots
dataFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_export';'_traces'}),...
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
videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,{'_VideoFrameTimes','vSync'}),...
    {videoFrameTimeFiles.name}));

% decide which file to use
% recName='vIRt22_1016_5100_50ms1Hz10mW';
spikeFileNum=1; dataFileNum=1; TTLFileNum=1; wTrackNumFile=1; %[1,2]; vFrameFileNum=1;
if sum(cellfun(@(flnm) contains(flnm,'vSync'),{videoFrameTimeFiles.name}))
    videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,'vSync'),...
    {videoFrameTimeFiles.name}));
else
    videoFrameTimeFiles=videoFrameTimeFiles(cellfun(@(flnm) contains(flnm,'_VideoFrameTimes'),...
    {videoFrameTimeFiles.name}));

end

%% Load spikes and recording traces
recDir=spikeSortingFiles(spikeFileNum).folder;
recName=spikeSortingFiles(spikeFileNum).name;
cd(recDir)
%also get some info about the recording if possible
% e.g., load rec_info from the spike file:
load(spikeSortingFiles(spikeFileNum).name,'rec_info')
if ~exist('rec_info','var')
    load([dataFiles(spikeFileNum).name(1:end-10) 'recInfo']);
else
    recInfo=rec_info; clear rec_info; 
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

% eval(['jrc ''load-bin'' ' bla ' int16'])
% trWav_raw = load_bin_(dataFileName, 'int16',[32 14 recDuration]);
% vcFile=fullfile(dataFileDir,dataFileName);
% fid=fopen(dataFileName, 'r');
% fclose(bla)
%         fid = fopen(vcFile, 'r');
% %         if header>0, fseek(fid, header, 'bof'); end
%         if isempty(dimm) % read all
%             S_file = dir(vcFile);
%             if numel(S_file)~=1, return; end % there must be one file
%             nData = floor((S_file(1).bytes - header) / bytesPerSample_(dataType));
%             dimm = [nData, 1]; %return column
%         end
% 
% fread_(foo,[32 14 187354],'int16');
% fclose(foo);
 
traces = memmapfile(fullfile(dataFileDir,dataFileName),'Format','int16');
allTraces=double(traces.Data); 
recDuration=int64(length(allTraces)/numElectrodes);
try
    allTraces=reshape(allTraces,[numElectrodes recDuration]);
%     allTraces=reshape(trWav_raw,[numElectrodes size(trWav_raw,2)*size(trWav_raw,3)]);
catch
    allTraces=reshape(allTraces',[recDuration numElectrodes]);
end
filterTraces=false;

spikes=LoadSpikeData(recName,traces);

%% Load TTLs (are they needed?)
try
    TTLDir=TTLFiles(TTLFileNum).folder;
    TTLFileName= TTLFiles(TTLFileNum).name;%[regexp(recName,'\S+?(?=_export)','match','once') '_TTLs.dat'];
    fid = fopen(fullfile(TTLDir,TTLFileName), 'r');
    TTLTimes = fread(fid,[2,Inf],'int32');
    fclose(fid);
    TTLs.times=TTLTimes(1,:);
    TTLs.samplingRate=1000;
catch
    TTLs=[];
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
    % videoFrameTimes=readVideoTTLData(dirListing);
end

%% Import whisker tracking data (aka "thetas")
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
                    ContinuityWhiskerID(whiskerTrackingData);
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
            whiskerTrackingData=WhiskerAngleSmoothFill(whiskerTrackingData(:,1),whiskerTrackingData(:,2));
            %     figure; hold on; plot(whiskerTrackingData)
            %adjust base angle if needed, e.g. 45degrees at full retraction:
%             whiskerTrackingData=whiskerTrackingData-min(whiskerTrackingData)+45;
        end
    elseif contains(whiskerTrackFileName,'.avi') %video file to extract whisker angle
        whiskerTrackingData=ExtractMultiWhiskerAngle_FFTonContours(fullfile(dirName,fileName));
        whiskerTrackingData=smoothdata(whiskerTrackingData,'rloess',20);
    else
        load([dirName fileName]);
    end
end

%% Recording start time (mostly for OE)
% Processor: Rhythm FPGA Id: 100 subProcessor: 0 start time: 27306000@30000Hz
if exist('recInfo','var') & isfield(recInfo,'recordingStartTime')
    startTime=double(recInfo.recordingStartTime); % 27306000; % !!!
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
if ~isempty(TTLs)
    TTLs.times=TTLs.times-double(startTime/(samplingRate/1000));
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

%% Find best units
spikes.unitID=double(spikes.unitID);
spikes.times=double(spikes.times);
% find most frequent units
[unitFreq,uniqueUnitIDs]=hist(spikes.unitID,unique(spikes.unitID));
[unitFreq,freqIdx]=sort(unitFreq','descend');
unitFreq=unitFreq./sum(unitFreq)*100; uniqueUnitIDs=uniqueUnitIDs(freqIdx);
bestUnitsIdx=find(unitFreq>0.1);
keepUnits=uniqueUnitIDs(bestUnitsIdx); keepUnits=sort(keepUnits(keepUnits~=0));
if isfield(spikes,'preferredElectrode')
    try
        titularChannels = unique(spikes.preferredElectrode(ismember(spikes.unitID,keepUnits)));
    catch
        titularChannels =find(~cellfun('isempty',spikes.preferredElectrode));
    end
end
% keepUnits=[1 2 3];
% titularChannels=[10 10 10];
keepTraces=titularChannels; %14; %[10 14 15];% keepTraces=1:16; %[10 14 15];

%% Keep selected recording trace and spike times,
recordingTraces=allTraces(keepTraces,:); %select the trace to keep
try
keepUnitsIdx=ismember(spikes.preferredElectrode,keepTraces);
unitID=spikes.unitID(keepUnitsIdx);
preferredElectrode=spikes.preferredElectrode(keepUnitsIdx);
try
    waveForms=spikes.waveforms(keepUnitsIdx,:);
catch
    waveForms=[];
end
spikeTimes=spikes.times(keepUnitsIdx);
catch
    unitID=spikes.unitID;
    spikeTimes=spikes.times;
    waveForms=spikes.waveforms;
    preferredElectrode=spikes.preferredElectrode;
end
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

% figure('Color','white');
% for traceNum=1:length(keepTraces)
% %     subplot(length(keepTraces),1,traceNum); 
%     figure('Color','white'); hold on;
%     plot(recordingTraces(traceNum,:));
%     correspondingUnits=unique(unitID(ismember(preferredElectrode,keepTraces(traceNum))));
%     for clusterNum=1:numel(correspondingUnits)
%         plot(spikeTimes(unitID==correspondingUnits(clusterNum)),...
%             ones(numel(spikeTimes(unitID==correspondingUnits(clusterNum))),1)*-300,'*');
%         %     plot(rs_spikeTimes,ones(size(rs_spikeTimes,1),1)*-300,'d')
%     end
%     midRec=round(size(recordingTraces,2)/2);
%     set(gca,'ylim',[-500 500],'xlim',[midRec-samplingRate midRec+samplingRate]);
% end

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
if numel(whiskerTrackingData)~=numel(vidTimes)
    whiskerTrackingData=whiskerTrackingData(1:numel(vidTimes)); %but check why that is
end
% periodBehavData=[whiskerTrackingData',vidTimes'];
% periodBehavData=[whiskerTrackingData(videoFrameTimes.TTLFrames(1):...
%                   size(videoFrameTimes.frameTime_ms,1)),... %Trace
%     videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1):end)-...
%     videoFrameTimes.frameTime_ms(videoFrameTimes.TTLFrames(1))+1]; % Time points
% figure; hold on
% plot(periodBehavData(:,2),periodBehavData(:,1))

%resample to 1ms precision
vidTimes_ms=vidTimes/samplingRate*1000;
% [periodBehavData(:,1),periodBehavData(:,2)] = resample(periodBehavData(:,1),periodBehavData(:,2),'pchip');
periodBehavData_ms=interp1(vidTimes_ms,whiskerTrackingData,...
    vidTimes_ms(1):vidTimes_ms(end));
% figure; plot(periodBehavData_ms);

%% Plot behavior data and find a period with whisking bouts
% no need to keep periods with no whisking
% figure; plot(periodBehavData); % select data point and export cursor info
% whiskingPeriod=1:cursor_info.Position(1); %in ms
peakWhisking_ms=diff(cummax(abs(diff(periodBehavData_ms))));
% peakWhiskingIdx=find(peakWhisking_ms==max(peakWhisking_ms));
% whiskingPeriod=peakWhiskingIdx-5000:peakWhiskingIdx+4999; %in ms

%% Filter periodic behavior traces into low-pass and bandpassed versions
LP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,0.3,'low'); %set-point
% figure; hold on
% plot(periodBehavData_ms); plot(LP_periodBehavData_ms,'LineWidth',2)

BP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,[0.3 20],'bandpass'); %whisking
% figure; hold on %plot(foo) % plot(periodBehavData_ms)
% plot(periodBehavData_ms-mean(periodBehavData_ms)); plot(BP_periodBehavData_ms,'LineWidth',1)

HP_periodBehavData_ms=FilterTrace(periodBehavData_ms,1000,0.3,'high')'; %whisking
% plot(HP_periodBehavData_ms,'LineWidth',1)

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
% Hilbert transform NEEDS ANGLE TO BE ZERO CENTERED !!!
% periodicSignal=smoothdata(BP_periodBehavData_ms,'robustfit');
% -mean(BP_periodBehavData_ms);
% baseSignal=FilterTrace(periodicSignal,1000,0.3,'low');
% baseSignal=LP_periodBehavData_ms-mean(LP_periodBehavData_ms);
% figure; hold on; plot(periodicSignal);
% periodicSignal(abs(periodicSignal)<2*std(abs(periodicSignal)))=0;plot(baseSignal);
HTBP_periodBehavData_ms=hilbert(BP_periodBehavData_ms);
% figure; plot(imag(HTBP_periodBehavData_ms));
% foo=abs(imag(HTBP_periodBehavData_ms))<mad((imag(HTBP_periodBehavData_ms)));
% bla=imag(HTBP_periodBehavData_ms); bla(foo)=0; plot(bla)
% plot(imag(HTBP_periodBehavData_ms));
whiskingPhase_ms=angle(HTBP_periodBehavData_ms);
% whiskingPhase_ms(foo)=0;
% plot(whiskingPhase_ms)
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






