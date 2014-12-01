%% Process and display data from TBSI recordings

%cd to folder

%% order files
DirListing=dir;
FileDates=cell2mat({DirListing(:).datenum});
[~,FDateIdx] = sort(FileDates,'descend');
DirListing = {DirListing(:).name};
DirListing=DirListing(FDateIdx);

%% get spike time stamps
prompt = 'Take spike timestamps from which channel ?';
ChannelSelect= input(prompt);
SpikeData = readNexFile(['SE-Spike-TS-Ch' num2str(ChannelSelect) '_.nex']);
SpikeTimes=round(SpikeData.neurons{1,1}.timestamps*1000); %converted to millisecondes

%% get Synchronizing TTL times
AnalogFile=DirListing(~cellfun('isempty',strfind(DirListing,'Analog')));
AnalogData = readNexFile(AnalogFile{:});
SyncTimes=AnalogData.contvars{1,1}.data;

SyncTimes(SyncTimes<3000)=0;
SyncTimes(SyncTimes>3000)=1;
SyncTimes=bwlabeln(SyncTimes);

if max(SyncTimes)
    for PulseNb=1:max(SyncTimes)
        if length(SyncTimes(SyncTimes==PulseNb))<=1000 %just artifact, not trigger
            SyncTimes(SyncTimes==PulseNb)=0;
        end
    end
end

SyncTimes=bwlabeln(SyncTimes);
NbPulses=max(SyncTimes);
SyncTimes=round(find(diff(SyncTimes)>1)./AnalogData.freq*1000); %converted to millisecondes

%% plot spikes aligned to sync times
RecDuration=ceil(size(AnalogData.contvars{1,1}.data,1)/AnalogData.contvars{1,1}.ADFrequency*1000);
SyncVect=zeros(1,RecDuration);SyncVect(SyncTimes)=1;
TrialInterval=diff(find(SyncVect,2));
SpikeVect=zeros(1,RecDuration);SpikeVect(SpikeTimes)=1;
% Rasters=nan(NbPulses,floor(TrialInterval/100)*100);
% get all samples (~trials)
Rasters=cellfun(@(x) SpikeVect(x-200:x+300),mat2cell(find(SyncVect)',ones(NbPulses,1)), 'UniformOutput',false);
% concatenate
Rasters=cat(1,Rasters{:});

cut_rasters = rasters(:,start:stop); % Isolate rasters of interest

[indy, indx] = ind2sub(size(cut_rasters),find(cut_rasters)); %find row and column coordinates of spikes

plot([indx';indx'],[indy';indy'+1],'color',cc(rastnum,:),'LineStyle','-'); % plot rasters
set(gca,'xlim',[1 length(start:stop)]);
axis(gca, 'off'); % axis tight sets the axis limits to the range of the data.



