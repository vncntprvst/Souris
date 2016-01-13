% compare offline spike sorting results with raw data and online spike sorting
function validate_OfflineSpikeSorting

%get most recently changed data folder
dataDir='C:\Data\';
dataDirListing=dir(dataDir);
[~,fDateIdx]=sort([dataDirListing.datenum],'descend');
recentDataFolder=[dataDir dataDirListing(fDateIdx(1)).name '\'];
%get most recent data folder in that folder
dataDirListing=dir(recentDataFolder);
[~,fDateIdx]=sort([dataDirListing.datenum],'descend');
recentDataFolder=[recentDataFolder dataDirListing(fDateIdx(1)).name '\'];

%% Get file path
[fname,dname] = uigetfile({'*.continuous;*.kwik;*.kwd;*.kwx;*.nex;*.ns*','All Data Formats';...
    '*.*','All Files' },'Most recent data',recentDataFolder);
cd(dname);
% expname=regexp(strrep(dname,'-','_'),'\\\w+','match');
disp(['loading ' dname fname]);
fname=regexp(fname,'^\w+\_|^\w+\.','match');fname=fname{:}(1:end-1);

%% find file format from directory listing
dirlisting = dir(dname);
% filedates=cell2mat({dirlisting(:).datenum}); % Order by date
% [~,fdateidx] = sort(filedates,'descend');
dirlisting = {dirlisting(:).name};
% dirlisting=dirlisting(fdateidx);
dirlisting=dirlisting(cellfun('isempty',cellfun(@(x) strfind('.',x(end)),dirlisting,'UniformOutput',false)));
fileformats={'continuous','kwe','kwik','nex'};
whichformat=cellfun(@(x) find(~cellfun('isempty',strfind(fileformats,x(end-2:end)))),dirlisting,'UniformOutput',false);
whichformat=fileformats(unique([whichformat{~cellfun('isempty',whichformat)}]));

%% raw data info
rawInfo=h5info([fname '_100.raw.kwd']);%'/recordings/0/data'
rawInfo=h5info([fname '_100.raw.kwd'],rawInfo.Groups.Name);

rec.dur=rawInfo.Groups.Datasets.Dataspace.Size;
rec.samplingRate=h5readatt([fname '_100.raw.kwd'],rawInfo.Groups.Name,'sample_rate');
rec.bitDepth=h5readatt([fname '_100.raw.kwd'],rawInfo.Groups.Name,'bit_depth');
rec.numRecChan=rawInfo.Groups.Datasets.Dataspace.Size-3;  %number of raw data channels.
% Last 3 are headstage's AUX channels (e.g accelerometer)
rawData=h5read([fname '_100.raw.kwd'],'/recordings/0/data',[1 1],[rec.numRecChan(1) Inf]);

%% online spike sorting
% KWX contains the spike data
% h5disp('experiment1.kwx')
ChanInfo=h5info([fname '.kwx']); %whichformat{:}

%Keep only channels with > 1Hz firing rate
GoodChans=cell(size(ChanInfo.Groups.Groups,1),1);
for chan=1:size(ChanInfo.Groups.Groups,1)
    if ChanInfo.Groups.Groups(chan).Datasets(2).Dataspace.Size/(rec.dur(2)/rec.samplingRate)>1
        GoodChans(chan)=regexp(ChanInfo.Groups.Groups(chan).Name,'\d+$','match');
    end
end

% get data from all channels with appropriate units
KeepChans=find(~cellfun('isempty',GoodChans));
ChanData=struct('Units',[],'SpikeTimes',[],'Waveforms',[]);
if ~isempty(KeepChans)
    for chan=1:length(KeepChans)
        %     figure('Position',[35 122 884 845]);
        %     MinMax=[0 0];
        ChanNum=str2double(GoodChans{KeepChans(chan)});
        ChanData(ChanNum+1).Units=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/recordings']);
        try
            ChanData(ChanNum+1).onssSpikeTimes=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/time_samples']);
            ChanData(ChanNum+1).onssWaveforms=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/waveforms_filtered']);
        catch
            ChanData(ChanNum+1).onssSpikeTimes=[];
            ChanData(ChanNum+1).onssWaveforms=[];
        end
%         % plot waveform
%         subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),ChanNum+1)
%         plot(sum(ChanData(ChanNum+1).Waveforms,3));
%         %assuming 48 samples here
%         set(gca,'XTick',0:16:48);
%         set(gca,'XTickLabel',round(-(1/rec.samplingRate)*16*10000)/10:0.5:round((1/rec.samplingRate)*32*10000)/10);
%         axis(gca,'tight'); box off;
%         set(gca,'Color','white','TickDir','out');
%         title(['Chan ' num2str(ChanNum+1)]);
%         
%         MinMax(1)=min(MinMax(1), min(get(gca,'ylim')));
%         MinMax(2)=max(MinMax(2), max(get(gca,'ylim')));
    end
%     for chan=1:length(KeepChans)
%         subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),chan)
%         set(gca,'ylim',MinMax)
%     end
end



%% raw data
%load only recording channels
rawData=h5read([fname '_100.raw.kwd'],'/recordings/0/data',[1 1],[rec.numRecChan(1) Inf]);

%% offline spike sorting
cd('C:\Data\export\PrV63_20151209_130657_OEphVnopp')
offssInfo=h5info('PrV63_20151209_130657_OEphVnopp.result.hdf5')

%%%%% need to arrange by channels, not just cluster id 
for id = 1:size(offssInfo.Groups(1).Datasets, 1)
    offssData = h5read('PrV63_20151209_130657_OEphVnopp.result.hdf5', ['/spiketimes/temp_' int2str(id - 1)]);
    if size(data, 1) == 0 %|| (same(offssData, zeros(1)) > 0)
        ChanData.offssSpikeTimes{id} = zeros(0, 1);
    else
        ChanData.offssSpikeTimes{id} = double(data)/(rec.samplingRate/1000);
    end
end


        %     figure('Position',[35 122 884 845]);
        %     MinMax=[0 0];
        
        %         % plot waveform
%         subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),ChanNum+1)
%         plot(sum(ChanData(ChanNum+1).Waveforms,3));
%         %assuming 48 samples here
%         set(gca,'XTick',0:16:48);
%         set(gca,'XTickLabel',round(-(1/rec.samplingRate)*16*10000)/10:0.5:round((1/rec.samplingRate)*32*10000)/10);
%         axis(gca,'tight'); box off;
%         set(gca,'Color','white','TickDir','out');
%         title(['Chan ' num2str(ChanNum+1)]);
%         
%         MinMax(1)=min(MinMax(1), min(get(gca,'ylim')));
%         MinMax(2)=max(MinMax(2), max(get(gca,'ylim')));

%     for chan=1:length(KeepChans)
%         subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),chan)
%         set(gca,'ylim',MinMax)
%     end

end