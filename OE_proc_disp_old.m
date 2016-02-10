%% Process and display data from Open-Ephys recordings

%% Get file path
% e.g., C:\Data\OpenEphys\PrV25_12_train2_2015-09-10_20-26-30
dname = uigetdir('C:\Data\OpenEphys\','Open Ephys Recordings');
cd(dname);

% find file format from directory listing
dirlisting = dir(dname);
% filedates=cell2mat({dirlisting(:).datenum}); % Order by date
% [~,fdateidx] = sort(filedates,'descend');
dirlisting = {dirlisting(:).name};
% dirlisting=dirlisting(fdateidx);
dirlisting=dirlisting(cellfun('isempty',cellfun(@(x) strfind('.',x(end)),dirlisting,'UniformOutput',false)));
fileformats={'continuous','kwe','kwik','nex'};
whichformat=cellfun(@(x) find(~cellfun('isempty',strfind(fileformats,x(end-2:end)))),dirlisting,'UniformOutput',false);
whichformat=fileformats(unique([whichformat{~cellfun('isempty',whichformat)}]));

switch whichformat{:}
    case 'continuous'
        % Get Spike data
        O_E binary format - see also 'load_open_ephys_data.m'
        fid = fopen('102_CH1.continuous');
        fid = fopen('all_channels.events');
        hdr = fread(fid, 1024, 'char*1');
        timestamp = fread(fid, 1, 'int64',0,'l');
        N = fread(fid, 1, 'uint16',0,'l');
        recordingNumber = fread(fid, 1, 'uint16', 0, 'l');
        samples = fread(fid, N, 'int16',0,'b');
        recordmarker = fread(fid, 10, 'char*1');
        fclose(fid);
    case {'kwe','kwik'}
        %kwik data
        % KWIK contains metadata about the experiment, as well as event data
        % h5disp('experiment1.kwik','/recordings/');
        SamplingRate=h5readatt(['experiment1.' whichformat{:}],'/recordings/0','sample_rate');
        
        % KWX contains the spike data
        % h5disp('experiment1.kwx')
        ChanInfo=h5info('experiment1.kwx');
        
        %KWD contains the continuous data for a given processor in the /recordings/#/data dataset, where # is the recording number, starting at 0
        RawInfo=h5info('experiment1_100.raw.kwd','/recordings/0/data');
        RecDur=RawInfo.Dataspace.Size;
        
        %Keep only channels with > 1Hz firing rate
        GoodChans=cell(size(ChanInfo.Groups.Groups,1),1);
        for chan=1:size(ChanInfo.Groups.Groups,1)
            if ChanInfo.Groups.Groups(chan).Datasets(2).Dataspace.Size/(RecDur(2)/SamplingRate)>1
                GoodChans(chan)=regexp(ChanInfo.Groups.Groups(chan).Name,'\d+$','match');
            end
        end
        
        % get data from all channels with appropriate units
        KeepChans=find(~cellfun('isempty',GoodChans));
        ChanData=struct('Units',[],'SpikeTimes',[],'Waveforms',[]);
        if ~isempty(KeepChans)
%             figure('Position',[35 122 884 845]);
            MinMax=[0 0];
            for chan=1:length(KeepChans)
                ChanNum=str2double(GoodChans{KeepChans(chan)});
                ChanData(ChanNum+1).Units=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/recordings']);
                try
                    ChanData(ChanNum+1).SpikeTimes=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/time_samples']);
                    ChanData(ChanNum+1).Waveforms=h5read('experiment1.kwx',['/channel_groups/' GoodChans{KeepChans(chan)} '/waveforms_filtered']);
                catch
                    ChanData(ChanNum+1).SpikeTimes=[];
                    ChanData(ChanNum+1).Waveforms=[];
                end
                % plot waveform
%                 subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),ChanNum+1)
%                 plot(sum(ChanData(ChanNum+1).Waveforms,3));
%                 %assuming 48 samples here
%                 set(gca,'XTick',0:16:48);
%                 set(gca,'XTickLabel',round(-(1/SamplingRate)*16*10000)/10:0.5:round((1/SamplingRate)*32*10000)/10);
%                 axis(gca,'tight'); box off;
%                 set(gca,'Color','white','TickDir','out');
%                 title(['Chan ' num2str(ChanNum+1)]);
                
%                 MinMax(1)=min(MinMax(1), min(get(gca,'ylim')));
%                 MinMax(2)=max(MinMax(2), max(get(gca,'ylim')));
            end
%             for chan=1:length(KeepChans)
%                 subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),chan)
%                 set(gca,'ylim',MinMax)
%             end
        end
        
        %% let user select Channels to plot?
        % default: all
        KeepChans=1:size(ChanData,2);
        
        %% get Trial structure
        % h5disp('experiment1.kwik','/event_types/TTL')
        % TTLinfo=h5info('experiment1.kwik','/event_types/TTL');
        
        TTL_ID = h5read('experiment1.kwe','/event_types/TTL/events/user_data/eventID');
        TTL_times = h5read('experiment1.kwe','/event_types/TTL/events/time_samples');
        % keep absolute time of TTL onset
        TTL_times=TTL_times(diff([0;TTL_ID])>0);
        % TTL sequence
        TTL_seq=diff(TTL_times);
        %we send 10ms TTLs, with 10ms interval. Two pulses for begining of trial
        %(e.g.,head through front panel). With sampling rate of 30kHz, that
        % interval should be 601 samples (20ms*30+1). Or 602 accounting for jitter.
        % TTL_seq should thus read as:
        %   601
        %   end of trial time
        %   inter-trial interval
        if TTL_seq(1)>=610 %missed first trial initiation, discard times
            TTL_seq(1)=999;
        end
        if TTL_seq(end)<=610 %unfinished last trial
            TTL_seq(end)=999;
        end
        
        Trials.start=TTL_times([TTL_seq<=610;false]);%Trials.start=Trials.start./uint64(SamplingRate/1000)
        Trials.end=TTL_times(find([TTL_seq<=610;false])+2);
        Trials.interval=TTL_times(find([TTL_seq<=610;false])+3)-TTL_times(find([TTL_seq<=610;false])+2);
        
    case 'nex'
        % do something
    otherwise
        warning('Unexpected format type. No data loaded')
end
%% gather data from selected channels
Rasters.channels=cell(length(KeepChans),2);
Rasters.epochnames={'BeginTrial','EndTrial'};
for chan=1:length(KeepChans)
    [Rasters.channels{chan,1},Rasters.channels{chan,2}]=deal(zeros(size(Trials.start,1),2000));
    Spkt=ChanData(KeepChans(chan)).SpikeTimes;
    for trialnb=1:size(Trials.start,1)
        %Collect spikes from 1st epoch (begining of trial)
        RastSWin=Trials.start(trialnb)-uint64(SamplingRate); % 1 sec before
        RastEWin=Trials.start(trialnb)+uint64(SamplingRate); % 1 sec afer
        SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        SpikeTimes=unique(SpikeTimes);
        Rasters.channels{chan,1}(trialnb,SpikeTimes)=1;
        %Collect spikes from 2nd epoch (end of trial)
        RastSWin=Trials.end(trialnb)-uint64(SamplingRate); % 1 sec before
        RastEWin=Trials.end(trialnb)+uint64(SamplingRate); % 1 sec afer
        SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        SpikeTimes=unique(SpikeTimes);
        Rasters.channels{chan,2}(trialnb,SpikeTimes)=1;
    end
end
%% plot raster showing all channels
figure('Position',[1050 120 750 790]);

MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,1),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
subplot(1,2,1)
imagesc(zscore(MeanChan,[],2)); %
% imagesc(MeanChan);
xlabel('Time');
ylabel('Channel','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim(2)/2);
set(gca,'XTick',[midl-500 midl midl+500]);
set(gca,'XTickLabel',[-500 0 500]);
patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response, begining of trial');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,2),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
subplot(1,2,2)
imagesc(zscore(MeanChan,[],2));
% imagesc(MeanChan);
xlabel('Time');
ylabel('Channel','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim(2)/2);
set(gca,'XTick',[midl-500 midl midl+500]);
set(gca,'XTickLabel',[-500 0 500]);
patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response, end of trial');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

%% plot sdf
BestChan=find(mean(MeanChan,2)==max(mean(MeanChan,2)));
start=1;
stop=size(Rasters.channels{BestChan,1},2);
conv_sigma=20;
alignmtt=1000;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters.channels{BestChan,1},conv_sigma,start,stop);
[sdf{2}, ~, rastsem{2}]=conv_raster(Rasters.channels{BestChan,2},conv_sigma,start,stop);
figure('Position',[1469 542 417 417]);
colormap default;
cmap = colormap(gcf);
hold on;

%plot sem
patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
patch([1:length(sdf{2}),fliplr(1:length(sdf{2}))],[sdf{2}-rastsem{2},fliplr(sdf{2}+rastsem{2})],cmap(22,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdf
plot(sdf{1},'Color',cmap(1,:),'LineWidth',1.8);
plot(sdf{2},'Color',cmap(22,:),'LineWidth',1.8);

set(gca,'XTick',200-(start+3*conv_sigma):200:(stop-start-6*conv_sigma));
set(gca,'XTickLabel',-(alignmtt-200):200:stop-(alignmtt+200));
axis(gca,'tight'); box off;
set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',12);
hylabel=ylabel(gca,'Firing rate (spikes/s)','FontName','Cambria','FontSize',12);

% draw alignment bar
currylim=get(gca,'YLim');
patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

%legend
legend('Front panel exploration','Reward port');
legend('boxoff')
text(200,currylim(2)-20,['Channel ' num2str(BestChan)],'FontName','Cambria');