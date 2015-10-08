%% Process and display data from Open-Ephys recordings  

%% Get file path
% e.g., C:\Data\OpenEphys\PrV25_12_train2_2015-09-10_20-26-30
dname = uigetdir('C:\Data\OpenEphys\','Open Ephys Reordings');    
cd(dname);

%% Get Spike data
%O_E binary format - see also 'load_open_ephys_data.m'
% fid = fopen('102_CH1.continuous'); 
% fid = fopen('all_channels.events');
% hdr = fread(fid, 1024, 'char*1');
% timestamp = fread(fid, 1, 'int64',0,'l');
% N = fread(fid, 1, 'uint16',0,'l');
% recordingNumber = fread(fid, 1, 'uint16', 0, 'l');
% samples = fread(fid, N, 'int16',0,'b');
% recordmarker = fread(fid, 10, 'char*1');
% fclose(fid);

%kwik data
% KWIK contains metadata about the experiment, as well as event data
% h5disp('experiment1.kwik','/recordings/');
SamplingRate=h5readatt('experiment1.kwik','/recordings/0','sample_rate');

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
    figure;
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
        subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),ChanNum+1)
        plot(sum(ChanData(ChanNum+1).Waveforms,3));
        %assuming 48 samples here
        set(gca,'XTick',0:16:48);
        set(gca,'XTickLabel',round(-(1/SamplingRate)*16*10000)/10:0.5:round((1/SamplingRate)*32*10000)/10);
        axis(gca,'tight'); box off;
        set(gca,'Color','white','TickDir','out');
        title(['Chan ' num2str(ChanNum+1)]);
        
        MinMax(1)=min(MinMax(1), min(get(gca,'ylim')));
        MinMax(2)=max(MinMax(2), max(get(gca,'ylim')));
    end
    for chan=1:length(KeepChans)
        subplot(ceil(sqrt(length(KeepChans))),ceil(length(KeepChans)/ceil(sqrt(length(KeepChans)))),chan)
        set(gca,'ylim',MinMax)
    end
end

%% let user select Channels to plot
% ....

%% get Trial structure
% h5disp('experiment1.kwik','/event_types/TTL')
% TTLinfo=h5info('experiment1.kwik','/event_types/TTL');

TTL_ID = h5read('experiment1.kwik','/event_types/TTL/events/user_data/eventID');
TTL_times = h5read('experiment1.kwik','/event_types/TTL/events/time_samples');
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

%% gather data
Spkt=ChanData(5).SpikeTimes;
Rasters{1}=zeros(size(Trials.start,1),2000); %begining of trial
Rasters{2}=zeros(size(Trials.end,1),2000); %end of trial
for trialnb=1:size(Trials.start,1)
%Collect spikes from 1st epoch (begining of trial)
RastSWin=Trials.start(trialnb)-uint64(SamplingRate); % 1 sec before
RastEWin=Trials.start(trialnb)+uint64(SamplingRate); % 1 sec afer
SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
SpikeTimes(SpikeTimes==0)=1; %no 0 indices
SpikeTimes=unique(SpikeTimes);
Rasters{1}(trialnb,SpikeTimes)=1;
%Collect spikes from 2nd epoch (end of trial)
RastSWin=Trials.end(trialnb)-uint64(SamplingRate); % 1 sec before
RastEWin=Trials.end(trialnb)+uint64(SamplingRate); % 1 sec afer
SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
SpikeTimes(SpikeTimes==0)=1; %no 0 indices
SpikeTimes=unique(SpikeTimes);
Rasters{2}(trialnb,SpikeTimes)=1;
end

%% plot raster

%% plot sdf
start=1;
stop=size(Rasters{1},2);
conv_sigma=20;
alignmtt=1000;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters{1},conv_sigma,start,stop);
[sdf{2}, ~, rastsem{2}]=conv_raster(Rasters{2},conv_sigma,start,stop);
figure;
colormap default;
cmap = colormap(gcf);
hold on;

%plot sem
patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
patch([1:length(sdf{2}),fliplr(1:length(sdf{2}))],[sdf{2}-rastsem{2},fliplr(sdf{2}+rastsem{2})],cmap(22,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdf
plot(sdf{1},'Color',cmap(1,:),'LineWidth',1.8);
plot(sdf{2},'Color',cmap(22,:),'LineWidth',1.8);

set(gca,'XTick',0:100:(stop-start-6*conv_sigma));
set(gca,'XTickLabel',-(alignmtt-(start+3*conv_sigma)):100:stop-(alignmtt+3*conv_sigma));
axis(gca,'tight'); box off;
set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
hylabel=ylabel(gca,'Firing rate (spikes/s)','FontName','Cambria','FontSize',10);

% draw alignment bar
currylim=get(gca,'YLim');
patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

%legend
legend('Front panel exploration','Reward port');