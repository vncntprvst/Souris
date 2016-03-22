%% Align and display rasters and sdf

%% Get file path
[fileName,dirName] = uigetfile({'*.mat; *.hdf5','Processed data';'*.dat','Flat data';...
    '*.*','All Files' },'Exported data','C:\Data\export');
cd(dirName);

% load file data
if strfind(fileName,'.mat')
    load(fileName);
    Behavior=processBehaviorData;
elseif strfind(fileName,'.hdf5')
    Spikes.Offline_SpkSort.data{16,1}=h5read(fileName,'/times_15');
end

% plot all spikes and trials
figure; hold on
plot(Spikes.Offline_Threshold.data{12, 1}  )

%Trials times are already converted in ms scale (
% -> should keep original timing
% -> change Spikes.clockTimes to Trials.clockTime

% Unless I screwed up and used StimPulse_kbControl, instead of
% Laser_Pulse_Control, stim params were 20ms High, 80ms Low.
% But Trials.end(1)-Trials.start(1) = 200ms. And looks like there may be
% upwards of 18ms delay on Digital IN detection. 
% -> check TTL_times in getOE_Trials
% -> check TTL rise and fall time. Plug TTL input to Analog IN.


for TTLNum=1:size(Trials.start,1)
patch([Trials.start(TTLNum)-Trials.startClockTime:Trials.end(TTLNum)-Trials.startClockTime,...
    fliplr(Trials.start(TTLNum)-Trials.startClockTime:Trials.end(TTLNum)-Trials.startClockTime)],...
    [zeros(1,Trials.end(TTLNum)-Trials.start(TTLNum)+1),...
    ones(1,Trials.end(TTLNum)-Trials.start(TTLNum)+1)*2],...
    [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
end

% switch whichformat
%     case
%     case
% end

%Select channels
% KeepChans=Spikes.channel;
KeepChans=15;

%Keep only common trials
[ephysCommonTrials, behaviorCommonTrials]=match_trials(Trials,Behavior);

%% gather data from selected channels
Rasters.channels=cell(length(KeepChans),2);
Rasters.epochnames={'BeginTrial','EndTrial'};
for chan=1:length(KeepChans)
    downSamplingRatio=uint64(Spikes.samplingRate(chan,1)/Spikes.samplingRate(chan,2));
    [Rasters.channels{chan,1},Rasters.channels{chan,2}]=deal(zeros(size(Trials.start,1),1501));
    Spkt=[zeros(1,round(Trials.startClockTime/downSamplingRatio)) Spikes.data{KeepChans(chan),2}(1,:)];
    if Trials.end(end)>size(Spkt,2)
        continue
    end
    for trialnb=1:size(Trials.start,1)
        
        %Collect spikes from 1st epoch (begining of trial)
        RastSWin=Trials.start(trialnb)-uint64(Spikes.samplingRate(chan,2)); % 1 sec before
        RastEWin=Trials.start(trialnb)+round(uint64(Spikes.samplingRate(chan,2))/2); % 1/2 sec afer
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,1}(trialnb,:)=SpikeTimes;
        %Collect spikes from 2nd epoch (end of trial)
        RastSWin=Trials.end(trialnb)-uint64(Spikes.samplingRate(chan,2)); % 1 sec before
        RastEWin=Trials.end(trialnb)+round(uint64(Spikes.samplingRate(chan,2))/2); % 1/2 sec afer
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,2}(trialnb,:)=SpikeTimes;

        %         %Collect spikes from 1st epoch (begining of trial)
        %         RastSWin=Trials.start(trialnb)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.start(trialnb)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,1}(trialnb,SpikeTimes)=1;
        %         %Collect spikes from 2nd epoch (end of trial)
        %         RastSWin=Trials.end(trialnb)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.end(trialnb)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,2}(trialnb,SpikeTimes)=1;
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