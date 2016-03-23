%% Align and display rasters and sdf

%% Get file path
[fileName,dirName] = uigetfile({'*.mat; *.hdf5','Processed data';'*.dat','Flat data';...
    '*.*','All Files' },'Exported data','C:\Data\export');
cd(dirName);

% load file data
if strfind(fileName,'.mat')
    load(fileName);
    try
        Behavior=processBehaviorData;
    catch
        Behavior=[];
    end
elseif strfind(fileName,'.hdf5')
    Spikes.Offline_SpkSort.data{16,1}=h5read(fileName,'/times_15');
end

%
Trials.start=Trials.start-Trials.startClockTime;
Trials.end=Trials.end-Trials.startClockTime;

% plot all spikes and trials
figure; hold on
plot(Spikes.Offline_Threshold.data{16, 1},'k')
plot(Spikes.Offline_SpkSort.data{16, 1},ones(1,size(Spikes.Offline_SpkSort.data{16, 1},1))*0.5,'sr')

for TTLNum=1:size(Trials.start,1)
patch([Trials.start(TTLNum):Trials.end(TTLNum),...
    fliplr(Trials.start(TTLNum):Trials.end(TTLNum))],...
    [zeros(1,Trials.end(TTLNum)-Trials.start(TTLNum)+1),...
    ones(1,Trials.end(TTLNum)-Trials.start(TTLNum)+1)],...
    [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
end

%Select channels
KeepChans=Spikes.Offline_Threshold.channel;
% KeepChans=16;

%Keep only verified trials
[ephysCommonTrials, behaviorCommonTrials]=match_trials(Trials,Behavior);

%% gather data from selected channels
Rasters.channels=cell(length(KeepChans),2);
Rasters.epochnames={'BeginTrial','EndTrial'};
preAlignWindow(1)=0.08; %80ms
preAlignWindow(2)=round(uint64(Spikes.Offline_Threshold.samplingRate(chan,2))/(1/preAlignWindow(1)));
postAlignWindow(1)=0.08; %80ms
postAlignWindow(2)=round(uint64(Spikes.Offline_Threshold.samplingRate(chan,2))/(1/postAlignWindow(1)));
for chan=1:length(KeepChans)
%     downSamplingRatio=uint64(Spikes.Offline_Threshold.samplingRate(chan,1)/Spikes.Offline_Threshold.samplingRate(chan,2));
    [Rasters.channels{chan,1},Rasters.channels{chan,2}]=deal(zeros(size(Trials.start,1),preAlignWindow(2)+postAlignWindow(2)+1));
    Spkt=Spikes.Offline_Threshold.data{KeepChans(chan),2}(1,:);
    if Trials.end(end)>size(Spkt,2)
        continue
    end
    for trialnb=1:size(Trials.start,1)
        
        %Collect spikes from 1st epoch (begining of trial)
        RastSWin=Trials.start(trialnb,2)-preAlignWindow(2); 
        RastEWin=Trials.start(trialnb,2)+postAlignWindow(2); 
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,1}(trialnb,:)=SpikeTimes;
        %Collect spikes from 2nd epoch (end of trial)
        RastSWin=Trials.end(trialnb,2)-preAlignWindow(2); % 1 sec before
        RastEWin=Trials.end(trialnb,2)+postAlignWindow(2); % 1/2 sec afer
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,2}(trialnb,:)=SpikeTimes;

        %         %Collect spikes from 1st epoch (begining of trial)
        %         RastSWin=Trials.start(trialnb,2)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.start(trialnb,2)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,1}(trialnb,SpikeTimes)=1;
        %         %Collect spikes from 2nd epoch (end of trial)
        %         RastSWin=Trials.end(trialnb,2)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.end(trialnb,2)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,2}(trialnb,SpikeTimes)=1;
    end
end
%% plot raster showing all channels
figure('Position',[1050 120 750 790]);
subplot(1,2,1)
colormap bone;
MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,1),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
% subplot(1,2,1)
imagesc(zscore(MeanChan,[],2)); %
% imagesc(MeanChan);
xlabel('Time (ms)');
ylabel('Channels','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim(2)/20)*10;
set(gca,'XTick',[midl-preAlignWindow(2)/2 midl midl+postAlignWindow(2)/2]);
set(gca,'XTickLabel',[-preAlignWindow(2)/2 0 postAlignWindow(2)/2]);
%opto stim patch
patch([repmat(midl,1,2) repmat(midl+Trials.end(1,2)-Trials.start(1,2),1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.1);

% patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
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
title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

%% plot sdf
BestChan=find(mean(MeanChan,2)==max(mean(MeanChan,2)));
start=1;
stop=size(Rasters.channels{BestChan,1},2);
conv_sigma=1;
alignmtt=preAlignWindow(2);
xTickSteps=round(preAlignWindow(2)/50)*10;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters.channels{BestChan,1},conv_sigma,start,stop);
[sdf{2}, ~, rastsem{2}]=conv_raster(Rasters.channels{BestChan,2},conv_sigma,start,stop);
% figure('Position',[1469 542 417 417]);
subplot(1,2,2)
colormap default;
cmap = colormap(gcf);
hold on;

%plot sem
startAlignPloth=gca; box off; %subplot(1,2,1);hold on; box off;
patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],...
    [0.16 0.38 0.27],'EdgeColor','none','FaceAlpha',0.2);
% endAlignPloth=subplot(1,2,2);hold on; box off;
% patch([1:length(sdf{2}),fliplr(1:length(sdf{2}))],[sdf{2}-rastsem{2},fliplr(sdf{2}+rastsem{2})],cmap(22,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdfs
FRploth=plot(startAlignPloth,sdf{1},'Color',[0.16 0.38 0.27],'LineWidth',1.8);

% set(startAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
set(startAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
set(startAlignPloth,'XTickLabel',-(alignmtt-xTickSteps):xTickSteps:stop-(alignmtt+xTickSteps));
axis(startAlignPloth,'tight');
set(startAlignPloth,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
hxlabel=xlabel(startAlignPloth,'Time (ms)','FontName','Cambria','FontSize',12);
hylabel=ylabel(startAlignPloth,'Firing rate (spikes/s)','FontName','Cambria','FontSize',12);

% plot(endAlignPloth,sdf{2},'Color',cmap(22,:),'LineWidth',1.8);
% 
% set(endAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
% set(endAlignPloth,'XTickLabel',-(alignmtt-xTickSteps):xTickSteps:stop-(alignmtt+xTickSteps));
% axis(endAlignPloth,'tight'); 
% set(endAlignPloth,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
% hxlabel=xlabel(endAlignPloth,'Time (ms)','FontName','Cambria','FontSize',12);
% hylabel=ylabel(endAlignPloth,'Firing rate (spikes/s)','FontName','Cambria','FontSize',12);

% draw alignment bar
currylim=get(startAlignPloth,'YLim');
axes(startAlignPloth)
% opto stim bar
OptoStimh=patch([repmat(alignmtt-(start+3*conv_sigma),1,2) repmat(alignmtt-(start+3*conv_sigma)+Trials.end(1,2)-Trials.start(1,2),1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
% "regular" alignement bar
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
% axes(endAlignPloth)
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

%legend
legend([FRploth,OptoStimh],{'Average firing rate','Optical stimulation'});
legend('boxoff')
% text(xTickSteps,currylim(2)-20,['Channel ' num2str(BestChan)],'FontName','Cambria');
title(['Channel ' num2str(BestChan)],'FontName','Cambria');