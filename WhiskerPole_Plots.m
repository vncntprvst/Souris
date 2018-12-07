% WhiskerPole_Plots

%% First, get data
[spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
    SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
    templates,templateIdx,stimType,stimTimes_ms] = WhiskerPole_GatherData;

%% plot units with stim
figure('Position',[1050 120 750 790]); hold on 
colormapSeed=lines; cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
patchColor=winter;
imagesc(zscore(SDFs_ms,[]));

for ROINum=1:length(stimTimes_ms)
    ROIstimTimes=stimTimes_ms{ROINum};
    for stimNum=1:length(ROIstimTimes)
        patch([ROIstimTimes(stimNum), ROIstimTimes(stimNum),...
            ROIstimTimes(stimNum)+50, ROIstimTimes(stimNum)+50], ...
            [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
            [0 0 0 0],patchColor(ROINum,:),'EdgeColor','none','FaceAlpha',0.7);
    end
end

%% align data to stimulus
preAlignWindow=200;
postAlignWindow=999;
stimTimeShift=400;
PSTH=nan(size(SDFs_ms,1),1200);
for unitNum=1:size(SDFs_ms,1)
    stimSDF=nan(numel(ROIstimTimes),preAlignWindow+postAlignWindow+1);
    for stimNum=1:numel(ROIstimTimes)
        try
            stimSDF(stimNum,:)=SDFs_ms(unitNum,...
                int32(ROIstimTimes(stimNum)+stimTimeShift-preAlignWindow:...
                ROIstimTimes(stimNum)+stimTimeShift+postAlignWindow));
            %smoothed:
            %             PSTH(stimNum,:)=convSpikeTime(...
            %                 ROIstimTimes(stimNum)-preAlignWindow:...
            %                 ROIstimTimes(stimNum)+postAlignWindow);
        catch
            continue
        end
    end
    PSTH(unitNum,:)=nanmean(stimSDF);
end

%% plot psth showing all units
keepUnits(9)=[];
figure('Position',[1050 120 750 790]); colormap bone;
% meanUnits=cellfun(@(x) conv_raster(x),{PSTH.startTrial},'UniformOutput',false);
% meanUnits=cell2mat(meanUnits');
imagesc(zscore(PSTH(keepUnits,:),[])); %
% imagesc(meanUnits);
xlabel('Time (ms)');
ylabel('Units','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
patch([repmat(preAlignWindow-2,1,2) repmat(preAlignWindow+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

set(gca,'XTick',0:100:sum(preAlignWindow+postAlignWindow));
set(gca,'XTickLabel',-preAlignWindow:100:postAlignWindow);

title('whisker touch responses');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';
