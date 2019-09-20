function NBC_Plots_SpikesWhisking(whiskingAngle,whiskingPhase,whiskingEpochsIdx,spikeRasters,whiskerTraces_breathFreq_ms)

% medWhiskAngleVal=median(periodBehavData_ms(1,:));
colormapSeed=lines;
cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
plotShift=6*std(whiskingAngle(1,:));
figure; hold on
% for whiskNum=1:3 %size(whiskingAngle,1)-1
%     plot(whiskingAngle(whiskNum,:),'Color',[whiskNum*0.1 whiskNum*0.1 whiskNum*0.1 0.5])
% end
whiskingEpochs=nan(1,numel(whiskingEpochsIdx));whiskingEpochs(whiskingEpochsIdx)=1;
plot(whiskingEpochs+plotShift,'linewidth',2,'Color','b')
% plot(whiskerTraces_breathFreq_ms(1,:)-plotShift)
plot(whiskingPhase(1,:));
plot(whiskingAngle(1,:),'Color','k','linewidth',1.8)
for spkNum=1:size(spikeRasters{1},1)
%% bars (pretty heavy on memory)
%     %find row and column coordinates of spikes
%     [indy, indx] = ind2sub(size(spikeRasters{1}(spkNum,:)),find(spikeRasters{1}(spkNum,:)));
%     % plot rasters
%     plot([indx;indx],[indy;indy+1]+min(whiskingAngle(1,:))-spkNum,'color',cmap(spkNum,:),'LineStyle','-');
%% diamonds
    plot(find(spikeRasters{1}(spkNum,:)),...
        ones(1,numel(find(spikeRasters{1}(spkNum,:))))-...
        plotShift-spkNum*1.5,'LineStyle','none',...
        'Marker','d','MarkerEdgeColor','none',...
        'MarkerFaceColor',cmap(spkNum,:),'MarkerSize',4);
end

for spkNum=1:size(spikeRasters{2},1)
    plot(find(spikeRasters{2}(spkNum,:)),...
        ones(1,numel(find(spikeRasters{2}(spkNum,:))))-...
        plotShift-spkNum*1.5-size(spikeRasters{1},1)*1.5,'LineStyle','none',...
        'Marker','d','MarkerEdgeColor','none',...
        'MarkerFaceColor','k','MarkerSize',4);
end
set(gca,'ylim',[-100 100]);