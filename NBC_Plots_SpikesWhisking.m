function NBC_Plots_SpikesWhisking(whiskingAngle,whiskingVelocity,whiskingPhase,whiskingEpochsIdx,spikeRasters)

% medWhiskAngleVal=median(periodBehavData_ms(1,:));
cmap=lines;
cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];
    sigma=2; %2ms smoothing
    sizeK = 6*sigma;
    
whiskingVelocity=diff(whiskingAngle);
whiskingVelocity=[whiskingVelocity(1);whiskingVelocity];
whiskingAngle=whiskingAngle(whiskingEpochsIdx);
whiskingVelocity=whiskingVelocity(whiskingEpochsIdx);
whiskingPhase=whiskingPhase(whiskingEpochsIdx);
spikeRasters=spikeRasters(:,whiskingEpochsIdx);

figure; hold on
plot(whiskingAngle); plot(whiskingVelocity*std(whiskingAngle)+mean(whiskingAngle))
plot(whiskingVelocity*std(whiskingAngle))

velCycles=find(diff(whiskingVelocity>0)==1)+1; %>0 <0
velCycles=velCycles(find(diff(velCycles)>50)+1); %keep under 20Hz (cleanest)
whiskingBoutExcerpts=nan(numel(velCycles)-1,50+(2*sizeK));
allRasters=nan(size(spikeRasters,1),numel(velCycles)-1,50+(2*sizeK));
for velCycleNum=1:numel(velCycles)-1
    whiskingBout=whiskingAngle(velCycles(velCycleNum):velCycles(velCycleNum+1));
    whiskingBoutPeak=find(whiskingBout==max(whiskingBout),1); %max min
    whiskingBoutExcerptIndex=(whiskingBoutPeak-25-sizeK:whiskingBoutPeak+24+sizeK)+velCycles(velCycleNum);
    whiskingBoutExcerpts(velCycleNum,:)=whiskingAngle(whiskingBoutExcerptIndex);
    allRasters(:,velCycleNum,:)=spikeRasters(:,whiskingBoutExcerptIndex);
end

meanWhiskingProfile=mean(whiskingBoutExcerpts); %mean(zscore(whiskingBoutExcerpts,[],2)); %
figure; plot(meanWhiskingProfile)

for cellNum=1:size(spikeRasters,1)
    cellRasters=squeeze(allRasters(cellNum,:,:));
    
    %% find raster indices
    [indy, indx] = ind2sub(size(cellRasters),find(cellRasters)); %find row and column coordinates of spikes
    %% calculate spike density function
    widthK = linspace(-sizeK / 2, sizeK / 2, sizeK);
    gaussFilter = exp(-widthK .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
%     gaussFilter(widthK<0)=0; % causal kernel
    convTrace = conv (mean(cellRasters), gaussFilter, 'same')*1000;
    convTrace = convTrace(1+sizeK:end-sizeK);
%     convTrace =EphysFun.MakeSDF(mean(cellRasters));
    %% plot
    figure;
    subplot(3,1,1)
    plot(meanWhiskingProfile)
    subplot(3,1,2)
    line([indx';indx'],[indy';indy'+1],'color','k',...
        'LineStyle','-','LineWidth',1.8); % plot rasters
    xlim([1+sizeK 50+sizeK])
%     yLims=get(gca,'YLim');
%     patch([24 24 25 25],[min(yLims) max(yLims) max(yLims) min(yLims)],...
%         'k','FaceAlpha',0.2,'EdgeColor','none');
    currAx = gca; % current axes
    currAx.XLim=[1 size(cellRasters,2)];
    currAx.YDir='reverse';
    currAx.Visible='off';
    subplot(3,1,3)
    plot(convTrace,'Color',cmap(1,:),'Linewidth',2); %[0.5 0.3 0.6]
%     axis tight;
    yLims=get(gca,'YLim');
    set(gca,'YLim',[0 max(yLims)]);
%     patch([10 10 60 60],[min(yLims) max(yLims) max(yLims) min(yLims)],...
%         'k','FaceAlpha',0.2,'EdgeColor','none');
    currAx = gca; % current axes
    currAx.FontSize = 12;
    currAx.TickDir = 'out';
    currAx.TickLength = [0.015 0.015];
    currAx.XTick=5:5:45;
    currAx.XTickLabels=[-20:5:20];
    currAx.FontName='Arial';
    currAx.Box='off';
    
end

%%



plotShift=max(whiskingAngle(1,:))+std(whiskingAngle(1,:))/2;
figure; hold on
% for whiskNum=1:3 %size(whiskingAngle,1)-1
%     plot(whiskingAngle(whiskNum,:),'Color',[whiskNum*0.1 whiskNum*0.1 whiskNum*0.1 0.5])
% end
whiskingEpochs=nan(1,numel(whiskingEpochsIdx));whiskingEpochs(whiskingEpochsIdx)=1;
plot(whiskingEpochs+plotShift,'linewidth',2,'Color','b')
% plot(whiskerTraces_breathFreq_ms(1,:)-plotShift)
plot(whiskingPhase(1,:));
plot(whiskingAngle(1,:),'Color','k','linewidth',1.5)

% plotShift=min(whiskingAngle(1,:))-std(whiskingAngle(1,:))/2;

if ~iscell(spikeRasters)
    for spkNum=1:size(spikeRasters,1)
        %% bars (pretty heavy on memory)
        %     %find row and column coordinates of spikes
        %     [indy, indx] = ind2sub(size(spikeRasters{1}(spkNum,:)),find(spikeRasters{1}(spkNum,:)));
        %     % plot rasters
        %     plot([indx;indx],[indy;indy+1]+min(whiskingAngle(1,:))-spkNum,'color',cmap(spkNum,:),'LineStyle','-');
        %% diamonds
        plot(find(spikeRasters(spkNum,:)),...
            ones(1,numel(find(spikeRasters(spkNum,:))))*...
            plotShift-spkNum*2,'LineStyle','none',...
            'Marker','d','MarkerEdgeColor','none',...
            'MarkerFaceColor',cmap(spkNum,:),'MarkerSize',4);
    end
else
    for spkNum=1:size(spikeRasters{1},1)
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
end