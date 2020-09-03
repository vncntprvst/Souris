function NBC_Plots_Overview(whiskerAngle,whiskerPhase,whiskerSetPoint,whiskingEpochs,ephys,TTLtimes,saveFig)
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% find base sertpoint and redress values if needed 
baseSP=mode(round(whiskerSetPoint/10)*10);
if baseSP<0
whiskerSetPoint=whiskerSetPoint+180;
whiskerAngle=whiskerAngle+180;
end

%% allocate
rasters=ephys.rasters(ephys.selectedUnits,:);
scaleF=ceil(max(whiskerAngle));
if ~exist('TTLtimes','var') | isempty(TTLtimes)
    TTLtimes=[]; else; pulseDur=mode(diff(TTLtimes)); end

%% open figure
figure('color','white','position',[62,167,1307,727]); 
angleAxH = axes('Position',[0.1 0.5 0.8 0.4]); hold on;
spikesAxH = axes('Position',[0.1 0.1 0.8 0.4]); hold on;

%% plot whisker angle and setPoint
axes(angleAxH)
pH{1}=plot(whiskerAngle,'color',cmap(1,:),'linewidth',1.2);
pH{2}=plot(whiskerSetPoint,'color',[cmap(2,:) 0.5],'linewidth',1.2);
for pulseNum=1:size(TTLtimes,2)
    patch([TTLtimes(1,pulseNum), TTLtimes(1,pulseNum),...
        TTLtimes(1,pulseNum)+pulseDur, TTLtimes(1,pulseNum)+pulseDur], ...
        [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
end

% plot Phase
patchColor=[222, 207, 103]./255; %180, 225, 228
protractPeriods=bwconncomp(whiskerPhase<=0 & whiskingEpochs);
proPIdx=cellfun(@(proP) [proP(1) proP(1) proP(end) proP(end)],...
    protractPeriods.PixelIdxList,'un',0); proPIdx=[proPIdx{:}];
pH{3}=patch('Faces',reshape(1:protractPeriods.NumObjects*4,[4,protractPeriods.NumObjects])',...
    'Vertices',[proPIdx',...
    repmat([scaleF-10;scaleF;scaleF;scaleF-10],protractPeriods.NumObjects,1)],...
    'FaceVertexCData',repmat([0;0;6;6],protractPeriods.NumObjects,1),...
    'FaceColor',patchColor,'EdgeColor','none','FaceAlpha',0.5);           
set(gca,'xtick',[],'xcolor','none','TickDir','out','box','off');
ylabel('whisker angle')
legend([pH{:}],{'Whisker angle','Whisker setpoint','Protraction phase'},'location','southeast')
legend('boxoff')

%% plot rasters 
axes(spikesAxH)
EphysFun.PlotRaster(rasters,'lines',[],'k')

% plot Phase
scaleF=size(rasters,1);
patch('Faces',reshape(1:protractPeriods.NumObjects*4,[4,protractPeriods.NumObjects])',...
    'Vertices',[proPIdx',...
    repmat([0;scaleF;scaleF;0],protractPeriods.NumObjects,1)],...
    'FaceVertexCData',repmat([0;0;6;6],protractPeriods.NumObjects,1),...
    'FaceColor',patchColor,'EdgeColor','none','FaceAlpha',0.5);   
for pulseNum=1:size(TTLtimes,2)
    patch(spikesAxH,[TTLtimes(1,pulseNum), TTLtimes(1,pulseNum),...
        TTLtimes(1,pulseNum)+pulseDur, TTLtimes(1,pulseNum)+pulseDur], ...
        [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
end
set(gca,'xtick',0:1000:size(rasters,2),...
        'xticklabel',(0:1000:size(rasters,2))/1000,...
        'ytick',0.5:size(rasters,1)-0.5,...
        'yticklabel',1:size(rasters,1),...
        'TickDir','out','box','off');
if isfield(ephys,'selectedUnits')
    set(gca,'yticklabel',ephys.selectedUnits);
end
xlabel('Time (s)')
ylabel('unit number')

% axes(spikesAxH)
% % yyaxis right
% axes('YAxisLocation','right')
% ylabel('X/Y Probe Geometry')
% set(gca,'ytick',0.5:size(rasters,1)-0.5,...
%         'yticklabel',ephys.unitCoordinates(:,2),...
%         'TickDir','out','box','off');

%link axes
linkaxes([angleAxH,spikesAxH],'x')

% zoom in whisking epochs
wEpochs=bwconncomp(whiskingEpochs);
viewStartPoint=round(wEpochs.PixelIdxList{1, 9}(1)/1000)*1000;
set(gca,'xlim',[viewStartPoint viewStartPoint+2000])

%% save figure
if saveFig
    savefig(gcf,['Overview_Whole_' ephys.recInfo.sessionName '.fig']) 
%     saveas(gcf,['Overview_Whole_' ephys.recName '.png'])
end