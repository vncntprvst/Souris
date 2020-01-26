function NBC_Plots_Overview(whiskerAngle,whiskerPhase,whiskerSetPoint,ephys,saveFig)
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

%% allocate
rasters=ephys.rasters(ephys.selectedUnits,:);
scaleF=ceil(max(whiskerAngle));

%% open figure
figure('color','white'); 
angleAxH = axes('Position',[0.1 0.5 0.8 0.4]); hold on;
spikesAxH = axes('Position',[0.1 0.1 0.8 0.4]); hold on;

%% plot whisker angle and setPoint
axes(angleAxH)
plot(whiskerAngle,'color',cmap(1,:),'linewidth',1.2);
plot(whiskerSetPoint,'color',[cmap(2,:) 0.5],'linewidth',1.2);

% plot Phase
patchColor=[180, 225, 228]./255;
protractPeriods=bwconncomp(whiskerPhase<=0);
proPIdx=cellfun(@(proP) [proP(1) proP(1) proP(end) proP(end)],...
    protractPeriods.PixelIdxList,'un',0); proPIdx=[proPIdx{:}];
patch('Faces',reshape(1:protractPeriods.NumObjects*4,[4,protractPeriods.NumObjects])',...
    'Vertices',[proPIdx',...
    repmat([0;scaleF;scaleF;0],protractPeriods.NumObjects,1)],...
    'FaceVertexCData',repmat([0;0;6;6],protractPeriods.NumObjects,1),...
    'FaceColor',patchColor,'EdgeColor','none','FaceAlpha',0.5);           
set(gca,'xtick',[],'xcolor','none','TickDir','out','box','off');
ylabel('whisker angle')
legend('Whisker angle','Whisker setpoint','Protraction phase')
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
set(gca,'xtick',0:1000:size(rasters,2),...
        'xticklabel',(0:1000:size(rasters,2))/1000,...
        'ytick',0.5:size(rasters,1)-0.5,...
        'yticklabel',1:size(rasters,1),...
        'TickDir','out','box','off');
xlabel('Time (s)')
ylabel('unit number')

%link axes
linkaxes([angleAxH,spikesAxH],'x')

%% save figure
if saveFig
    savefig(phaseTuningSummaryFig,[ephys.recName '_overview.fig'])
    saveas(phaseTuningSummaryFig,[ephys.recName '_overview.png'])
end