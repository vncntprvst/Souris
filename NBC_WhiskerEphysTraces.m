function NBC_WhiskerEphysTraces(whiskers,ephys,opt)
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% if there's data about contralateral whiskers, keep best
switch opt.xpType
    case 'GFE3' % for GFE3 experiments
        opoWhiskerIdx=[whiskers.bestWhisker]==true & contains({whiskers.side},'left');
        opoWhisker=whiskers(opoWhiskerIdx);
        whisker=whiskers([whiskers.bestWhisker]==true & contains({whiskers.side},'right'));
        wLabels={'Contra whisker angle','Ipsi whisker angle'};
    case 'asymmetry'
        opoWhiskerIdx=[whiskers.bestWhisker]==true & contains({whiskers.side},'right');
        opoWhisker=whiskers(opoWhiskerIdx);
        whisker=whiskers([whiskers.bestWhisker]==true & contains({whiskers.side},'left'));
        wLabels={'Ipsi whisker angle','Contra whisker angle'};
    otherwise
        whisker=whiskers([whiskers.bestWhisker]==true & contains({whiskers.side},'left'));
end

% find base setpoint and redress values if needed
baseSP=mode(round(whisker.setPoint/10)*10);
if baseSP<0
    whisker.setPoint=whisker.setPoint+180;
    whisker.angle=whisker.angle+180;
end

%% open figure
figure('color','white','position',[62,167,1307,727],...
    'name',['Overview_' ephys.recInfo.baseName],'NumberTitle', 'off');
angleAxH = axes('Position',[0.1 0.5 0.8 0.4]); hold on;
ephysAxH = axes('Position',[0.1 0.1 0.8 0.4]); hold on;

%% plot whisker angle and setpoint
axes(angleAxH)
pH{1}=plot(whisker.timestamp,whisker.angle,'color',cmap(1,:),'linewidth',1.2);
if exist('opoWhisker','var') && ~isempty(opoWhisker)
    pH{2}=plot(opoWhisker.timestamp,opoWhisker.angle,'color',cmap(4,:),'linewidth',1.2);
else
    pH{2}=plot(whisker.timestamp,whisker.setPoint,'color',[cmap(2,:) 0.5],'linewidth',1.2);
end

set(gca,'xtick',[],'xcolor','none','TickDir','out','box','off');
ylabel('whisker angle')
if exist('opoWhisker','var') && ~isempty(opoWhisker)
    legend([pH{:}],{wLabels{1},wLabels{2},'Protraction phase'},'location','southeast')
else
    legend([pH{:}],{'Whisker angle','Whisker setpoint','Protraction phase'},'location','southeast')
end
legend('boxoff')

%% plot ephys trace
axes(ephysAxH)
plot([1:size(ephys.traces,2)]/30000,ephys.traces(opt.chNum,:),'color','k','linewidth',1.2);
set(gca,'TickDir','out','box','off');
xlabel('Time (s)')


%link axes
linkaxes([angleAxH,ephysAxH],'x')

% zoom in whisking epochs
if opt.zoomin
    wEpochs=bwconncomp(whiskingEpochs);
    viewStartPoint=round(wEpochs.PixelIdxList{1, 9}(1)/1000)*1000;
    set(gca,'xlim',[viewStartPoint viewStartPoint+2000])
end
%% save figure
if opt.saveFig
    savefig(gcf,['Overview_Whole_' ephys.recInfo.sessionName '.fig'])
    %     saveas(gcf,['Overview_Whole_' ephys.recName '.png'])
end