function phaseTuning=NBC_Plots_PhaseTuning_PolarPlots(whiskingPhase,...
    dataMask,spikeRate,saveFig,recName)
% whiskingPhase in radians
% spikeRate in Hz

if nargin<5
    recName='PhaseTuning_polarPlot';
end
if ~contains(recName,'subplot')
    phaseTuningSummaryFig=figure('name',recName); %'position',[969    49   944   948],
end
numBins=30;
whiskingPhase = whiskingPhase(dataMask);
edges = linspace(min(whiskingPhase), max(whiskingPhase), numBins+1);
centers = mean([ edges(1:end-1); edges(2:end) ]);
[ ~, binInd ] = histc(whiskingPhase, edges);
    
colormapSeed=lines;
cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
phaseTuning=nan(size(spikeRate,1),1);

%% comparing with UnitExplorer GetTuning

for unitNum=1:size(spikeRate,1)
    clearvars unitSpikeRate 
    % Data masking
    unitSpikeRate= spikeRate(unitNum,dataMask);
    [binMeanSpikeRate,binSteSpikeRate]=deal(nan(numBins,1));
    % Bining firing rate by stimulus value range
    for binNum = length(edges)-1 : -1 : 1
        ratesVect = unitSpikeRate(binInd == binNum);
        numSample = numel(ratesVect);
        if numSample == 0
            meanSpikeRate = 0;
            steSpikeRate = 0;
        else
            meanSpikeRate = nanmean(ratesVect);
            steSpikeRate = MMath.StandardError(ratesVect);
        end
        binMeanSpikeRate(binNum) = meanSpikeRate;
%         binSteSpikeRate(binNum) = steSpikeRate;
    end
    
%     subplot(1,2,1)
%     set(gcf,'Color','w')
%     box off
%     shadedErrorBar(centers, binMeanSpikeRate,binSteSpikeRate, 'lineprops','k');
%     hold on;
%     plot(centers, binMeanSpikeRate, 'k', 'LineWidth',2);
%     xlabel({'Phase  (rad)'; '0 = Protraction'}, 'FontSize', 18);
%     ylabel('Firing rate (Spk/s)', 'FontSize', 18);
%     axis tight
%     yl = ylim;
%     ylim([0 yl(2)]);
%     set(gca,'TickDir','out')
%     box off
%     xlim([-pi pi])
%     set(gca,'xdir', 'reverse'); %, 'ydir', 'reverse')
    
    
    %% convert to thetas
    thetas=cell(numel(centers),1);
    for binNum=1:numel(centers)
        thetas{binNum}=ones(round(binMeanSpikeRate(binNum)),1)*centers(binNum);
    end
    thetas=vertcat(thetas{:});
    if isempty(thetas)
        disp('not enough spikes')
        continue
    end
    % stats
    phaseStats=circ_stats(thetas);
    if  circ_rtest(thetas)<0.05 %((phaseStats.kurtosis>0.04 || phaseStats.skewness<-0.02) || ...
        phaseTuning(unitNum,1)=rad2deg(phaseStats.mean);
    end
    if ~isnan(phaseTuning(unitNum,1))
        phEdgeColor=cmap(4,:);phFaceColor=cmap(4,:);
%         if phaseTuning(unitNum,1)<90
%             phEdgeColor=[1 0 0];phFaceColor=[0.85 0.15 0.15];
%         else
%             phEdgeColor=[1 0 0];phFaceColor=[0.7 0.3 0.2];
%         end
%     elseif phaseTuning(unitNum,1)<0
%         if phaseTuning(unitNum,1)>-90
%             phEdgeColor=[0 1 0];phFaceColor=[0.27 0.6 0.2];
%         else
%             phEdgeColor=[1 0 0];phFaceColor=[0.3 0.7 0.2];
%         end
    else
        phEdgeColor='k';phFaceColor='k'; %EdgeAlpha=0.5;
    end
    
    if size(spikeRate,1)>1
        subplot(ceil(size(spikeRate,1)/4),4,unitNum)
    end
    
%     rightPanel=subplot(1,2,2);
%     polaraxes(rightPanel);
    polarhistogram(thetas,binNum,'Displaystyle','bar',...
        'Normalization','count','LineWidth',2,...
        'EdgeColor',phEdgeColor,'FaceColor',phFaceColor,...
        'EdgeAlpha',0);
    paH = gca;
    % For Kyle's convention: 
    paH.ThetaZeroLocation='left';
    paH.ThetaDir = 'clockwise';
        paH.ThetaTickLabel={'mP','','','P','','',...
        'mR','','','R','',''};
%     paH.ThetaTickLabel={'Protracted','','','Protraction','','',...
%         'Retracted','','','Retraction','',''};
end

% % colormapSeed=lines;
% % cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
% phaseTuning=nan(size(spikeRasters,1),1);
%
% for unitNum=1:size(spikeRasters,1)
%     % get spikes
%     unitSpikes=spikeRasters(unitNum,:);
%     % restrict to whisking periods
%     unitSpikesWE=unitSpikes(dataMask);
%     whiskingPhaseWE=whiskingPhase(1,dataMask);
%     spikeOnWPhase=whiskingPhaseWE(logical(unitSpikesWE));
%     spikeOnWPhase=spikeOnWPhase+pi;
%     spikeOnWPhaseStats=circ_stats(spikeOnWPhase');
%
% %     breathingPhaseWE=breathingPhase(1,dataMask);
% %     spikeOnBPhase=breathingPhaseWE(logical(unitSpikesWE));
% %     spikeOnBPhase=spikeOnBPhase+pi;
% %     spikeOnBPhaseStats=circ_stats(spikeOnBPhase');
%
% %     if ((spikeOnWPhaseStats.kurtosis>0.04 || spikeOnWPhaseStats.skewness<-0.02) &&...
% %             circ_rtest(spikeOnWPhase')<0.05) %|| ...
% % %             ((spikeOnBPhaseStats.kurtosis>0.04 || spikeOnBPhaseStats.skewness<-0.02) &&...
% % %             circ_rtest(spikeOnBPhase')<0.05)
%         phaseTuning(unitNum,1)=rad2deg(spikeOnWPhaseStats.mean);
% %         phaseTuning(unitNum,2)=rad2deg(spikeOnBPhaseStats.mean);
% %     end
% end
% tunedIdx=find(~isnan(phaseTuning(:,1)));
% for plotNum=1:numel(tunedIdx)
%     unitSpikes=spikeRasters(tunedIdx(plotNum),:);
%     % restrict to whisking periods
%     unitSpikesWE=unitSpikes(dataMask);
%     whiskingPhaseWE=whiskingPhase(1,dataMask);
%     spikeOnWPhase=whiskingPhaseWE(logical(unitSpikesWE));
%     spikeOnWPhase=spikeOnWPhase+pi;
% %     breathingPhaseWE=breathingPhase(1,dataMask);
% %     spikeOnBPhase=breathingPhaseWE(logical(unitSpikesWE));
% %     spikeOnBPhase=spikeOnBPhase+pi;
%
%         %     figure;
%         if size(spikeRasters,1)>1
%             subplot(ceil(numel(tunedIdx)/4),4,plotNum)
%         end
%         wpph=polarhistogram(spikeOnWPhase,numBins,'Displaystyle','bar',...
%             'Normalization','count','LineWidth',2); %hold on
%         wpph.EdgeAlpha=0; %wpph.DisplayStyle='bar';
%         if phaseTuning(tunedIdx(plotNum),1)>0
%             if phaseTuning(tunedIdx(plotNum),1)<90
%                 wpph.phEdgeColor=[1 0 0];wpph.phFaceColor=[0.85 0.15 0.15];
%             else
%                 wpph.phEdgeColor=[1 0 0];wpph.phFaceColor=[0.7 0.3 0.2];
%             end
%         elseif phaseTuning(tunedIdx(plotNum),1)<0
%             if phaseTuning(tunedIdx(plotNum),1)>-90
%                 wpph.phEdgeColor=[0 1 0];wpph.phFaceColor=[0.27 0.6 0.2];
%             else
%                 wpph.phEdgeColor=[1 0 0];wpph.phFaceColor=[0.3 0.7 0.2];
%             end
%         else
%             wpph.phEdgeColor='k';wpph.phFaceColor='k'; wpph.EdgeAlpha=0.5;
%         end
%
% %         bpph=polarhistogram(spikeOnBPhase,20,'Displaystyle','stairs',...
% %             'Normalization','count','LineWidth',0.5);
% %         bpph.phEdgeColor='k';  bpph.EdgeAlpha=0.3;
% %       bpph.phFaceColor='none';
% %
% %         bpph.EdgeAlpha=0;
% %         if phaseTuning(tunedIdx(plotNum),1)>0
% %             if phaseTuning(tunedIdx(plotNum),1)<90
% %                 bpph.phEdgeColor=[1 0 0];bpph.phFaceColor=[0.85 0.15 0.15];
% %             else
% %                 bpph.phEdgeColor=[1 0 0];bpph.phFaceColor=[0.7 0.3 0.2];
% %             end
% %         elseif phaseTuning(tunedIdx(plotNum),1)<0
% %             if phaseTuning(tunedIdx(plotNum),1)>-90
% %                 bpph.phEdgeColor=[0 1 0];bpph.phFaceColor=[0.27 0.6 0.2];
% %             else
% %                 bpph.phEdgeColor=[1 0 0];bpph.phFaceColor=[0.3 0.7 0.2];
% %             end
% %         else
% %             bpph.phEdgeColor='k';bpph.phFaceColor='k'; bpph.EdgeAlpha=0.5;
% %         end
% %
%         hold off
% end
if saveFig
    savefig(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.fig'])
    saveas(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.png'])
end
% % unitsOfInterestIdx=[5,1,6,10];unitsOfInterest=keepUnits(unitsOfInterestIdx);
