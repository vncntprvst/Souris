function phaseTuning=NBC_Plots_PhaseTuning(whiskingPhase,breathingPhase,...
    whiskingEpochs,spikeRasters,saveFig,recName)
% Polar plots of the spectral coherence between spiking activity and vibrissa motion
% at the peak frequency of whisking;

if nargin<5
    recName='PhaseTuning_polarPlot';
end
phaseTuningSummaryFig=figure('name',recName); %'position',[969    49   944   948],
% colormapSeed=lines;
% cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
phaseTuning=nan(size(spikeRasters,1),2);

for unitNum=1:size(spikeRasters,1)
    % get spikes
    unitSpikes=spikeRasters(unitNum,:);
    % restrict to whisking periods
    unitSpikesWE=unitSpikes(whiskingEpochs);
    whiskingPhaseWE=whiskingPhase(1,whiskingEpochs);
    spikeOnWPhase=whiskingPhaseWE(logical(unitSpikesWE));
    spikeOnWPhase=spikeOnWPhase+pi;
    spikeOnWPhaseStats=circ_stats(spikeOnWPhase');
    
    breathingPhaseWE=breathingPhase(1,whiskingEpochs);
    spikeOnBPhase=breathingPhaseWE(logical(unitSpikesWE));
    spikeOnBPhase=spikeOnBPhase+pi;
    spikeOnBPhaseStats=circ_stats(spikeOnBPhase');
    
    if ((spikeOnWPhaseStats.kurtosis>0.04 || spikeOnWPhaseStats.skewness<-0.02) &&...
            circ_rtest(spikeOnWPhase')<0.05) || ...
            ((spikeOnBPhaseStats.kurtosis>0.04 || spikeOnBPhaseStats.skewness<-0.02) &&...
            circ_rtest(spikeOnBPhase')<0.05)
        phaseTuning(unitNum,1)=rad2deg(spikeOnWPhaseStats.mean);
        phaseTuning(unitNum,2)=rad2deg(spikeOnBPhaseStats.mean);
    end
end
tunedIdx=find(~isnan(phaseTuning(:,1)));
for plotNum=1:numel(tunedIdx)
    unitSpikes=spikeRasters(tunedIdx(plotNum),:);
    % restrict to whisking periods
    unitSpikesWE=unitSpikes(whiskingEpochs);
    whiskingPhaseWE=whiskingPhase(1,whiskingEpochs);
    spikeOnWPhase=whiskingPhaseWE(logical(unitSpikesWE));
    spikeOnWPhase=spikeOnWPhase+pi;
    breathingPhaseWE=breathingPhase(1,whiskingEpochs);
    spikeOnBPhase=breathingPhaseWE(logical(unitSpikesWE));
    spikeOnBPhase=spikeOnBPhase+pi;

        %     figure;
        if size(spikeRasters,1)>1
            subplot(ceil(numel(tunedIdx)/4),4,plotNum)
        end
        wpph=polarhistogram(spikeOnWPhase,20,'Displaystyle','bar',...
            'Normalization','count','LineWidth',2); hold on
        wpph.EdgeAlpha=0; %wpph.DisplayStyle='bar';
        if phaseTuning(tunedIdx(plotNum),1)>0
            if phaseTuning(tunedIdx(plotNum),1)<90
                wpph.EdgeColor=[1 0 0];wpph.FaceColor=[0.85 0.15 0.15];
            else
                wpph.EdgeColor=[1 0 0];wpph.FaceColor=[0.7 0.3 0.2];
            end
        elseif phaseTuning(tunedIdx(plotNum),1)<0
            if phaseTuning(tunedIdx(plotNum),1)>-90
                wpph.EdgeColor=[0 1 0];wpph.FaceColor=[0.27 0.6 0.2];
            else
                wpph.EdgeColor=[1 0 0];wpph.FaceColor=[0.3 0.7 0.2];
            end
        else
            wpph.EdgeColor='k';wpph.FaceColor='k'; wpph.EdgeAlpha=0.5;
        end
        
        bpph=polarhistogram(spikeOnBPhase,20,'Displaystyle','stairs',...
            'Normalization','count','LineWidth',0.5); 
        bpph.EdgeColor='k';  bpph.EdgeAlpha=0.3;
%       bpph.FaceColor='none';
%         
%         bpph.EdgeAlpha=0;
%         if phaseTuning(tunedIdx(plotNum),1)>0
%             if phaseTuning(tunedIdx(plotNum),1)<90
%                 bpph.EdgeColor=[1 0 0];bpph.FaceColor=[0.85 0.15 0.15];
%             else
%                 bpph.EdgeColor=[1 0 0];bpph.FaceColor=[0.7 0.3 0.2];
%             end
%         elseif phaseTuning(tunedIdx(plotNum),1)<0
%             if phaseTuning(tunedIdx(plotNum),1)>-90
%                 bpph.EdgeColor=[0 1 0];bpph.FaceColor=[0.27 0.6 0.2];
%             else
%                 bpph.EdgeColor=[1 0 0];bpph.FaceColor=[0.3 0.7 0.2];
%             end
%         else
%             bpph.EdgeColor='k';bpph.FaceColor='k'; bpph.EdgeAlpha=0.5;
%         end
% 
        hold off
end
if saveFig
    savefig(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.fig'])
    saveas(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.png'])
end
% unitsOfInterestIdx=[5,1,6,10];unitsOfInterest=keepUnits(unitsOfInterestIdx);
