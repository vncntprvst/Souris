function NBC_Plots_PhaseTuning(whiskingPhase_ms,whiskingPeriod,spikeRasters_ms,saveFig,recName)
% Polar plots of the spectral coherence between spiking activity and vibrissa motion 
% at the peak frequency of whisking;




if nargin<5
    recName='PhaseTuning_polarPlot';
end
phaseTuningSummaryFig=figure('position',[969    49   944   948],'name',recName);
colormapSeed=lines;
cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];

for unitNum=1:size(spikeRasters_ms,1)
    % get spikes
    unitSpikes=spikeRasters_ms(unitNum,:);
    % restrict to whisking periods
    unitSpikesWP=unitSpikes(whiskingPeriod);
    whiskingPhaseWP_ms=whiskingPhase_ms(1,whiskingPeriod);
    spikeOnWPhase=whiskingPhaseWP_ms(logical(unitSpikesWP));
    spikeOnWPhase=spikeOnWPhase+pi;
    %     figure;
    subplot(ceil(size(spikeRasters_ms,1)/4),4,unitNum)
    pph=polarhistogram(spikeOnWPhase,20,'Displaystyle','stairs',...
        'Normalization','count','LineWidth',2);
    %mean angle value
    %     rad2deg(circ_mean((spikeOnWPhase')))+360
    phaseStats=circ_stats(spikeOnWPhase');
    if phaseStats.skewness>0.02 && circ_rtest(spikeOnWPhase')<0.05
        phaseTuning=rad2deg(phaseStats.mean);
        pph.DisplayStyle='bar'; pph.EdgeAlpha=0;
        if phaseTuning>0 && phaseTuning<180
            pph.EdgeColor='r';pph.FaceColor='r';
        elseif phaseTuning>180 && phaseTuning~=0
            pph.EdgeColor='g';pph.FaceColor='g';
        else
            pph.EdgeColor='k';pph.FaceColor='k';
        end
    else
        phaseTuning=[];
    end
end
if saveFig
    savefig(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.fig'])
    saveas(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.png'])
end
% unitsOfInterestIdx=[5,1,6,10];unitsOfInterest=keepUnits(unitsOfInterestIdx);
