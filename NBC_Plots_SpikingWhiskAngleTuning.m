function NBC_Plots_SpikingWhiskAngleTuning(whiskerAngleData,whiskingPhaseData,...
    whiskingEpochs,spikeRasters,saveFig,recName)
if nargin<6
    recName='SpikingWhiskAngleCC';
end
SpikeWAngleCorrFigure=figure('name',recName); %'position',[969    49   944   948]
for unitNum=1:size(spikeRasters,1)
    % get spikes
    unitSpikes=spikeRasters(unitNum,:);
    % restrict to whisking periods
    unitSpikesWP=unitSpikes(whiskingEpochs);
    whiskingAngleWP=whiskerAngleData(1,whiskingEpochs);
    whiskingPhaseWP=whiskingPhaseData(1,whiskingEpochs);
    spikeOnAngle{1}=whiskingAngleWP(logical(unitSpikesWP) & logical(whiskingPhaseWP>0));
    spikeOnAngle{2}=whiskingAngleWP(logical(unitSpikesWP) & logical(whiskingPhaseWP<0));
%     if size(spikeRasters,1)>1
    subplot(ceil(size(spikeRasters,1)/4),4,unitNum); hold on
%     end
    histogram(spikeOnAngle{1},20,'Displaystyle','stairs',...
        'Normalization','count','LineWidth',2,'FaceColor','none','EdgeColor',[0.2 0.5 0.7]);
%     subplot(ceil(size(spikeRasters,1)/4)*2,4,unitNum*2)
    histogram(spikeOnAngle{2},20,'Displaystyle','stairs',...
        'Normalization','count','LineWidth',2,'FaceColor','none','EdgeColor',[0.7 0.5 0.2]);
    %% add randomized values
    randomEpochsSample=randsample(1:numel(whiskingPhaseData(1,:)),sum(whiskingEpochs));
    randomEpochs=false(1,numel(whiskingEpochs)); randomEpochs(randomEpochsSample)=true;
    unitSpikesRand=unitSpikes(randomEpochs);
    whiskingAngleRand=whiskerAngleData(1,randomEpochs);
    spikeOnRandAngle=whiskingAngleRand(logical(unitSpikesRand));
        histogram(spikeOnRandAngle,20,...
        'Normalization','count','FaceColor',[0 0 0],...
        'EdgeColor','none','FaceAlpha',0.1);

end
if saveFig
    savefig(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.fig'])
    saveas(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.png'])
end

