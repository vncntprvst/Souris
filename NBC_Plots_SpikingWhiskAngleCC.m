function NBC_Plots_SpikingWhiskAngleCC(periodBehavData,whiskingPhase,...
    whiskingPeriodIdx,spikeRasters,SDFs,saveFig,recName)
if nargin<6
    recName='SpikingWhiskAngleCC';
end
% only consider time periods when whisking occurs angle and phase
timeVector=1:numel(periodBehavData(1,:)); timeVector(~whiskingPeriodIdx)=NaN;
% plot(timeVector,periodBehavData_ms-median(periodBehavData_ms))
oscillationPattern=cos(whiskingPhase(1,:));
% cross correlation
SpikeWAngleCorrFigure=figure('name',recName); %'position',[969    49   944   948]
for unitNum=1:size(SDFs,1) %find(keepUnits==15); %bestUnit=2; %4;
    unitSDF=SDFs(unitNum,:); % unitSpikes=spikeRasters_ms(unitNum,:);
    [acor,lag] = xcorr(unitSDF(whiskingPeriodIdx),...
        oscillationPattern(whiskingPeriodIdx),150,'coeff');
    figure(SpikeWAngleCorrFigure)
    if size(spikeRasters,1)/4 > 1
        subplot(ceil(size(spikeRasters,1)/4),4,unitNum);
    end
    ccph=plot(lag,acor,'color','k','LineWidth',2);set(gca,'ylim',[-0.4 0.4]); %xlabel('Lag (ms)')
%     title({['Cross correlation for vIRt unit ' num2str(keepUnits(unitNum))];...
%         'Spike density function vs. Whisking angle'});
    %retraction unit if it shoots up around 0, protraction if before
    if abs(acor(50:150))>0.2
        ccph.Color='r'; %cmap(keepUnits(unitNum),:);
        burstyWRCell=true;
    else
        burstyWRCell=false;
    end
    if burstyWRCell
%         % plot spike times, angle and phase
%         figure; hold on
%         plot(timeVector,BP_periodBehavData_ms(1,:)); %-median(BP_periodBehavData_ms));
%         plot(timeVector,whiskingPhase_ms(1,:))
%         % for unitNum=1:size(spikeRasters_ms,1) %find(keepUnits==15)
%         unitSpikes=spikeRasters_ms(unitNum,:);
%         %     unitSpikesExcerpt=unitSpikes(whiskingPeriod);
%         unitSpikes(isnan(timeVector) | unitSpikes==0)=nan;
% %         unitSpikes(unitSpikes==1)=unitNum;
%         plot(timeVector,unitSpikes,'d')
%         set(gca, 'xlim', [longestWhiskingPeriodIdx(1),longestWhiskingPeriodIdx(end)])
%         if isempty(phaseTuning)
%             tuningLabel='none';
%         else
%            tuningLabel= num2str(phaseTuning);
%         end
%         title([ephys.recName ' unit ' num2str(keepUnits(unitNum))  ' Tuning ' tuningLabel],'interpreter','latex')
%         savefig(gcf,[ephys.recName '_BurstUnit' keepUnits(unitNum) '.fig'])
%         saveas(gcf,[ephys.recName '_BurstUnit' keepUnits(unitNum) '.png'])
%         close(gcf);
    end
end
if saveFig
    savefig(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.fig'])
    saveas(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.png'])
end

