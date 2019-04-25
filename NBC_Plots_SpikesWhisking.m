function NBC_Plots_SpikesWhisking(BP_periodBehavData_ms,whiskingPhase_ms,spikeRasters_ms)

% medWhiskAngleVal=median(periodBehavData_ms(1,:));

figure; hold on
plot(BP_periodBehavData_ms(1,:))
plot(whiskingPhase_ms(1,:))
for spkNum=1:size(spikeRasters_ms,1)
    plot(find(spikeRasters_ms(spkNum,:)),...
        ones(1,numel(find(spikeRasters_ms(spkNum,:))))*...
        (spkNum-round(size(spikeRasters_ms,1)/2)),'d'); % +medWhiskAngleVal
end
set(gca,'ylim',[-100 100]);