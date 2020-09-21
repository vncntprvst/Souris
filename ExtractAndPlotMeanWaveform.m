figure; hold on 

for wfNum=1:32
    waveForms=NaN(size(spikesTimes,1),100);
%         electrodesId=unique(spikes.preferredElectrode);
    waveForms=ExtractChunks(ephysData.traces(wfNum,:),...
            spikesTimes,100,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
    % scale to resolution
    waveForms=waveForms.*ephysData.recInfo.bitResolution;
    subplot(2,16,wfNum)
    plot(mean(waveForms))
end


figure
for wfNum=1:size(phyData.templates,3)
    waveForm=squeeze(phyData.templates(10,:,wfNum));
    subplot(2,16,wfNum) %double(phyData.channel_map(
    plot(waveForm);
    ylim([-0.6 0.2])
end