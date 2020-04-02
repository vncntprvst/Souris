function NBC_Plots_Continuity(ephysData,selectedUnits,whiskerAngle,dataMask,saveFig)
% plot waveforms, ACG ISI across wisking epochs

colormapSeed=lines; cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(copper))/2;autumn];
if saveFig
    recName=regexp(cd,['(?<=\' filesep ')\w+$'],'match','once');
end

samplingRate=ephysData.samplingRate;
wEpochs=bwconncomp(dataMask);
durationThd=cellfun(@(x) length(x),wEpochs.PixelIdxList)>4000;
dataMask(vertcat(wEpochs.PixelIdxList{~durationThd}))=false;
wEpochs.PixelIdxList=wEpochs.PixelIdxList(durationThd);
wEpochs.NumObjects=sum(durationThd);

for unitNum=1:numel(selectedUnits)
    unitIDs=ephysData.unitID==selectedUnits(unitNum);
    unitSpikeTimes=double(ephysData.times(unitIDs))/samplingRate*1000;
    unitWaveforms=ephysData.waveforms(unitIDs,:);
    
    numWepochs=wEpochs.NumObjects;
    if wEpochs.NumObjects>1
        ContinuityFig=figure('Color','white','position',[72         205        1120         730]);
    else
        ContinuityFig=figure('Color','white','position',[1278         200         634         711]);
    end
    subplot(4,numWepochs,1:numWepochs); hold on; axis tight
    plot(whiskerAngle,'k'); 
    if wEpochs.NumObjects>1
        plot(dataMask*nanstd(whiskerAngle)+nanmean(whiskerAngle),'r','linewidth',1.5)
    end
    title(['Unit ' num2str(selectedUnits(unitNum))]);
    for wEpochNum=1:numWepochs
        %% look at each whisking epoch
        eIdx=wEpochs.PixelIdxList{wEpochNum};
        unitEIdx=unitSpikeTimes>=eIdx(1) & unitSpikeTimes<=eIdx(end);
        
        wavH=subplot(4,numWepochs,numWepochs+wEpochNum); hold on
        plot(wavH,mean(unitWaveforms(unitEIdx,:)),'linewidth',2,'color',cmap(selectedUnits(unitNum),:));
        
        
        wfSEM=std(unitWaveforms(unitEIdx,:))/ sqrt(size(unitWaveforms(unitEIdx,:),2)); %standard error of the mean
        wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(unitWaveforms(unitEIdx,:))-wfSEM,fliplr(mean(unitWaveforms(unitEIdx,:))+wfSEM)],...
            cmap(selectedUnits(unitNum),:),'EdgeColor','none','FaceAlpha',0.2); %cmap(cellNum,:)
        
        set(gca,'xcolor','none','ycolor','none');axis('tight');box off;
        %     xlabel('Time (ms)');%     ylabel('Voltage (\muV)');
        
        acgH=subplot(4,numWepochs,numWepochs*2+wEpochNum);
        try
            EphysFun.PlotACG(ones(1,sum(unitEIdx))*selectedUnits(unitNum),unitSpikeTimes(unitEIdx),selectedUnits(unitNum),1000,acgH,cmap)
        catch
        end
        if wEpochs.NumObjects>1
            set(gca,'xcolor','none','ycolor','none');
        end
        box off;
        
        isiH=subplot(4,numWepochs,numWepochs*3+wEpochNum);
        try
            EphysFun.PLotISI(ones(1,sum(unitEIdx))*selectedUnits(unitNum),unitSpikeTimes(unitEIdx),selectedUnits(unitNum),1000,isiH,cmap)
        catch
        end
        if wEpochs.NumObjects>1
            set(gca,'xcolor','none','ycolor','none');
        end
        box off;
    end

    if saveFig
        savefig(ContinuityFig,[recName '_Unit' num2str(selectedUnits(unitNum))  '_Continuity.fig'])
        saveas(ContinuityFig,[recName '_Unit' num2str(selectedUnits(unitNum)) '_Continuity.png'])
        close
    end
end
