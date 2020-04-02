function phaseTuning=NBC_Plots_PhaseTuning(whiskerAngle,whiskerPhase,ephysData,dataMask,splitEpochs,saveFig)
% whiskingPhase in radians
% spikeRate in Hz
% Kyle's plot adapted from UnitExplorer GetTuning

colormapSeed=lines; cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(copper))/2;autumn];

if isfield(ephysData.recInfo,'sessionName'); recName=ephysData.recInfo.sessionName;
else; recName='PhaseTuning_polarPlot'; end

%% Data masking: global
% whiskerPhase = whiskerPhase(dataMask);
% whiskerAngle = whiskerAngle(dataMask);
% spikeRasters = ephysData.rasters(ephysData.selectedUnits,dataMask);
% spikeRate=ephysData.spikeRate(ephysData.selectedUnits,dataMask);
%% Data masking: look at each whisking epoch
spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
spikeRate=ephysData.spikeRate(ephysData.selectedUnits,:);
wEpochs=bwconncomp(dataMask);
durationThd=cellfun(@(x) length(x),wEpochs.PixelIdxList)>4000;
dataMask(vertcat(wEpochs.PixelIdxList{~durationThd}))=false;
wEpochs.PixelIdxList=wEpochs.PixelIdxList(durationThd);
if splitEpochs
wEpochs.NumObjects=sum(durationThd);
else
    wEpochs.PixelIdxList={vertcat(wEpochs.PixelIdxList{:})};
    wEpochs.NumObjects=1;
end
% [~,sortDuration]=sort(cellfun(@(x) length(x),wEpochs.PixelIdxList),'descend');
% wEpochs.PixelIdxList=wEpochs.PixelIdxList(sortDuration);

%% normalized whisker angle
if false
    cycleAngle=cell(numel(ptWhisks.PixelIdxList),1);
    if numel(ptWhisks.PixelIdxList)>1000
        %take sample
        whiskList=1:round(numel(ptWhisks.PixelIdxList)/1000):numel(ptWhisks.PixelIdxList);
    else
        whiskList=1:numel(ptWhisks.PixelIdxList);
    end
    for whiskCycle=1:numel(whiskList)
        peakWhiskIdx=ptWhisks.PixelIdxList{whiskList(whiskCycle)}(1);
        %     if whiskCycle>1
        wCycleSIdx=peakWhiskIdx-find(whiskerPhase(peakWhiskIdx-1:-1:1)>0,1)+1;
        %     else
        %         wCycleSIdx=peakWhiskIdx-find(whiskerPhase(peakWhiskIdx-1:-1:1)<0,1,'last'); %should be one
        %     end
        if ~isempty(wCycleSIdx)
            wCycleEIdx=ptWhisks.PixelIdxList{whiskList(whiskCycle)}(end);
            [cycleAngle{whiskCycle},phaseValues]=...
                resample(whiskerAngle(wCycleSIdx:wCycleEIdx),...
                whiskerPhase(wCycleSIdx:wCycleEIdx));
            cycleAngle{whiskCycle} = interp1(phaseValues,cycleAngle{whiskCycle},linspace(-pi,pi,100));
            %normalize
            cycleAngle{whiskCycle} =rescale(cycleAngle{whiskCycle});
        else
            continue
        end
    end
    cycleAngle=vertcat(cycleAngle{:});
%             plot(linspace(-pi,pi,100),nanmean(cycleAngle),'k')
%             set(gca,'ytick',[0 0.5 1],'xtick',[-pi 0 pi],'xticklabel',{'0','\pi','2\pi'})
%             set(gca,'tickdir','out')
%             % xlabel('Phase \phi (radians)')
%             ylabel('Normalized whisker angle')
end


for unitNum=1:size(spikeRasters,1)
    numWepochs=wEpochs.NumObjects;
    if wEpochs.NumObjects>1
        figure('Color','white','position',[165         153        1143         911])
    else
        figure('Color','white','position',[1278         200         634         711]);
    end
    subplot(4,numWepochs,1:numWepochs); hold on; 
    plot(whiskerAngle,'k'); 
%     if wEpochs.NumObjects>1
        plot(dataMask*nanstd(whiskerAngle)+nanmean(whiskerAngle),'r','linewidth',1.5)
%     end
    set(gca,'tickdir','out'); axis tight
    title(['Unit ' num2str(ephysData.selectedUnits(unitNum))]);
    for wEpochNum=1:numWepochs
        clearvars eWhiskerAngle eWhiskerPhase sp2H sp3H sp4H
        eWhiskerAngle=whiskerAngle(wEpochs.PixelIdxList{wEpochNum});
        eWhiskerPhase=whiskerPhase(wEpochs.PixelIdxList{wEpochNum});
        
%         ptWhisks=bwconncomp(eWhiskerPhase>0);
        %% probability density function of phase for spiking events
        numBins=32; % each bin = pi/16 radians
        edges = linspace(min(eWhiskerPhase), max(eWhiskerPhase), numBins*2+1);
%         edges = linspace(-pi-pi/numBins,pi+pi/numBins, numBins+1);
        centers = mean([ edges(1:end-1); edges(2:end) ]);
        [ ~, ~, phaseBins ] = histcounts(eWhiskerPhase, edges);
        samplingRate=1000; %change in case this isn't at 1kHz SR
        phaseTuning=nan(size(spikeRate,1),numWepochs);

        try
            unitSpikeEvent=spikeRasters(unitNum,wEpochs.PixelIdxList{wEpochNum});
            attribPhaseBin = phaseBins(logical(unitSpikeEvent));
            %     number of spikes in each phase bin N(?k|spike)
            spikePhaseBinCount=histcounts(attribPhaseBin,[1 1+unique(phaseBins)]);%         [spikePhaseBinCount,uniqueSpikePhaseBins]=hist(phaseVals,unique(phaseVals));
            %     probability density function P(?k|spike)
            spikePhasePDF=spikePhaseBinCount/sum(spikePhaseBinCount);
            spikePhasePDF=sum(reshape([spikePhasePDF(2:end),spikePhasePDF(1)],2,numBins))/2;
            spikePhasePDF=[spikePhasePDF(end) spikePhasePDF];
            %         spikePhasePDF=movsum(spikePhasePDF,6);spikePhasePDF=spikePhasePDF(6:6:end);
            %     number of phase occurence for each phase bin N(?k)
            phaseBinCount=histcounts(phaseBins,[1 1+unique(phaseBins)]);
            %     probability density function P(?k)
            phasePDF=phaseBinCount/sum(phaseBinCount);
            phasePDF=sum(reshape([phasePDF(2:end),phasePDF(1)],2,numBins))/2;
            phasePDF=[phasePDF(end) phasePDF];
            %         phasePDF=movsum(phasePDF,2);phasePDF=phasePDF(2:2:end);
            % mean spike rate for each phase bin ?[?k] = SR*N(?k|spike)/N(?k)
            meanPhaseSpikeRate=samplingRate*spikePhaseBinCount./phaseBinCount;
            %     fit sine wave
            %     modulation depth of the averaged whisking response
            
            %% plot
            % plot epoch angle trace
            if wEpochs.NumObjects>1
                subplot(4,numWepochs,1:numWepochs);
                plot(wEpochs.PixelIdxList{wEpochNum},eWhiskerAngle,'b')
            end

            %plot PDF
            sp2H=subplot(4,numWepochs,numWepochs+wEpochNum); hold on
            plot(linspace(-pi,pi, numBins+1),spikePhasePDF,'linewidth',1.2,'Color', [0 0 0]); %centers
            plot(linspace(-pi,pi, numBins+1),phasePDF,'linewidth',1.2,'Color', [0 0 0 0.5]); %centers
            set(gca,'ytick',0:0.05:1,...
                'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'0','\pi','2\pi'},...
                'tickdir','out'); 
            axis tight
            if wEpochs.NumObjects==1
%                 'ylim',[0 0.1]
              ylabel('Probability density')
            end
            
            %plot spike rate
%             subplot(numEpochs,5,(wEpochNum-1)*5+3); hold on
%             %bar / histogram of mean spike rates
%             bar(centers,meanPhaseSpikeRate)
%             % sine fit
%             yLim=get(gca,'ylim');
%             set(gca,'ytick',0:50:max(yLim),'xtick',[-pi 0 pi],'xticklabel',{'0','\pi','2\pi'},'tickdir','out')
% %             ylabel('Spike rate (Hz)')
% %             xlabel('Phase \phi (radians) ')
%             % end
            
            % for unitNum=1:size(spikeRate,1)
            clearvars unitSpikeRate
            eSpikeRate= spikeRate(unitNum,wEpochs.PixelIdxList{wEpochNum});
            [binMeanSpikeRate,binSESpikeRate]=deal(nan(numBins*2,1));
%             chunkLength=2000;
%             numChunks=floor(numel(unitSpikeRate)/chunkLength);
            % Kyle's plot
%             figure('Color','w');  box off
            sp3H=subplot(4,numWepochs,numWepochs*2+wEpochNum);hold on;
%             for chunkNum=4:numChunks
                clearvars binMeanSpikeRate binSESpikeRate
%                 chunkIndex=(chunkNum-1)*chunkLength+1:chunkNum*chunkLength;
                % Bining firing rate and spikes
                for binNum = 1:length(edges)-1 % : -1 : 1
%                     chunkSpikeRate=unitSpikeRate(chunkIndex);
                    ratesVect= eSpikeRate(phaseBins == binNum);%(chunkIndex)
                    numSample = numel(ratesVect);
                    if numSample == 0
                        meanSpikeRate = 0;
                        steSpikeRate = 0;
                    else
                        meanSpikeRate = nanmean(ratesVect);
                        steSpikeRate = MMath.StandardError(ratesVect);
                    end
                    binMeanSpikeRate(binNum) = meanSpikeRate;
                    binSESpikeRate(binNum) = steSpikeRate;
                end
                
            binMeanSpikeRate=sum(reshape([binMeanSpikeRate(2:end),binMeanSpikeRate(1)],2,numBins))/2;
            binMeanSpikeRate=[binMeanSpikeRate(end) binMeanSpikeRate];               
                
                %             shadedErrorBar(centers, binMeanSpikeRate,binSESpikeRate, 'lineprops','k');
            plot(linspace(-pi,pi, numBins+1), binMeanSpikeRate, 'LineWidth',2) %centers %,'color',cmap(unitNum,:));%'k' 
%             end
%             xlabel({'Phase  (rad)'; '0 = Max Protraction'}, 'FontSize', 18);
%             ylabel('Firing rate (Spk/s)', 'FontSize', 18);
            axis tight
            yl = ylim;
            ylim([0 yl(2)]);
            set(gca,'TickDir','out')
            box off
            xlim([-pi pi])
            if wEpochs.NumObjects==1
                xlabel('Phase (radians)')
                ylabel('Spike rate (Hz)')
            end
            
            % set(gca,'xdir', 'reverse'); %, 'ydir', 'reverse')
            %% convert to thetas: make as many phase # as FR for that phase #
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
                phEdgeColor=cmap(unitNum,:);phFaceColor=cmap(unitNum,:);
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
            
            %     if size(spikeRate,1)>1
            %         subplot(ceil(size(spikeRate,1)/4),4,unitNum)
            %     end
            
            sp4H=subplot(4,numWepochs,numWepochs*3+wEpochNum); %subplot(2,1,2);
            polarhistogram(thetas,binNum,'Displaystyle','bar',...
                'Normalization','count','LineWidth',2,...
                'EdgeColor',phEdgeColor,'FaceColor',phFaceColor,...
                'EdgeAlpha',0);
            paH = gca;
            paH.ThetaZeroLocation='left';
            paH.ThetaTickLabel={'maxP','','','R','','',...
                'maxR','','','P','',''};
            paH.ThetaDir = 'counterclockwise';
            % For Kyle's convention:
            %     paH.ThetaDir = 'clockwise';

        catch
            continue
        end
    end
%                 spPos=get(sp2H,'position');
%             spPos(1)=0.1+((wEpochNum-1)*(1/numWepochs));
%             spPos(3)=1/numWepochs;
%             set(sp2H,'position',spPos);
% 
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
