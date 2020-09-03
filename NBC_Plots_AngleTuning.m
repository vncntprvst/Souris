function angleTuning=NBC_Plots_AngleTuning(whiskerAngle,ephysData,dataMask,...
    labels,splitEpochs,saveFig)
% whiskingAngle in radians
% spikeRate in Hz

colormapSeed=lines; cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(copper))/2;autumn];

if isfield(ephysData.recInfo,'sessionName'); recName=ephysData.recInfo.sessionName;
else; recName='AngleTuning_polarPlot'; end

% find base setpoint
baseSP=mode(round(whiskerAngle/10)*10);
if baseSP<0
    whiskerAngle=whiskerAngle+180;
end
baseSP=deg2rad(baseSP);

%% Data masking: global
% whiskerAngle = whiskerAngle(dataMask);
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

for unitNum=1:size(spikeRasters,1)
    numWepochs=wEpochs.NumObjects;
    if wEpochs.NumObjects>1
        figure('Color','white','position',[165         153        1143         911])
    else
        figure('Color','white','position',[1035         224         802         672]);
    end
    %     subplot(4,numWepochs,1:numWepochs); hold on;
    %     plot(whiskerAngle,'k');
    %     if wEpochs.NumObjects>1
    %         plot(dataMask*nanstd(whiskerAngle)+nanmean(whiskerAngle),'r','linewidth',1.5)
    %     end
    %     set(gca,'tickdir','out'); axis tight
    %     title({['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' recName];
    %         ['Tuning to ' labels ' phase']},'interpreter','none');
    for wEpochNum=1:numWepochs
        clearvars eWhiskerAngle sp2H sp3H sp4H
        eWhiskerAngle=whiskerAngle(wEpochs.PixelIdxList{wEpochNum});
        eWhiskerAngle=deg2rad(eWhiskerAngle);
        %% probability density function of phase for spiking events
        numBins=16; %
        edges = linspace(min(eWhiskerAngle), max(eWhiskerAngle), numBins*2+1);
        centers = mean([ edges(1:end-1); edges(2:end) ]);
        [ ~, ~, angleBins ] = histcounts(eWhiskerAngle, edges);
        samplingRate=1000; %change in case this isn't at 1kHz SR
        angleTuning=nan(size(spikeRate,1),numWepochs);
        
        try
            unitSpikeEvent=spikeRasters(unitNum,wEpochs.PixelIdxList{wEpochNum});
            attribAngleBin = angleBins(logical(unitSpikeEvent));
            %     number of spikes in each angle bin N(theta_k|spike)
            spikeAngleBinCount=histcounts(attribAngleBin,1:max(unique(angleBins))+1);%         [spikeAngleBinCount,uniqueSpikeAngleBins]=hist(phaseVals,unique(phaseVals));
            %     probability density function P(theta_k|spike)
            spikeAnglePDF=spikeAngleBinCount/sum(spikeAngleBinCount);
            spikeAnglePDF=sum(reshape([spikeAnglePDF(2:end),spikeAnglePDF(1)],2,numBins))/2;
            spikeAnglePDF=[spikeAnglePDF(end) spikeAnglePDF];
            %         spikeAnglePDF=movsum(spikeAnglePDF,6);spikeAnglePDF=spikeAnglePDF(6:6:end);
            %     number of phase occurence for each phase bin N(?k)
            angleBinCount=histcounts(angleBins,1:max(unique(angleBins))+1);
            %     probability density function P(theta_k)
            anglePDF=angleBinCount/sum(angleBinCount);
            anglePDF=sum(reshape([anglePDF(2:end),anglePDF(1)],2,numBins))/2;
            anglePDF=[anglePDF(end) anglePDF];
            %         phasePDF=movsum(phasePDF,2);phasePDF=phasePDF(2:2:end);
            % mean spike rate for each angle bin ?[?k] = SR*N(theta_k|spike)/N(theta_k)
            meanAngleSpikeRate=samplingRate*spikeAngleBinCount./angleBinCount;
            %     fit sine wave
            %     modulation depth of the averaged whisking response
            
            %% plots
            
            %plot normalized angle traces
            %             subplot(4,numWepochs,wEpochNum); hold on %horizontal arrangment:subplot(numEpochs,4,(wEpochNum-1)*4+1); %vertical arrangment: subplot(4,numWepochs,wEpochNum);
            %             subplot(4,numWepochs,1:numWepochs);
            %             plot(wEpochs.PixelIdxList{wEpochNum},eWhiskerAngle,'b')
            %             set(gca,'ytick',[0 0.5 1],'xtick',[-pi 0 pi],'xticklabel',{'0','\pi','2\pi'})
            %% plot PDF
            sp2H=subplot(2,2,4); hold on ; %(4,numWepochs,numWepochs+wEpochNum)
            plot(linspace(edges(1),edges(end), numBins+1),spikeAnglePDF,'linewidth',1.2,'Color', [0 0 0]); %centers
            plot(linspace(edges(1),edges(end), numBins+1),anglePDF,'linewidth',1.2,'Color', [0 0 0 0.5]); %centers
            xlabels={num2str(round(baseSP,2)),'0','\pi/2','\pi'};
            [sortLabels,sortLabelsIdx]=sort([baseSP 0 pi/2 pi]);
            set(gca,'ytick',0:0.05:1,...
                'xlim',[0 pi],'xtick',sortLabels,'xticklabel',xlabels(sortLabelsIdx),...
                'tickdir','out');
            
            %             axis tight
            if wEpochs.NumObjects==1
                %                 'ylim',[0 0.1]
                %               ylabel('PDF')
            end
            legend('P(\theta_k|spike)','P(\theta_k)','location','southeast')
            legend('boxoff')
            title({'Probability density function'; 'of angle for spiking events'})
            
            %plot spike rate
            %             subplot(numEpochs,5,(wEpochNum-1)*5+3); hold on
            %             %bar / histogram of mean spike rates
            %             bar(centers,meanAngleSpikeRate)
            %             % sine fit
            %             yLim=get(gca,'ylim');
            %             set(gca,'ytick',0:50:max(yLim),'xtick',[-pi 0 pi],'xticklabel',{'0','\pi','2\pi'},'tickdir','out')
            % %             ylabel('Spike rate (Hz)')
            % %             xlabel('Angle \phi (radians) ')
            %             % end
            
            % for unitNum=1:size(spikeRate,1)
            
            %% Plot spike rate
            clearvars unitSpikeRate binMeanSpikeRate binSESpikeRate
            eSpikeRate= spikeRate(unitNum,wEpochs.PixelIdxList{wEpochNum});
            % Bining firing rate and spikes
            [binMeanSpikeRate,binSESpikeRate]=deal(nan(numBins*2,1));
            for binNum = 1:numBins*2 % : -1 : 1
                %                     chunkSpikeRate=unitSpikeRate(chunkIndex);
                ratesVect= eSpikeRate(angleBins == binNum);%(chunkIndex)
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
            
            sp3H=subplot(2,2,3);hold on; %(4,numWepochs,numWepochs*2+wEpochNum)
            plot(linspace(edges(1),edges(end), numBins*2), binMeanSpikeRate, 'LineWidth',2) %centers %,'color',cmap(unitNum,:));%'k'
            axis tight
            yl = ylim;
            ylim([0 yl(2)]);
            text(baseSP,  yl(2)/10, 'resting state','HorizontalAlignment','Center')
            set(gca,'xlim',[0 pi],'xtick',sortLabels,'xticklabel',xlabels(sortLabelsIdx),...
                'tickdir','out');
            box off
            if wEpochs.NumObjects==1
                xlabel('Angle (radians)')
                ylabel('Spike rate (Hz)')
            end
            title({'Average spike rate'; ['across ' labels ' angles']})
            % set(gca,'xdir', 'reverse'); %, 'ydir', 'reverse')
            
            %% Polar plot
            % convert to thetas: make as many angle # as FR for that angle #
            thetas=cell(numel(centers),1);
            for binNum=1:numel(centers)
%                 round(binMeanSpikeRate(binNum))
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
                angleTuning(unitNum,1)=rad2deg(phaseStats.mean);
            end
            if ~isnan(angleTuning(unitNum,1))
                phEdgeColor=cmap(unitNum,:);phFaceColor=cmap(unitNum,:);
            else
                phEdgeColor='k';phFaceColor='k'; %EdgeAlpha=0.5;
            end
            
            sp4H=subplot(2,2,[1 2]);
            polarhistogram(thetas,edges,'Displaystyle','bar',...
                'Normalization','count','LineWidth',2,...
                'EdgeColor',phEdgeColor,'FaceColor',phFaceColor,...
                'EdgeAlpha',0);%binNum
            paH = gca;
            paH.ThetaLim = [0 180];
            paH.ThetaZeroLocation='right';
            paH.ThetaTickLabel={'max Retraction','','','','','',...
                'max Protraction'};
            paH.ThetaDir = 'counterclockwise';
            paH.RAxis.Label.String='Mean spike rate';
            
            title(['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' strrep(recName,'_','') ...
                ' - Tuning to ' labels ' angle'],'interpreter','none');
            
        catch
            continue
        end
    end
end

if saveFig
    savefig(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.fig'])
    saveas(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.png'])
end
