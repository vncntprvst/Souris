%% Population / Summary plots
clearvars;
processedFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*processedData.mat'},'UniformOutput', false); %whiskerTrackingData
colormapSeed=lines;
cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
samplingRate=30000;
for fileNum=1:numel(processedFiles{:})
    clearvars -except processedFiles fileNum cmap samplingRate %behav ephys
    load(fullfile(processedFiles{1}(fileNum).folder,processedFiles{1}(fileNum).name));
    %% narrow number of units
    bestUnits=EphysFunctions.FindBestUnits(ephys.spikes.unitID);
    [spikeRasters_ms,unitList]=EphysFunctions.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
        ephys.samplingRate,int32(size(ephys.traces,2)/ephys.samplingRate*1000));
    keepUnits=ismember(unitList,bestUnits); %keepUnits=unitList % decide which units to keep
    keepTraces=nan(1,numel(bestUnits));
    for unitNum=1:numel(bestUnits)
        keepTraces(unitNum)=mode(ephys.spikes.preferredElectrode(...
            ephys.spikes.unitID==bestUnits(unitNum)));
    end
    %     keepTraces=unique(keepTraces);
    % SDFs_ms=EphysFunctions.MakeSDF(spikeRasters_ms); % compute spike density functions
    spikeRasters_ms=spikeRasters_ms(keepUnits,:);
    % SDFs_ms=SDFs_ms(keepUnits,:);
    %% whisking data resampling.
    periodBehavData_ms=WhiskingAnalysisFunctions.ResampleBehavData... % use ephys sampling rate, if video times are based on TTLs
        (behav.whiskerTrackingData,behav.vidTimes,ephys.samplingRate);
    if size(spikeRasters_ms,2)~=size(periodBehavData_ms,2) % make sure behavior and spike traces have same length
        if size(spikeRasters_ms,2)<size(periodBehavData_ms,2)
            periodBehavData_ms=periodBehavData_ms(:,1:size(spikeRasters_ms,2));
        else
            %             SDFs_ms=SDFs_ms(:,1:size(periodBehavData_ms,2));
            spikeRasters_ms=spikeRasters_ms(:,1:size(periodBehavData_ms,2));
        end
    end
    %% filtered whisking traces
    whiskerTraces_whiskFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(periodBehavData_ms,1000,[4 20]);
    whiskerTraces_breathFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(periodBehavData_ms,1000,[1 4]);
    %% find phase
    whiskingPhase_ms=WhiskingAnalysisFunctions.ComputePhase(whiskerTraces_whiskFreq_ms); %,whiskingEpochs);
    
    %% Find whisking periods of at least 500ms
    [whiskingEpochsIdx,wAmplitude,setPoint]=WhiskingAnalysisFunctions.FindWhiskingEpochs(...
        whiskerTraces_whiskFreq_ms(1,:),whiskingPhase_ms(1,:),500);
    whiskingEpochsInfo= bwconncomp(whiskingEpochsIdx) ;
    
    whiskBoutDuration=cellfun(@(whiskBout) numel(whiskBout), whiskingEpochsInfo.PixelIdxList);
    longestWhiskBout=whiskBoutDuration==max(whiskBoutDuration);
    longWhiskBoutInit=whiskingEpochsInfo.PixelIdxList{longestWhiskBout}(1);
    
    %% get average waveform
    if isempty(ephys.spikes.waveforms)
        continue
    end
    avWaveform=nan(numel(bestUnits), size(ephys.spikes.waveforms,2));
    %     figure; hold on
    for unitNum=1:numel(bestUnits)
        avWaveform(unitNum,:)=nanmean(ephys.spikes.waveforms(...
            ephys.spikes.unitID==bestUnits(unitNum),:));
        %         plot(avWaveform(unitNum,:))
    end
    % neuron with biggest waveform
    biggestUnit=sum(abs(diff(avWaveform,[],2)),2)==max(sum(abs(diff(avWaveform,[],2)),2));
    displayWindow=max([longWhiskBoutInit-1000 1]):min([longWhiskBoutInit+58999 numel(whiskingEpochsIdx)-1]) ;
    traces=PreProcData(ephys.traces(:,displayWindow*30),ephys.spikes.samplingRate,{'CAR','all'});
    bestTrace=traces(keepTraces(biggestUnit),:);
    
    % plot excerpts
    traceExcerptFigure=figure('position',[ 334 703 1264 272]); hold on
    plot(bestTrace,'k');
    plot(linspace(1,numel(bestTrace),numel(displayWindow)),...
        10*whiskerTraces_whiskFreq_ms(1,displayWindow)+max(bestTrace)+100);
    plot(linspace(1,numel(bestTrace),numel(displayWindow)),...
        10*whiskerTraces_breathFreq_ms(1,displayWindow)+max(bestTrace)+100);
    savefig(traceExcerptFigure,[ephys.recName(1:end-4) '_traceExcerpt.fig'])
    saveas(traceExcerptFigure,[ephys.recName(1:end-4) 'traceExcerpt.png'])
    close(traceExcerptFigure)
    % plot unit summary
    UnitSummaryFigure=figure('position',[469  428 1216 597]); hold on
    %% Plot average waveform
    for unitNum=1:length(bestUnits)
        subplot(4,length(bestUnits),unitNum)
        unitWaveforms=single(ephys.spikes.waveforms(...
            ephys.spikes.unitID==bestUnits(unitNum),:));
        if ~isnan(mean(unitWaveforms))
            lineh(unitNum)=plot(mean(unitWaveforms),'linewidth',2,'Color',...
                [cmap(bestUnits(unitNum),:),0.7]);
            wfSEM=std(unitWaveforms)/ sqrt(size(unitWaveforms,2)); %standard error of the mean
            wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
            patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
                [mean(unitWaveforms)-wfSEM,fliplr(mean(unitWaveforms)+wfSEM)],...
                cmap(bestUnits(unitNum),:),'EdgeColor','none','FaceAlpha',0.2);
        end
        set(gca,'XTick',linspace(0,size(unitWaveforms,1),5),...
            'XTickLabel',round(linspace(-round(size(unitWaveforms,1)/2),...
            round(size(unitWaveforms,1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
        axis('tight');box off;
        xlabel('Time (ms)');
        ylabel('Voltage (\muV)');
        set(gca,'Color','white','FontSize',10,'FontName','Calibri');
        hold off
        
        %% Plot autocorrelogram
        subplot(4,length(bestUnits),unitNum+length(bestUnits))
        %get unit spike times
        unitST=ephys.spikes.times(ephys.spikes.unitID==bestUnits(unitNum));
        % change to ms timescale
        unitST=unitST/(samplingRate/1000);
        unitST=unitST(unitST>0);
        spikeTimeIdx=zeros(1,unitST(end));
        spikeTimeIdx(unitST)=1;
        binSize=1;
        numBin=ceil(size(spikeTimeIdx,2)/binSize);
        binUnits = histcounts(double(unitST), linspace(0,size(spikeTimeIdx,2),numBin));
        binUnits(binUnits>1)=1; %no more than 1 spike per ms
        % compute autocorrelogram
        [ACG,lags]=xcorr(double(binUnits),200,'unbiased'); %'coeff'
        ACG(lags==0)=0;
        ACGh=bar(lags,ACG);
        ACGh.FaceColor = cmap(bestUnits(unitNum),:);
        ACGh.EdgeColor = 'none';
        % axis('tight');
        box off; grid('on'); %set(gca,'yscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
        xlabel('Autocorrelogram (1 ms bins)')
        set(gca,'xlim',[-20 20],... %'ylim',[0 max([max(get(gca,'ylim')) 10^1])]
            'Color','white','FontSize',10,'FontName','Calibri','TickDir','out');
        hold off
        
        %plot angle tuning
        subplot(4,length(bestUnits),unitNum+length(bestUnits)*2); hold on
        NBC_Plots_SpikingWhiskAngleTuning(whiskerTraces_whiskFreq_ms,whiskingPhase_ms,...
            whiskingEpochsIdx,spikeRasters_ms(unitNum,:),false,'subplot') %periodBehavData_ms
        %plot phase tuning
        subplot(4,length(bestUnits),unitNum+length(bestUnits)*3);
        NBC_Plots_PhaseTuning(whiskingPhase_ms,whiskingEpochsIdx,...
            spikeRasters_ms(unitNum,:),false,'subplot')
        
    end
    savefig(UnitSummaryFigure,[ephys.recName(1:end-4) '_UnitSummary.fig'])
    saveas(UnitSummaryFigure,[ephys.recName(1:end-4) '_UnitSummary.png'])
    close(UnitSummaryFigure)
    
    save([ephys.recName(1:end-4) '_procdata'],'whiskerTraces_whiskFreq_ms','whiskingPhase_ms',...
            'whiskingEpochsIdx','spikeRasters_ms','bestTrace','biggestUnit','whiskerTraces_breathFreq_ms')
end