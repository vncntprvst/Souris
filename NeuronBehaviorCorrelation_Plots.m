% NeuronBehaviorCorrelation_Plots

%% First, get data
dirFiles=dir;
processedData=cellfun(@(x) contains(x,'processedData'), {dirFiles.name});
if sum(processedData)
    load(fullfile(dirFiles(processedData).folder,dirFiles(processedData).name));
else
    [ephys,behav]=NeuronBehaviorCorrelation_GatherData;
    save('processedData','ephys','behav','-v7.3');
end

% for whiskerTraceNum=1:size(whiskerTrackingData,1)
% end

bestUnits=EphysFunctions.FindBestUnits(ephys.spikes.unitID);
spikeRasters_ms=MakeRasters(spikeTimes,unitID,samplingRate,traceLength)
SDFs_ms=MakeSDF(spikeRasters_ms,unitNum)
    % make sure behavior and spike traces have same length
    if size(SDFs_ms,2)~=numel(periodBehavData_ms(whiskerTraceNum,:))
        % check what
    end
    
bestUnits=find(sum(ephys.spikeRasters_ms,2)/size(ephys.spikeRasters_ms,2)*10^5>100);
buSpikeRasters_ms=ephys.spikeRasters_ms(bestUnits,:);
buSDFs_ms=ephys.SDFs_ms(bestUnits,:);
keepTraces=unique(ephys.preferredElectrode(ismember(ephys.unitID,bestUnits)));
% Find whisking periods
peakWhiskingIdx=find(behav.peakWhisking_ms(1,:)==max(behav.peakWhisking_ms(1,:)));
% whiskingPeriod=peakWhiskingIdx-5000:peakWhiskingIdx+4999; %in ms
whiskingPeriodIdx=[0 abs(diff(behav.BP_periodBehavData_ms(1,:)))>=3*mad(diff(behav.BP_periodBehavData_ms(1,:)))] &...
                    abs(behav.BP_periodBehavData_ms(1,:))>=3*mad(behav.BP_periodBehavData_ms(1,:));
whiskingPeriodIdx=movsum(whiskingPeriodIdx,500)>0;
whiskingPeriodIdxInfo=bwconncomp(whiskingPeriodIdx);
whiskingPeriodDuration=cellfun(@(x) numel(x),whiskingPeriodIdxInfo.PixelIdxList);
longestWhiskingPeriodIdx=whiskingPeriodIdxInfo.PixelIdxList{...
    whiskingPeriodDuration==max(whiskingPeriodDuration)};

%% whisking vs spikes plot
NBC_Plots_SpikesWhisking(behav.BP_periodBehavData_ms,behav.whiskingPhase_ms,buSpikeRasters_ms);

%% 
NBC_Plots_WhiskingPhaseVideoFrame(behav.periodBehavData_ms,behav.vidTimes_ms,...
    behav.whiskingPhase_ms,buSpikeRasters_ms,buSDFs_ms);

%% overlay whisking and spikes on video 
NBC_Plots_SpikesWhiskingOverlayOnVideo(behav.periodBehavData_ms,behav.vidTimes_ms,...
    buSpikeRasters_ms,buSDFs_ms); %cursor_1,cursor_2

% Definitions

%% mock up units to test code
% mock up protraction unit
% fakeUnitSpikes=zeros(1,numel(behav.periodBehavData_ms));
% fakeUnitSpikes([0 diff(behav.BP_periodBehavData_ms)]>0 &...
% [0,0,diff(diff(behav.BP_periodBehavData_ms))]<0 &...
% behav.whiskingPhase_ms<0)=1;
% % mock up retraction unit
% fakeUnitSpikes=zeros(1,numel(behav.periodBehavData_ms));
% fakeUnitSpikes([0 diff(behav.BP_periodBehavData_ms)]<0 &...
% [0,0,diff(diff(behav.BP_periodBehavData_ms))]>0 &...
% behav.whiskingPhase_ms>0)=1;
% unitSpikes=fakeUnitSpikes; unitSDF=GaussConv(unitSpikes,5)*1000;

%% Is it tuned to Phase? - Polar plot summary 
% bestWhiskingPeriod=5; %if based on previous time chunk plots
% timeIndex=(bestWhiskingPeriod-1)*60+1;
% whiskingPeriod=int32(linspace(timeIndex*1000,(timeIndex+60)*1000-1,60*1000));
% whiskingPeriod=true(1,length(behav.whiskingPhase_ms)); 
whiskingPeriod=whiskingPeriodIdx; %based on defined whisking periods
phaseTuningSummaryFig=figure('position',[969    49   944   948],'name',ephys.recName);
colormapSeed=lines;
cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];

for unitNum=1:size(buSpikeRasters_ms,1) %find(keepUnits==15) 
    unitSpikes=buSpikeRasters_ms(unitNum,:);
    unitSpikesExcerpt=unitSpikes(whiskingPeriod);
%     vibrissaAngle=behav.BP_periodBehavData_ms;   
%     vibrissaAngleExcerpt=vibrissaAngle(whiskingPeriod);
    % coordConversion=90; %adjust depending on camera position
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt+coordConversion; % *180/pi;
    % vibrissaAngleExcerpt=vibrissaAngleExcerpt/180*pi; %convert back to radians
    
    % Hilbert transform NEEDS ANGLE TO BE ZERO CENTERED !!!
%     HTvibrissaAngleExcerpt=HTBP_periodBehavData_ms(whiskingPeriod);
    % figure; plot(whiskingPhase);
%     whiskingPhaseExcerpt_ms=angle(HTvibrissaAngleExcerpt);
    whiskingPhaseExcerpt_ms=behav.whiskingPhase_ms(1,whiskingPeriod);
    spikeOnWPhase=whiskingPhaseExcerpt_ms(logical(unitSpikesExcerpt));
    spikeOnWPhase=spikeOnWPhase+pi;
%     figure;
    subplot(ceil(size(buSpikeRasters_ms,1)/4),4,unitNum)
    pph=polarhistogram(spikeOnWPhase,20,'Displaystyle','stairs',...
        'Normalization','count','LineWidth',2);
    %mean angle value
%     rad2deg(circ_mean((spikeOnWPhase')))+360
    phaseStats=circ_stats(spikeOnWPhase');
    if phaseStats.skewness>0.02 & circ_rtest(spikeOnWPhase')<0.05
        phaseTuning=rad2deg(phaseStats.mean);
           pph.DisplayStyle='bar'; pph.EdgeAlpha=0
        if phaseTuning>0 & phaseTuning<180
           pph.EdgeColor='r';pph.FaceColor='r';
        elseif phaseTuning>180 & phaseTuning~=0
           pph.EdgeColor='g';pph.FaceColor='g';
        else
           pph.EdgeColor='k';pph.FaceColor='k';
        end
    else
        phaseTuning=[];
    end
end
savefig(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.fig'])
saveas(phaseTuningSummaryFig,[ephys.recName '_phaseTuningSummary.png'])
% unitsOfInterestIdx=[5,1,6,10];unitsOfInterest=keepUnits(unitsOfInterestIdx);

%% Is WR bursting reliable? - Plot spike times with whisking angle
% only consider time periods when whisking occurs angle and phase
timeVector=1:numel(behav.periodBehavData_ms(1,:)); timeVector(~whiskingPeriodIdx)=NaN;
% plot(timeVector,behav.periodBehavData_ms-median(behav.periodBehavData_ms))
oscillationPattern=cos(behav.whiskingPhase_ms(1,:));
% cross correlation
SpikeWAngleCorrFigure=figure('position',[969    49   944   948],'name',ephys.recName);
for unitNum=1:size(buSDFs_ms,1) %find(keepUnits==15); %bestUnit=2; %4;
    unitSDF=buSDFs_ms(unitNum,:); % unitSpikes=spikeRasters_ms(unitNum,:);
    [acor,lag] = xcorr(unitSDF(whiskingPeriodIdx),...
        oscillationPattern(whiskingPeriodIdx),150,'coeff');
    figure(SpikeWAngleCorrFigure)
    subplot(ceil(size(buSpikeRasters_ms,1)/4),4,unitNum);
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
%         plot(timeVector,behav.BP_periodBehavData_ms(1,:)); %-median(behav.BP_periodBehavData_ms));
%         plot(timeVector,behav.whiskingPhase_ms(1,:))
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
%         title([ephys.recName ' unit ' num2str(bestUnits(unitNum))  ' Tuning ' tuningLabel],'interpreter','latex')
%         savefig(gcf,[ephys.recName '_BurstUnit' bestUnits(unitNum) '.fig'])
%         saveas(gcf,[ephys.recName '_BurstUnit' bestUnits(unitNum) '.png'])
%         close(gcf);
    end
end
savefig(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.fig'])
saveas(SpikeWAngleCorrFigure,[ephys.recName '_SpikeWAngleCorr.png'])

