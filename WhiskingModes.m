    %Looking at whisking mode separation
    
    

    % first def of peak Idx
    figure; hold on 
    plot(whiskerAmplitude); plot(whiskerAngle)
    peakIdx = whiskerPhase(1:end-1)<0 & whiskerPhase(2:end)>=0;
    plot(find(peakIdx),zeros(sum(peakIdx),1),'kd')
    wfreq=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);
    plot(wfreq)
    fAmp=WhiskingFun.GetAmplitude(whiskerAngle,whiskerPhase);
    
    % second def of peak Idx, based on velocity 
% Because velocity is highest during whisking, however, it is more common
% the estimate the whisking frequency from the angular velocity (??) which 
% contains a single harmonic typically denoted f0 (Berg and Kleinfeld 2003).
    figure; hold on
    plot(whiskerVelocity)
    figure; hold on
    histogram(abs(whiskerVelocity))

    whisksIdx = bwconncomp(whiskerVelocity<0);
    peakIdx = zeros(length(whiskerVelocity),1);
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
    plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    wfreq=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);
    figure; hold on
    plot(wfreq)
    
    figure; hold on 
    plot(whiskerAngle)
    plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    plot(fAmp)
    
    figure; hold on
    histogram(wfreq(whiskerAmplitude<=5))
    histogram(wfreq(whiskerAmplitude>5 & whiskerAmplitude<15))
    histogram(wfreq(whiskerAmplitude>=15))
    legend('low amp','mid amp','high amp')
    title('Frequency distribution for different whisking amplitudes')
        
    figure; hold on
    histogram(whiskerAmplitude(wfreq<=5))
    histogram(whiskerAmplitude(wfreq>5 & wfreq<15))
    histogram(whiskerAmplitude(wfreq>=15))
    legend('low freq','mid freq','high freq')
    title('Amplitude distribution for different whisking frequencies')
    
    figure; hold on
    histogram(whiskerAmplitude(whiskerSetPoint<70))
    histogram(whiskerAmplitude(whiskerSetPoint>=70))
    legend('low setpoint','high setpoint')
    title('Amplitude distribution for different set-points')
    
    figure; hold on
    histogram(wfreq(whiskerSetPoint<70))
    histogram(wfreq(whiskerSetPoint>=70))
    legend('low setpoint','high setpoint')
    title('Frequency distribution for different set-points')
    
    figure; hold on
    histogram(whiskerSetPoint(whiskerAmplitude<=5))
    histogram(whiskerSetPoint(whiskerAmplitude>5 & whiskerAmplitude<15))
    histogram(whiskerSetPoint(whiskerAmplitude>=15))
    legend('low amp','mid amp','high amp')
    title('SetPoint distribution for different whisking amplitudes')

    figure; hold on
    histogram(whiskerAmplitude(whiskerSetPoint<70 & whiskerAmplitude<15))
    histogram(whiskerAmplitude(whiskerSetPoint<70 & whiskerAmplitude>=15))
    histogram(whiskerAmplitude(whiskerSetPoint>70 & whiskerAmplitude>=15))
    legend('low SP low amp','low SP high amp','high SP high amp')
    title('Amplitude distribution for different whisking modes')

    
    %% plot angle and velocity by whisk duration
    % Figure from Knutsen's Whisking Kinematics
%     https://link.springer.com/chapter/10.2991/978-94-6239-133-8_46

    whiskerSetPoint=behav.whiskerTrackingData.SetPoint;
    whiskerAngle_nf=behav.whiskerTrackingData.Angle;
    whiskerAngle_nf=fillmissing(whiskerAngle_nf,'linear');
    whiskerAngle_LP=WhiskingFun.LowPassBehavData(whiskerAngle_nf,1000,40);
    whiskerVelocity_LP=diff(whiskerAngle_LP); whiskerVelocity=[whiskerVelocity(1) whiskerVelocity]; 
    figure; hold on
    plot(whiskerAngle_LP)
    plot(whiskerSetPoint)

    figure; hold on; plot(whiskerPhase); plot(diff(whiskerPhase))
    
    whisksIdx = bwconncomp(diff(whiskerPhase)>0); %troughIdx = find(whiskerPhase(1:end-1)>=pi/2 & whiskerPhase(2:end)<=-pi/2);
    whiskDur=cellfun(@(whisk) length(whisk), whisksIdx.PixelIdxList);
    whisksIdx.PixelIdxList=whisksIdx.PixelIdxList(whiskDur>1); %prune
    whiskDur=cellfun(@(whisk) length(whisk), whisksIdx.PixelIdxList);

    [durBinCount,binEdges,binIdx]=histcounts(whiskDur,40);   
    binCat=unique(binIdx);
    [allCatMeanAngle,allCatMeanVel]=deal(nan(numel(binCat),100));
    for whiskDurCat=1:numel(binCat)
        catWhiskList=whisksIdx.PixelIdxList(binIdx==binCat(whiskDurCat));
%         catWhiskList=catWhiskList(cellfun(@numel, catWhiskList)>1); %prune out
        catMeanAngle=cellfun(@(catWhisk) interp1(whiskerAngle_LP(catWhisk),...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanAngle(whiskDurCat,:)=median(vertcat(catMeanAngle{:}));
        catMeanVel=cellfun(@(catWhisk) interp1(whiskerVelocity_LP(catWhisk)*1000,...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanVel(whiskDurCat,:)=median(vertcat(catMeanVel{:}));
    end
    binlabels=round(binEdges+mode(diff(binEdges)/2));
    figure; 
    subplot(1,2,1)
    imagesc(allCatMeanAngle); colorbar
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    subplot(1,2,2)
    imagesc(allCatMeanVel); colorbar
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    
