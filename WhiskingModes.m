    %% Looking at whisking mode separation
    %% Four main modes found: 
    % foveal: high frequency > 10Hz, medium amplitude >25 <35, high setpoint/angular values >70 at start of whisk
    % exploratory: lower frequency < 10, high amplitude >35, medium setpoint/angular values ?
    % resting: lower frequency < 10Hz, low/medium amplitude <25, low setpoint/angular values <60
    % twiches: high frequency > 25Hz, low amplitude <10, low setpoint/angular values <70 at start of whisk
    
    
    %% Find session files and load
    whiskerFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*_wMeasurements*'},'UniformOutput', false);
    whiskerFiles=vertcat(whiskerFiles{~cellfun('isempty',whiskerFiles)});
    for fileNum=1:numel(whiskerFiles)
        wData(fileNum)=load(fullfile(whiskerFiles(fileNum).folder,whiskerFiles(fileNum).name));
    end
    
    wAngle=[wData(:).Angle];
    wAngle_nf=[wData(:).Angle_raw];
    wPhase=[wData(:).Phase];
    wAmplitude=[wData(:).Amplitude];
    wVelocity=[wData(:).Velocity];
    wSetPoint=[wData(:).SetPoint];
    wFrequency=[wData(:).Freq];
    samplingRate=wData(1).samplingRate;
    
    isNanIdx=isnan(wAngle) | isnan(wVelocity) | isnan(wAmplitude) | isnan(wFrequency) | isnan(wSetPoint);
    wAngle=wAngle(~isNanIdx); wVelocity=wVelocity(~isNanIdx); wAmplitude=wAmplitude(~isNanIdx); wFrequency=wFrequency(~isNanIdx); wSetPoint=wSetPoint(~isNanIdx);
    
    %% compute frequency (not inst. freq)
    whisksIdx = bwconncomp(wVelocity<0);
    peakIdx = zeros(length(wVelocity),1);
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
%     plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    wfreq=movsum(peakIdx,samplingRate);
    
    %% clustering by mean amp, freq ... values
    threshold = 4/(samplingRate/2);
    [coeffB,coeffA] = butter(3,threshold,'low');
    velFilter = filtfilt(coeffB, coeffA,abs(wVelocity)) ;
    velEpochsIdx = velFilter>0.2; % threshold by velocity
    
    whiskBoutList = bwconncomp(velEpochsIdx) ;
    
    wBoutVals=cellfun(@(wBout) [mean(wAmplitude(wBout)) mean(wfreq(wBout))...
        mean(wSetPoint(wBout)) length(wBout)] ,...
        whiskBoutList.PixelIdxList, 'UniformOutput' ,false);
    wBoutVals=vertcat(wBoutVals{:});
    wBoutVals=wBoutVals(wBoutVals(:,4)>100,:);
%     [~,wBoutValsSortIdx]=sort(wBoutVals(:,2));
%     wBoutVals=wBoutVals(wBoutValsSortIdx,:);
    labels=DimReducPlot(wBoutVals(:,1:3),1:size(wBoutVals,1),{'tsne','plot'});
    labelCat=unique(labels);
    for catNum=1:numel(labelCat)
        meanCatVals(catNum,:)=[catNum, round(mean(wBoutVals(labels==labelCat(catNum),:))), sum(labels==labelCat(catNum))];
    end    
    meanCatVals=table(meanCatVals(:,1),meanCatVals(:,2),meanCatVals(:,3),meanCatVals(:,4),meanCatVals(:,5),meanCatVals(:,6),...
        'VariableNames',{'LabelNum','wAmplitude','wFrequency','wSetPoint','Duration','Number'});

%         'wAmplitude','wFrequency','wSetPoint','Duration'
% with > 100 dur only / PCA
%       1     22          23          62          471     638   mid range in amp freq and dur
%       2     12          26          55          298     1804  short low amp high freq
%       3     36          19          87          7988    202   big amp high SP low freq
% with > 100 dur only / tsne
%       1	  10          27          53          282     879
%       2	  13          27          59          287     847
%       3	  25          22          66          2137    918

%     % with >500 dur only
%     cat1    18          16          57          528
%     cat2    19          15          57          725
%     cat3    30          15          76          6173
%     
%     
wBoutVals=cellfun(@(wBout) [min(wAmplitude(wBout)) max(wAmplitude(wBout))...
        min(wfreq(wBout))  max(wfreq(wBout))...
        min(wSetPoint(wBout)) max(wSetPoint(wBout)) length(wBout)] ,...
        whiskBoutList.PixelIdxList, 'UniformOutput' ,false);
    wBoutVals=vertcat(wBoutVals{:});
    wBoutVals=wBoutVals(wBoutVals(:,7)>100,:);
%     [~,wBoutValsSortIdx]=sort(wBoutVals(:,2));
%     wBoutVals=wBoutVals(wBoutValsSortIdx,:);
    labels=DimReducPlot(wBoutVals(:,1:6),1:size(wBoutVals,1),{'tsne','plot'});
    labelCat=unique(labels);
    for catNum=1:numel(labelCat)
        meanCatVals(catNum,:)=[catNum, round(mean(wBoutVals(labels==labelCat(catNum),:))), sum(labels==labelCat(catNum))];
    end    
    meanCatVals=table(meanCatVals(:,1),meanCatVals(:,2),meanCatVals(:,3),...
        meanCatVals(:,4),meanCatVals(:,5),meanCatVals(:,6),...
        meanCatVals(:,7),meanCatVals(:,8),meanCatVals(:,9),...
        'VariableNames',{'LabelNum','MinwAmplitude','MaxwAmplitude',...
        'MinwFrequency','MaxwFrequency','MinwSetPoint','MaxwSetPoint',...
        'Duration','Number'});
    
    %     'MinwAmplitude','MaxwAmplitude','MinwFrequency','MaxwFrequency','MinwSetPoint','MaxwSetPoint'
% short dur   8           15              25              27              60              62  twitch (<200ms dur)
% mid amp     6           25              24              27              52              57  exploratory
% big amp     6           54              19              25              54              82  long foveal
% small amp   7           12              25              28              51              54  short retracted
    
    
    % first def of peak Idx
    figure; hold on 
    plot(wAmplitude); plot(wAngle)
    peakIdx = wPhase(1:end-1)<0 & wPhase(2:end)>=0;
    plot(find(peakIdx),zeros(sum(peakIdx),1),'kd')
    wfreq=movsum(peakIdx,samplingRate);
    plot(wfreq)
    fAmp=WhiskingFun.GetAmplitude(wAngle,wPhase);
    
    % second def of peak Idx, based on velocity 
% Because velocity is highest during whisking, however, it is more common
% the estimate the whisking frequency from the angular velocity (??) which 
% contains a single harmonic typically denoted f0 (Berg and Kleinfeld 2003).
    figure; hold on
    plot(wVelocity)
    figure; hold on
    histogram(abs(wVelocity))

    whisksIdx = bwconncomp(wVelocity<0);
    peakIdx = zeros(length(wVelocity),1);
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
    plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    wfreq=movsum(peakIdx,samplingRate);
    figure; hold on
    plot(wfreq)
    
    figure; hold on 
    plot(wAngle)
    plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    plot(fAmp)
    
    figure; hold on
    histogram(wfreq(wAmplitude<=5))
    histogram(wfreq(wAmplitude>5 & wAmplitude<15))
    histogram(wfreq(wAmplitude>=15))
    legend('low amp','mid amp','high amp')
    title('Frequency distribution for different whisking amplitudes')
        
    figure; hold on
    histogram(wAmplitude(wfreq<=5))
    histogram(wAmplitude(wfreq>5 & wfreq<15))
    histogram(wAmplitude(wfreq>=15))
    legend('low freq','mid freq','high freq')
    title('Amplitude distribution for different whisking frequencies')
    
    figure; hold on
    histogram(wAmplitude(wSetPoint<70))
    histogram(wAmplitude(wSetPoint>=70))
    legend('low setpoint','high setpoint')
    title('Amplitude distribution for different set-points')
    
    figure; hold on
    histogram(wfreq(wSetPoint<70))
    histogram(wfreq(wSetPoint>=70))
    legend('low setpoint','high setpoint')
    title('Frequency distribution for different set-points')
    
    figure; hold on
    histogram(wSetPoint(wAmplitude<=5))
    histogram(wSetPoint(wAmplitude>5 & wAmplitude<15))
    histogram(wSetPoint(wAmplitude>=15))
    legend('low amp','mid amp','high amp')
    title('SetPoint distribution for different whisking amplitudes')

    figure; hold on
    histogram(wAmplitude(wSetPoint<70 & wAmplitude<15))
    histogram(wAmplitude(wSetPoint<70 & wAmplitude>=15))
    histogram(wAmplitude(wSetPoint>70 & wAmplitude>=15))
    legend('low SP low amp','low SP high amp','high SP high amp')
    title('Amplitude distribution for different whisking modes')

    
    %% plot angle and velocity by whisk duration
    % Figure from Knutsen's Whisking Kinematics
%     https://link.springer.com/chapter/10.2991/978-94-6239-133-8_46


    wAngle_nf=fillmissing(wAngle_nf,'linear');
    whiskerAngle_LP=WhiskingFun.LowPassBehavData(wAngle_nf,1000,40);
    whiskerVelocity_LP=diff(whiskerAngle_LP); wVelocity=[wVelocity(1) wVelocity]; 
    figure; hold on
    plot(whiskerAngle_LP)
    plot(wSetPoint)

    figure; hold on; plot(wPhase); plot(diff(wPhase))
    
    whisksIdx = bwconncomp(diff(wPhase)>0); %troughIdx = find(whiskerPhase(1:end-1)>=pi/2 & whiskerPhase(2:end)<=-pi/2);
    whiskDur=cellfun(@(whisk) length(whisk), whisksIdx.PixelIdxList);
    whisksIdx.PixelIdxList=whisksIdx.PixelIdxList(whiskDur>1); %prune
    whiskDur=cellfun(@(whisk) length(whisk), whisksIdx.PixelIdxList);

    [durBinCount,binEdges,binIdx]=histcounts(whiskDur,40);   
    binCat=unique(binIdx);
    [allCatMeanAngle,allCatMeanVel]=deal(nan(numel(binCat),100));
    for whiskDurCat=1:numel(binCat)
        catWhiskList=whisksIdx.PixelIdxList(binIdx==binCat(whiskDurCat));
%         catWhiskList=catWhiskList(cellfun(@numel, catWhiskList)>1); %prune out
        catMeanAngle=cellfun(@(catWhisk) interp1(wAngle(catWhisk),...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanAngle(whiskDurCat,:)=mean(vertcat(catMeanAngle{:}));
        catMeanVel=cellfun(@(catWhisk) interp1(wVelocity(catWhisk)*1000,...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanVel(whiskDurCat,:)=mean(vertcat(catMeanVel{:}));
    end
    binlabels=round(1000./round(binEdges+mode(diff(binEdges)/2)),1);
    figure; 
    subplot(1,2,1); 
    imagesc(allCatMeanAngle); cH = colorbar; set(cH,'Ylim',[40 110])
    ylabel('Whisk Frequency (Hz)'); xlabel('Normalized Time')
    title('\theta')
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    subplot(1,2,2); 
    imagesc(allCatMeanVel);  cH = colorbar; set(cH,'Ylim',[-10^3 10^3])
    ylabel('Whisk Frequency (Hz)'); xlabel('Normalized Time')
    title('\theta''')
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    
    
    %% using whisker frequency
%     whisksIdx = bwconncomp(whiskerVelocity<0);
    whisksIdx = bwconncomp(diff(wPhase)>0);
    peakIdx = zeros(length(wVelocity),1);
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
%     plot(find(peakIdx),zeros(sum(peakIdx),1),'gd')
    wfreq=movsum(peakIdx,samplingRate);
%     figure;  hold on; plot(wfreq)
    [freqBinCount,binEdges,binIdx]=histcounts(wfreq,20); 
    binCat=unique(binIdx);
    % index frequency to whisk
    whiskBinIdx=binIdx(logical(peakIdx));
    [allCatMeanAngle,allCatMeanVel]=deal(nan(numel(binCat),100));
    for whiskDurCat=1:numel(binCat)
        catWhiskList=whisksIdx.PixelIdxList(whiskBinIdx==binCat(whiskDurCat));
        catWhiskList=catWhiskList(cellfun(@numel, catWhiskList)>1); %prune out
        catMeanAngle=cellfun(@(catWhisk) interp1(whiskerAngle_LP(catWhisk),...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanAngle(whiskDurCat,:)=mean(vertcat(catMeanAngle{:}));
        catMeanVel=cellfun(@(catWhisk) interp1(whiskerVelocity_LP(catWhisk)*1000,...
            linspace(1,length(catWhisk),100)),catWhiskList, 'UniformOutput', false);
        allCatMeanVel(whiskDurCat,:)=mean(vertcat(catMeanVel{:}));
    end
    binlabels=round(binEdges+mode(diff(binEdges)/2),1);
    binlabels=fliplr(binlabels);
    figure; 
    subplot(1,2,1); 
    imagesc(allCatMeanAngle); cH = colorbar; set(cH,'Ylim',[40 110])
    ylabel('Whisk Frequency (Hz)'); xlabel('Normalized Time')
    title('\theta')
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    subplot(1,2,2); 
    imagesc(allCatMeanVel);  cH = colorbar; set(cH,'Ylim',[-10^3 10^3])
    ylabel('Whisk Frequency (Hz)'); xlabel('Normalized Time')
    title('\theta''')
    set(gca,'ydir','normal','yticklabel',binlabels(get(gca,'ytick')))
    

    %% find peak whisking frequency
        whiskFreqSpectrum=abs(fft(whiskerTraces_whiskFreq(whiskNum,whiskingEpochs)));
        numVals=numel(whiskerTraces_whiskFreq(whiskNum,whiskingEpochs));
        whiskFreqSpectrum = smoothdata(whiskFreqSpectrum(1:floor(numVals/2)+1),...
            'movmean',round(numel(whiskFreqSpectrum)/2/200));
        % peakWhiskFreq(2:end-1) = 2*peakWhiskFreq(2:end-1);
        freqArray = 1000*(0:(numVals/2))/numVals;
%       figure;  plot(freqArray,whiskFreqSpectrum)
        % title('Single-Sided Amplitude Spectrum of X(t)')
        peakWhiskFreq=freqArray(whiskFreqSpectrum==max(whiskFreqSpectrum));

