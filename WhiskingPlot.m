% Plotting whisking angle trace

whiskerTrackDir='D:\Data\vIRt35\vIRt35_0724';
whiskerTrackFileName='vIRt35_0724_4300_20190724-140215_HSCam_WhiskersBaseCentroidTip.csv';
videoTimeStamps='vIRt35_0724_4300_20190724-140217_HSCam.csv';

whiskerTrackingData=readtable(fullfile(whiskerTrackDir,whiskerTrackFileName));
whiskerTrackingData=ContinuityWhiskerID(whiskerTrackingData);

% figure; hold on
% plot(whiskerTrackingData(:,1),whiskerTrackingData(:,2),'.')
% ax = gca;
% ax.YDir = 'reverse';

whiskerAngle=WhiskerAngleSmoothFill(whiskerTrackingData);
% whiskerAngle_ms=WhiskingAnalysisFunctions.ResampleBehavData...
%     (whiskerAngle,behav.vidTimes,ephys.samplingRate);

whiskerTraces_whiskFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(whiskerAngle,500,[4 20]);
whiskerTraces_breathFreq_ms=WhiskingAnalysisFunctions.BandPassBehavData(whiskerAngle,500,[1 4]);

timeStamps=readtable(fullfile(whiskerTrackDir,videoTimeStamps));
timeStamps=timeStamps.Var1(:,1)-timeStamps.Var2(:,1);
timeStamps=timeStamps/1000;

figure
timeWindow=8000:23000;
subplot(2,2,1:2); hold on
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerAngle(timeWindow),'k')
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerTraces_whiskFreq_ms(timeWindow),'b')
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerTraces_breathFreq_ms(timeWindow),'r')
xlabel('Time (s)')
ylabel('Angle (degrees)')
axis tight; box off
legend('raw trace','BP filtered 4-20Hz','BP filtered 1-4Hz','location','best')
subplot(2,2,3);hold on %zoom in #1 
timeWindow=11500:13000;
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerAngle(timeWindow),'k')
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerTraces_breathFreq_ms(timeWindow)+mean(whiskerAngle(timeWindow)),'r')
xlabel('Time (s)'); set(gca,'ylim',[75 120],'xticklabel',7:0.5:10); %'xlim',[7 10]
ylabel('Angle (degrees)')
axis tight; box off
subplot(2,2,4); hold on %zoom in #2
timeWindow=18000:19500;
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerAngle(timeWindow),'k')
plot(timeStamps(timeWindow)-timeStamps(timeWindow(1)),whiskerTraces_breathFreq_ms(timeWindow)+mean(whiskerAngle(timeWindow)),'r')
xlabel('Time (s)'); set(gca,'ylim',[75 120],'xticklabel',20:0.5:23); %'xlim',[20 23],
ylabel('Angle (degrees)')
axis tight; box off



