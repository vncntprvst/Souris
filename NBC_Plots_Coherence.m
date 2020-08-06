function NBC_Plots_Coherence(whiskerAngle,whiskerPhase,ephys,whiskingEpochs,selectedUnits,cmap)

fs = 1000; winL= 2000; overlapL= 1900; % nsc = floor(numel(whiskerAngle_m)/4.5); nov = floor(nsc/2);

whiskerAngle_m=fillmissing(whiskerAngle,'nearest');
wspectrumVals=spectrogram(whiskerAngle_m,winL,overlapL,[],fs,'yaxis');
wspectrumVals=abs(wspectrumVals);
wf_wps=sum(wspectrumVals(1:15,:)); %power spectrum at whisker frequency range
figure; axis2=subplot(2,1,1); hold on
plot(wf_wps)
%plot whisking epochs
misTvalSize=mod(numel(whiskingEpochs),numel(wf_wps));
epochTarray=linspace(1,numel(wf_wps),numel(whiskingEpochs)-misTvalSize);
epochVarray=whiskingEpochs(misTvalSize/2+1:end-misTvalSize/2)*nanstd(wf_wps)+nanmean(wf_wps);
plot(epochTarray,epochVarray,'r','linewidth',1.5)
xlabel('Whisking power spectrum')

for unitNum=1:numel(selectedUnits)
    unitRasters=ephys.rasters(selectedUnits(unitNum),:);
    spectrumVals_s=spectrogram(unitRasters,winL,overlapL,[],fs,'yaxis');
    spectrumVals_s=abs(spectrumVals_s);
    sum_spectrumVals_s(unitNum,:)=sum(spectrumVals_s(1:15,:));
end
wf_sps=mean(sum_spectrumVals_s);%[7,10,12,13,15] %sum(sspectrumVals(1:15,:)); %power spectrum at whisker frequency range
axis1=subplot(2,1,2);  hold on
plot(wf_sps)
epochVarray=whiskingEpochs(misTvalSize/2+1:end-misTvalSize/2)*nanstd(wf_sps)+nanmedian(wf_sps);
plot(epochTarray,epochVarray,'r','linewidth',1.5)
xlabel('Spiking power spectrum')
linkaxes([axis1,axis2],'x')

figure; %imagesc(sspectrumVals)
spectrogram(unitRasters,winL,overlapL,[],fs,'yaxis');
set(gca,'ylim',[1 30]); %'YDir','normal'); 
caxis([-60 -25]);hold on
epochTarray=linspace(1,max(get(gca,'xlim')),numel(whiskingEpochs));
epochVarray=whiskingEpochs*30;
plot(epochTarray,epochVarray,'r','linewidth',1.5)

% coherence within each whisking epoch
whiskerPhase_m=fillmissing(whiskerPhase,'nearest');
whiskBoutList = bwconncomp(whiskingEpochs) ; 
durationThd=cellfun(@(x) length(x),whiskBoutList.PixelIdxList)>4000; %keep only whisking bouts longer than 4s 
whiskBoutList.PixelIdxList=whiskBoutList.PixelIdxList(durationThd);
whiskBoutList.NumObjects=sum(durationThd);
figure('color','w','name','Spike / Angle Coherence','position',[1000 70 400 1200]); hold on
histHand = gobjects(2,whiskBoutList.NumObjects);
cxy=nan(257,whiskBoutList.NumObjects);
binNum=10; shiftF=0.6;
for wEpoch=whiskBoutList.NumObjects:-1:1
    signalVals=whiskerPhase_m(whiskBoutList.PixelIdxList{wEpoch}); %whiskerAngle_m whiskerPhase_m
    signalVals=signalVals-mean(signalVals);
%     signalVals=detrend(unwrap(signalVals))-mean(detrend(unwrap(signalVals)));
    cxy(:,wEpoch)=mscohere(unitRasters(whiskBoutList.PixelIdxList{wEpoch}),...
        signalVals,500,250,[],fs);
    freqVals=linspace(0,fs/2,numel(cxy(:,wEpoch)));
    xVals = [freqVals freqVals(end) 0]; 
    yVals = [cxy(:,wEpoch)',0,0]; 
    fill(gca,xVals,yVals+(wEpoch*shiftF)-1,cmap(wEpoch,:),'FaceAlpha',0.4,'EdgeColor','w','LineWidth',1);
end
% delete(histHand)
set(gca,'xlim',[0 50],'ytick',(shiftF:shiftF:(shiftF*whiskBoutList.NumObjects))-0.5,...
    'ylim',[-0.5 whiskBoutList.NumObjects*shiftF],'yticklabel',1:whiskBoutList.NumObjects,...
    'tickdir','out','box','off','XGrid', 'on');
ylabel('Whisking epoch')
xlabel('Frequency (Hz)')
title('Spike / Phase Magnitude-squared Coherence')
% legend(num2str((1:whiskBoutList.NumObjects)'))

%ridge plot
lgd=mat2cell(num2str((1:whiskBoutList.NumObjects)'),ones(whiskBoutList.NumObjects,1),2);
ridgeplot(mat2cell(cxy',ones(whiskBoutList.NumObjects,1),257),binNum,[0 0.2],lgd,cmap)


