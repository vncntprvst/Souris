
numElectrodes=32;
recDuration=allRecInfo{6, 1}.dur;
allTraces=reshape(ephys.traces,[numElectrodes recDuration]);
filtTraces=PreProcData(allTraces,ephys.spikes.samplingRate,{'CAR','all'});

w1Centroid=behav.whiskerTrackingData.Whisker1_Centroid_Y;
w2Centroid=behav.whiskerTrackingData.Whisker2_Centroid_Y;
w3Centroid=behav.whiskerTrackingData.Whisker3_Centroid_Y;
% numel(w2Centroid)*60
keepUnits =[7     8    15    16    21 52    63    66];

fps=500;
sampRate=30000;
recMidPoint=round(size(filtTraces,2)/2);
wMidPoint=round(size(w2Centroid,1)/2);


for unitNum=1:numel(keepUnits)
    figure('position',[36 767 1541 331]);
    
    unitIdx=ephys.spikes.unitID==keepUnits(unitNum);
    ChNum=mode(ephys.spikes.preferredElectrode(unitIdx));
    spikeTimes=ephys.spikes.times(unitIdx);
    subplot(2,1,1);    hold on
    plot(w2Centroid(wMidPoint-30*fps:wMidPoint+30*fps))
    convST=spikeTimes/(sampRate/fps);
    tmConvST=convST(convST>=wMidPoint-30*fps & convST<=wMidPoint+30*fps)-wMidPoint+30*fps;
    plot(tmConvST,ones(numel(tmConvST,1))*nanmean(w2Centroid),'dk')
    
    tmST=spikeTimes(spikeTimes>=recMidPoint-30*sampRate &...
        spikeTimes<=recMidPoint+30*sampRate)-recMidPoint+30*sampRate;
    subplot(2,1,2);    hold on
    plot(filtTraces(ChNum,recMidPoint-30*sampRate:recMidPoint+30*sampRate))
    plot(tmST,ones(numel(tmST,1))*-300,'dk')
end
