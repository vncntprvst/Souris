function whiskerTrackingData=ContinuityWhiskerID(whiskerTrackingData)
%% use base points to derive whisker id
emptyVariables=isnan(nanmean(whiskerTrackingData.Variables));
variableNames=whiskerTrackingData.Properties.VariableNames(~emptyVariables);
whiskerData=whiskerTrackingData(:,~emptyVariables).Variables;
baseData=whiskerData(:,cellfun(@(varName) contains(varName,'Base'),variableNames));

% base continuity on base x or y, whichever has lower std
baseVar=[nanmean(std(baseData(:,1:2:size(baseData,2)))),...
nanmean(std(baseData(:,2:2:size(baseData,2))))];
contDirection=find(baseVar==min(baseVar));

% plots
figure; hold on
for whiskBase=1:2:size(baseData,2)
    plot(baseData(:,whiskBase),baseData(:,whiskBase+1),'.')
end
figure; hold on 
for whiskBase=contDirection:2:size(baseData,2)
    plot(baseData(:,whiskBase))
end

3*std(diff(baseData(:,2)))


changeIdx=ischange(baseData(:,contDirection:2:size(baseData,2)),'variance','Threshold',100*var(baseData(:,2)));
plot(changeIdx(:,2)+170,'k')


for whiskBase=contDirection:2:size(baseData,2)
    whiskNum=(whiskBase-contDirection)/2+1;

    if whiskBase==contDirection %first
        baseNMO=[[baseData(:,whiskBase); baseData(end,whiskBase)],...
            [baseData(1,whiskBase+2); baseData(:,whiskBase+2)]];
                baseNMO=baseData(:,2:2:whiskBase+2);

        NMOdifferential=diff(baseNMO,[],2);
        [pks,locs] =findpeaks((diff(NMOdifferential)),'MinPeakHeight',3*std((diff(NMOdifferential))));
        figure; hold on; 
        plot((diff(NMOdifferential)))
        plot(locs,pks,'dk')
        
        for peakNum=1:numel(pks)
            if abs(baseNMO(locs(peakNum)+1,2)- baseNMO(locs(peakNum),1)) <=  3*std(diff(baseNMO(:,1)))
                while
                
                
    elseif whiskBase==size(baseData,2) %last
        baseNMO=[baseData(1,whiskBase-2); baseData(:,whiskBase-2)],...
            [[baseData(:,whiskBase); baseData(end,whiskBase)]];
    else
        baseNMO=[[baseData(1,whiskBase-2); baseData(:,whiskBase-2)],...
            [baseData(:,whiskBase); baseData(end,whiskBase)],...
            [baseData(1,whiskBase+2); baseData(:,whiskBase+2)]];
    end
%     ...
%     diff([[baseData(:,contDirection); baseData(end,contDirection)],...
%     [baseData(1,contDirection+4); baseData(:,contDirection+4)]],[],2),...
%     diff([[baseData(:,contDirection); baseData(end,contDirection)],...
%     [baseData(1,contDirection+6); baseData(:,contDirection+6)]],[],2)];

figure; hold on 
for whiskBaseDiffNMO=1:3
    plot(baseNMO(:,whiskBaseDiffNMO)+200)
end

