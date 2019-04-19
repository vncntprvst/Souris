function whiskerTrackingData=ContinuityWhiskerID(whiskerTrackingData)
% Finds continuous values for each identified whisker
% Outputs whisker centroids (x/y in alternate columns) for each frame
% (rows). 3rd dimension contains whisker "base" x/y values, to derive
% angle.

%% use base points to derive whisker id
emptyVariables=isnan(nanmean(whiskerTrackingData.Variables));
variableNames=whiskerTrackingData.Properties.VariableNames(~emptyVariables);
whiskerData=whiskerTrackingData(:,~emptyVariables).Variables;
baseData=whiskerData(:,cellfun(@(varName) contains(varName,'Base'),variableNames));
centroidXData=whiskerData(:,cellfun(@(varName) contains(varName,'Centroid_X'),variableNames));
centroidYData=whiskerData(:,cellfun(@(varName) contains(varName,'Centroid_Y'),variableNames));
baseXData=baseData(:,1:2:size(baseData,2));
baseYData=baseData(:,2:2:size(baseData,2));

% base continuity on base x or y, whichever has lower std
baseVar=[mean(nanstd(baseData(:,1:2:size(baseData,2)))),...
mean(nanstd(baseData(:,2:2:size(baseData,2))))];
contDirection=find(baseVar==min(baseVar));

cdBaseData=baseData(:,contDirection:2:size(baseData,2));
[corrCdBaseData,baseIdx]=sort(cdBaseData,2,'ascend'); %descend
% rearrange centroid values the same way
for rowNum=1:size(cdBaseData,1)
    centroidXData(rowNum,:)=centroidXData(rowNum,baseIdx(rowNum,:));
    centroidYData(rowNum,:)=centroidYData(rowNum,baseIdx(rowNum,:));
    baseXData(rowNum,:)=baseXData(rowNum,baseIdx(rowNum,:));
    baseYData(rowNum,:)=baseYData(rowNum,baseIdx(rowNum,:));
end
% 
lowestInitialW=find(~isnan(corrCdBaseData(1,:)),1,'last');
% find(isnan(corrCdBaseData(:,lowestInitialW)),1)
missingValPeriods=bwconncomp(bwlabel(isnan(corrCdBaseData(:,lowestInitialW))));
baseDiff=abs(diff(corrCdBaseData(:,lowestInitialW)));
for nanPeriod=1:missingValPeriods.NumObjects
    preNaNPVal=corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod}(1)-1,lowestInitialW);
    valDiff=abs(corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod}(1),:)-preNaNPVal);
    closeVals=valDiff<(3*nanstd(baseDiff));
    if logical(sum(closeVals))
        actualBase=find(valDiff==min(valDiff)); %find(closeVals,1,'last')
        baseGap=lowestInitialW-actualBase;
        if baseGap>0
            corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),baseGap) ...
                corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod},1:end-baseGap)];
            centroidXData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),baseGap) ...
                centroidXData(missingValPeriods.PixelIdxList{nanPeriod},1:end-baseGap)];
            centroidYData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),baseGap) ...
                centroidYData(missingValPeriods.PixelIdxList{nanPeriod},1:end-baseGap)];
            baseXData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),baseGap) ...
                baseXData(missingValPeriods.PixelIdxList{nanPeriod},1:end-baseGap)];
            baseYData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),baseGap) ...
                baseYData(missingValPeriods.PixelIdxList{nanPeriod},1:end-baseGap)];
        elseif baseGap<0
            corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [corrCdBaseData(missingValPeriods.PixelIdxList{nanPeriod},1+abs(baseGap):end) ...
                nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),abs(baseGap))];
            centroidXData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [centroidXData(missingValPeriods.PixelIdxList{nanPeriod},1+abs(baseGap):end) ...
                nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),abs(baseGap))];
            centroidYData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [centroidYData(missingValPeriods.PixelIdxList{nanPeriod},1+abs(baseGap):end) ...
                nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),abs(baseGap))];
            baseXData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [baseXData(missingValPeriods.PixelIdxList{nanPeriod},1+abs(baseGap):end) ...
                nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),abs(baseGap))];
            baseYData(missingValPeriods.PixelIdxList{nanPeriod},:)=...
                [baseYData(missingValPeriods.PixelIdxList{nanPeriod},1+abs(baseGap):end) ...
                nan(numel(missingValPeriods.PixelIdxList{nanPeriod}),abs(baseGap))];
        end
    end
end

%order them by the number of elements
[~,sortNonNan]=sort(sum(~isnan(centroidXData)),'descend');
whiskerTrackingData=[centroidXData(:,sortNonNan),centroidYData(:,sortNonNan)];
whiskerTrackingData=cat(3,whiskerTrackingData,[baseXData(:,sortNonNan),baseYData(:,sortNonNan)]);
whiskerTrackingData=whiskerTrackingData(:,reshape([1:size(baseData,2)/2;...
    size(baseData,2)/2+1:size(baseData,2)],1,size(baseData,2)),:);

%% plots
% figure; hold on
% for whiskBase=1:2:size(baseData,2)
%     plot(baseData(:,whiskBase),baseData(:,whiskBase+1),'.')
% end
% figure; hold on 
% for whiskBase=contDirection:2:size(baseData,2)
%     plot(baseData(:,whiskBase))
% end
% 
% figure; hold on 
% for whiskBase=6:-1:1
%     plot(corrCdBaseData(:,whiskBase))
% end
% 
% 
% baseIdx=nan(size(baseData,1),numel(contDirection:2:size(baseData,2)));
% corrBaseData=baseData;%(contDirection:2:size(baseData,2));
% for whiskBase=contDirection:2:contDirection*3 %size(corrBaseData,2)
%     whiskNum=(whiskBase-contDirection)/2+1;
%     baseIdx(1,whiskNum)=whiskNum; 
%     % find peaks
%     baseDiff=abs(diff(corrBaseData(:,whiskBase)));
% %     [pks,locs]=findpeaks(baseDiff,'MinPeakHeight',1.5*nanstd(baseDiff));
%     peakLoc=[0;0;abs(diff(baseDiff))]>2*nanstd(baseDiff);
%     
%     for peakNum=1:2:numel(locs)
%         frameIdx=locs(peakNum);
%         while frameIdx<=locs(peakNum+1)
%             if abs(corrBaseData(frameIdx,whiskBase)-corrBaseData(frameIdx-1,whiskBase))...
%                     <2*nanstd(diff(corrBaseData(:,whiskBase)))
%                 baseIdx(frameIdx,whiskNum)=whiskNum;
%             else %try to find another base
%                 correctBase=find(abs(...
%                     corrBaseData(frameIdx,contDirection:2:size(corrBaseData,2))-...
%                     corrBaseData(frameIdx-1,whiskBase))...
%                     <nanstd(diff(corrBaseData(:,whiskBase))),1);
%                 if ~isempty(correctBase)
%                     baseIdx(frameIdx,whiskNum)=correctBase;
%                     corrBaseData(frameIdx,whiskBase)=corrBaseData(frameIdx,correctBase*contDirection);
%                 else
%                     baseIdx(frameIdx,whiskNum)=NaN;
%                     corrBaseData(frameIdx,whiskBase)=corrBaseData(frameIdx-1,whiskBase);
%                 end
%             end
%             frameIdx=frameIdx+1;
%         end
%     end
% end
% 
% 3*std(diff(baseData(:,2)))
% 
% 
% changeIdx=ischange(baseData(:,contDirection:2:size(baseData,2)),'variance','Threshold',100*var(baseData(:,2)));
% plot(changeIdx(:,2)+170,'k')
% 
% plot([0; diff(baseData(:,contDirection))]+200);
% 
% for whiskBase=contDirection:2:size(baseData,2)
%     whiskNum=(whiskBase-contDirection)/2+1;
% 
%     if whiskBase==contDirection %first
%         baseNMO=[[baseData(:,whiskBase); baseData(end,whiskBase)],...
%             [baseData(1,whiskBase+2); baseData(:,whiskBase+2)]];
%                 baseNMO=baseData(:,2:2:whiskBase+2);
% 
%         NMOdifferential=diff(baseNMO,[],2);
%         [pks,locs] =findpeaks((diff(NMOdifferential)),'MinPeakHeight',3*std((diff(NMOdifferential))));
%         figure; hold on; 
%         plot((diff(NMOdifferential)))
%         plot(locs,pks,'dk')
%         
%         for peakNum=1:numel(pks)
%             if abs(baseNMO(locs(peakNum)+1,2)- baseNMO(locs(peakNum),1)) <=  3*std(diff(baseNMO(:,1)))
% %                 while
%             end
%         end
%                 
%                 
%     elseif whiskBase==size(baseData,2) %last
%         baseNMO=[baseData(1,whiskBase-2); baseData(:,whiskBase-2)],...
%             [[baseData(:,whiskBase); baseData(end,whiskBase)]];
%     else
%         baseNMO=[[baseData(1,whiskBase-2); baseData(:,whiskBase-2)],...
%             [baseData(:,whiskBase); baseData(end,whiskBase)],...
%             [baseData(1,whiskBase+2); baseData(:,whiskBase+2)]];
%     end
% end
% %     ...
% %     diff([[baseData(:,contDirection); baseData(end,contDirection)],...
% %     [baseData(1,contDirection+4); baseData(:,contDirection+4)]],[],2),...
% %     diff([[baseData(:,contDirection); baseData(end,contDirection)],...
% %     [baseData(1,contDirection+6); baseData(:,contDirection+6)]],[],2)];
% 
% figure; hold on 
% for whiskBaseDiffNMO=1:3
%     plot(baseNMO(:,whiskBaseDiffNMO)+200)
% end

