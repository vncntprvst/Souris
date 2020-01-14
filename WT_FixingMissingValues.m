%% Fixing missing values
whiskFrame = [whiskerData.fid];
frameUID=unique(whiskFrame);
missingFrames=find(diff(frameUID)>1)+1;
% frameUID(~ismember(frameUID,[currWhiskerData.fid]));
whiskID = double([whiskerData.wid]);
[whiskFreq,whiskUID]=hist(whiskID,unique(whiskID));

whiskFreq=whiskFreq(whiskUID>0); 
numWhisker = sum(whiskFreq/double(frameUID(end))>0.1); %need to be in more than 10% of frames
whiskUID=whiskUID(whiskUID>0);
   
locThd=2;

% start with the most frequent besides 0
whiskNum=1;
% currWhiskID=whiskUID(whiskNum);
% currWhiskerData=whiskerData([whiskerData.wid]==currWhiskID,:);

%%%%%%% do not use below %%%%%%%%
%% change ID by ordering whisker
% for frameNum=frameUID
%     if sum([whiskerData.fid]==frameNum)==numWhisker
%         %order according to follicle position
%         currFrameWhiskerData=whiskerData([whiskerData.fid]==frameNum);
%         [~,orderIdx]=sort([currFrameWhiskerData.follicle],'descend');
%         orderIdx=num2cell(orderIdx);
%         [whiskerData([whiskerData.fid]==frameNum).wid]=deal(orderIdx{:});
%     end
% end
% 
% %% Check and fix missing ones
% missingCWDataIdx=frameUID(~ismember(frameUID,...
%     [whiskerData([whiskerData.wid]==whiskUID(whiskNum),:).fid]));
% 
% for frameNum=fliplr(missingCWDataIdx) %(end) %go backward
%     currFrameWhiskerData=whiskerData([whiskerData.fid]==frameNum,:);
%     for postFN=frameNum+1:max(frameUID)
%         if ismember(postFN,frameUID) && ~ismember(postFN,missingCWDataIdx), break, end, end
%     try
%         postFrameWhiskerData=whiskerData([whiskerData.fid]==postFN,:); %&...
%         %             [whiskerData.wid]==currWhiskID,:);
%         fDist=nan(numel(postFrameWhiskerData),numel([currFrameWhiskerData.follicle]));
%         for compNum=1:numel(postFrameWhiskerData)
%             fDist(compNum,:)=abs([currFrameWhiskerData.follicle]-postFrameWhiskerData(compNum).follicle);
%         end
%         if sum(sum(fDist<locThd))
%             foundMatch=true;
%             % priority to potential match
%             %         supposedIDMatchIdx=[postFrameWhiskerData.wid]==currWhiskID;
%             %         putWhiskIdx=fDist(supposedIDMatchIdx,:)==min(fDist(supposedIDMatchIdx,:));
%             %         if fDist(supposedIDMatchIdx,putWhiskIdx)<locThd
%             
%             putWhiskIdx=currFrameWhiskerData(logical(sum(fDist==min(min(fDist))))).wid;
%             if putWhiskIdx>0
%                 %keep that in store and find which one is that one?
%             end
%             whiskerData([whiskerData.fid]==frameNum &...
%                 [whiskerData.wid]==putWhiskIdx,:).wid=whiskUID(whiskNum);
%         else
%             foundMatch=false; disp(['fail 1 - ' num2str(frameNum)]);
%         end
%     catch
%         foundMatch=false; disp(['fail 2 - ' num2str(frameNum)]);
%     end
%       
%     if foundMatch==false
%         for preFN=frameNum-1:-1:min(frameUID)
%             if ismember(preFN,frameUID) && ~ismember(preFN,missingCWDataIdx), break, end, end
%         try
%             preFrameWhiskerData=whiskerData([whiskerData.fid]==preFN,:);% &...
%             %                 [whiskerData.wid]==currWhiskID,:);
%             
%             fDist=nan(numel(preFrameWhiskerData),numel([currFrameWhiskerData.follicle]));
%             for compNum=1:numel(preFrameWhiskerData)
%                 fDist(compNum,:)=abs([currFrameWhiskerData.follicle]-preFrameWhiskerData(compNum).follicle);
%             end
%             if sum(sum(fDist<locThd))
%                 foundMatch=true;
%                 putWhiskIdx=currFrameWhiskerData(logical(sum(fDist==min(min(fDist))))).wid;
%                 if putWhiskIdx>0
%                     %keep that in store and find which one is that one?
%                 end
%                 whiskerData([whiskerData.fid]==frameNum &...
%                     [whiskerData.wid]==putWhiskIdx,:).wid=whiskUID(whiskNum);
%             else
%                 foundMatch=false; disp(['fail 3 - ' num2str(frameNum)]);
%             end
%         catch
%             foundMatch=false; disp(['fail 4 - ' num2str(frameNum)]);
%         end
%     end
%     if foundMatch==false
%         % see if needed
%     end
% end
%     