
% Start in WhiskerTracking folder (contains measurement files)
sessionDir=cd;
ext = '.measurements';
include_files = arrayfun(@(x) x.name(1:(end-length(ext))),...
    dir([sessionDir filesep '*' ext]),'UniformOutput',false);

if ~exist('whiskerPadCoordinates','var')
    whiskerPadCoordinates=NaN;
else
    [faceSideInImage,protractionDirection,linkingDirection]=...
        GetWhiskerPadParams(whiskerPadCoordinates,whiskerPadRelativeLocation);
end


% Apply whisker linker to better sort out whisker id
for fileNum=1:numel(include_files)
    try
        %Load measurements
        measurementsPath=[include_files{fileNum} '.measurements'];
        trialMeasurements = Whisker.LoadMeasurements(measurementsPath);
        %Load whisker
        whiskerPath=[include_files{fileNum} '.whiskers'];
        trialWhiskers = Whisker.LoadWhiskers(whiskerPath);
        
    catch
        continue
    end
    
        copyfile(fullfile(cd,measurementsPath),fullfile(cd,'OriginalMeasurements',measurementsPath));

        % Link across frames
    linkedMeasurements = WhiskerLinkerLite(trialMeasurements,...
        'whiskerpadROI',whiskerPadCoordinates,...
        'linkingDirection',linkingDirection,...
        'faceSideInImage',faceSideInImage,...
        'protractionDirection',protractionDirection);
    
    %% plot individual frames
%     video=VideoReader([include_files{fileNum} '.mp4']);
%     vidFrame = readFrame(video);
%     figure; hold on
%     image(vidFrame); set(gca, 'YDir', 'reverse');
%     rectangle('Position',whiskerPadCoordinates,'EdgeColor','w')
%     
%     fIdx=[trialMeasurements.fid] == 2;
%     frMeasures = trialMeasurements(fIdx);
%     
%     % figure; hold on
%     plot([frMeasures.follicle_x],[frMeasures.follicle_y],'.')
%     plot([frMeasures.tip_x],[frMeasures.tip_y],'x')
    
    %% whole trial operations
    blacklist = ... %[ trialMeasurements.length ] < whiskerLengthThresh | ...
        [ trialMeasurements.follicle_x ] > whiskerPadCoordinates(1)+whiskerPadCoordinates(3) | ...
        [ trialMeasurements.follicle_x ] < whiskerPadCoordinates(1) | ...
        [ trialMeasurements.follicle_y ] > whiskerPadCoordinates(2)+whiskerPadCoordinates(4) | ...
        [ trialMeasurements.follicle_y ] < whiskerPadCoordinates(2);
    
    wpTrialMeasurements = trialMeasurements(~blacklist);
    wIdx=[wpTrialMeasurements.label];
    
    %% plot original clusters
%     figure; hold on
%     image(vidFrame); set(gca, 'YDir', 'reverse');
%     plot([wpTrialMeasurements.follicle_x],[wpTrialMeasurements.follicle_y],'.')
%     plot([wpTrialMeasurements.tip_x],[wpTrialMeasurements.tip_y],'x')
%     
%     figure; hold on
%     image(vidFrame); set(gca, 'YDir', 'reverse');
%     for wNum=0:2
%         plot([wpTrialMeasurements(wIdx==wNum).follicle_x],[wpTrialMeasurements(wIdx==wNum).follicle_y],'.')
%         plot([wpTrialMeasurements(wIdx==wNum).tip_x],[wpTrialMeasurements(wIdx==wNum).tip_y],'x')
%     end
%     
%     figure; hold on
%     image(vidFrame); set(gca, 'YDir', 'reverse');
%     plot([wpTrialMeasurements(wIdx==-1).follicle_x],[wpTrialMeasurements(wIdx==-1).follicle_y],'.')
%     plot([wpTrialMeasurements(wIdx==-1).tip_x],[wpTrialMeasurements(wIdx==-1).tip_y],'x')
            
    %% get distance to original cluster
    initCLuster=[[wpTrialMeasurements.follicle_x]',...
        [wpTrialMeasurements.follicle_y]',...
        [wpTrialMeasurements.tip_x]',...
        [wpTrialMeasurements.tip_y]'];
    unsortedObs=initCLuster(wIdx==-1,:);
    [~,ic_comps] = pca(zscore(initCLuster));
    [~,uo_comps] = pca(zscore(unsortedObs));
    for wNum=0:2
        %     mDist=squareform(pdist([mean(initCLuster(wIdx==wNum,:));unsortedObs]));
        %     mDist2Clus(:,wNum+1)=mDist(2:end,1);
        mDist=squareform(pdist([mean(ic_comps(wIdx==wNum,:));uo_comps]));
        mDist2Clus(:,wNum+1)=mDist(2:end,1);
    end
    clusAlloc=mod(find((mDist2Clus==min(mDist2Clus,[],2))'),3)';
    clusAlloc(clusAlloc==0)=3; clusAlloc=clusAlloc-1;
    
    % plot outcome
%     figure; hold on
%     image(vidFrame); set(gca, 'YDir', 'reverse');
%     for wNum=0:2
%         plot([unsortedObs(clusAlloc==wNum,1)],[unsortedObs(clusAlloc==wNum,2)],'.')
%         plot([unsortedObs(clusAlloc==wNum,3)],[unsortedObs(clusAlloc==wNum,4)],'x')
%     end
    
    %% allocate to frames that need it
    unAllocObs=find(wIdx==-1);
    wpIdx=find(~blacklist);
    for wNum=0:2
        fidList=unAllocObs(clusAlloc==wNum);
        nonRedundent_fidList=~ismember([wpTrialMeasurements(fidList).fid],...
            [trialMeasurements(wpIdx(wIdx==wNum)).fid]);
        fidList=fidList(nonRedundent_fidList);
        [trialMeasurements(wpIdx(fidList)).label]=deal(wNum);
    end
    
    %% Save output
    Whisker.SaveMeasurements(measurementsPath,linkedMeasurements.outMeasurements);
    
end