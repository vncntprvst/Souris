% All timestamps are in seconds
% Video frame timestamps are stored in behav.vidTimes
% Whisker data is in behav.whiskers or whiskers 
% Spike timestamps are in ephys.timestamps

% bWhisk is an index for best whisker trace

%%  whisking bouts are previously detected as follow:
% ampThd=18; %12; %amplitude threshold
% freqThld=1; %frequency threshold
% minBoutDur=1000; %500; % minimum whisking bout duration: 1s
% whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
%     whiskers(bWhisk).amplitude,whiskers(bWhisk).frequency,...
%     ampThd, freqThld, minBoutDur);
% whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
% whiskingEpochsList=bwconncomp(whiskingEpochs);
% [~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList.PixelIdxList),'descend');
% whiskingEpochsList.PixelIdxListSorted=whiskingEpochsList.PixelIdxList(wBoutDurSort);


%% select whisking bout and cell number
boutNum=7; %19 7 4
cellNum=37; %37 26 54

%% Get behavior video frames
traceIndex=whiskingEpochsList.PixelIdxList{boutNum}/1000+behav.vidTimes(1)-0.001;
[boutFrames,frameIndex]=WhiskingBoutVideo(ephys.recInfo.likelyVideoFile,ephys.recInfo.dirName,...
    traceIndex,behav.vidTimes,false);
%remove any drift by adding frames when fps > intended fps  (crude)
boutTimeLine=behav.vidTimes(frameIndex);boutTimeLine=boutTimeLine-boutTimeLine(1);
extraFrameIdx=find(diff(round(boutTimeLine))>mode(diff(round(boutTimeLine))));
for frameIdx=1:numel(extraFrameIdx)
    boutFrames=[boutFrames(1:extraFrameIdx(frameIdx)+frameIdx-1),...
        boutFrames(extraFrameIdx(frameIdx)+frameIdx-1),...
        boutFrames(extraFrameIdx(frameIdx)+frameIdx:end)];
end
boutFrames=boutFrames(1:end-numel(extraFrameIdx));
vidDims=size(boutFrames(1).cdata);

%% Get corresponding spikes and overlay
traceTS=behav.whiskers(bWhisk).timestamp(whiskingEpochsList.PixelIdxList{boutNum}+2);
spikeIndex=ephys.timestamps>=traceTS(1) & ephys.timestamps<=traceTS(end);
spikeTimes = movmean(ephys.rasters(cellNum,spikeIndex(1:size(ephys.rasters,2))),2);
spikeTimes = logical(spikeTimes(1:2:end));
figure('position',[500 450  vidDims(2) vidDims(1)],'color','k');
boutFrames=FrameByFrame_Overlay(boutFrames,vertcat(whiskers(bWhisk).phase...
    (whiskingEpochsList.PixelIdxList{boutNum}(1:2:end)+2),spikeTimes)); %whisker.angle

%% add audio for a given cell, given speed
slowFactor=20; %500/25
wBoutAudio = MakeSpikeAudio(ephys.rasters(cellNum,spikeIndex(1:size(ephys.rasters,2))),...
    slowFactor*10, [1 2500], 1000);

%% write video
frameRate=500/slowFactor;
videoFWriter = vision.VideoFileWriter(fullfile(cd,...
    [ephys.recInfo.likelyVideoFile(1:end-4) '_Bout' num2str(boutNum)...
    'x' num2str(slowFactor) '_Unit' num2str(cellNum) '_FP_PhaseTuning.avi']));
videoFWriter.FrameRate =frameRate ;
videoFWriter.AudioInputPort = true;
videoFWriter.VideoCompressor = 'None (uncompressed)'; % 'MJPEG Compressor';

for frameNum = 1 : 2500 %numel(boutFrames)
    videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(:,frameNum)); %if matrix form factor
end
release(videoFWriter);


function spikeAudio = MakeSpikeAudio(spikeRasters, slowFactor, frameRange, sampleRate)

spikeBinary = spikeRasters(:);
waveforms = repmat(-gausswin(8), 1, sum(spikeBinary));
spikeTimesScaled = find(spikeBinary) * slowFactor;
spikeAudio = zeros(length(spikeBinary)*slowFactor, 1);


for j = 1 : size(waveforms,2)
    spikeAudio(spikeTimesScaled(j):spikeTimesScaled(j)+size(waveforms,1)-1) = waveforms(:,j);
end

spikeAudio = interp1(1:length(spikeAudio),spikeAudio,1:sampleRate*5*slowFactor);
spikeAudio = reshape(spikeAudio, [ round(length(spikeAudio)/2500), 2500 ]);
spikeAudio = spikeAudio(:, frameRange(1):frameRange(2));
spikeAudio = int16(spikeAudio / max(abs(spikeAudio(:))) * double(intmax('int16')));

end