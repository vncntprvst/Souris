function slopeCoeff=ExtractMultiWhiskerAngle_FFTonContours(fileName)

whiskerCountoursVideo=VideoReader(fileName);
vHeight=whiskerCountoursVideo.Height;
vWidth=whiskerCountoursVideo.Width;
numReadFrames=30000;
% currAxes = axes;
% vidFrame=zeros(vHeight,vWidth,numReadFrames);
slopeCoeff=zeros(numReadFrames,1);
frameNum=1;
while hasFrame(whiskerCountoursVideo) %frameNum<=numReadFrames
    thatFrame = readFrame(whiskerCountoursVideo);
    vidFrame=thatFrame(:,:,1); %(:,:,frameNum)
    %     image(vidFrame, 'Parent', currAxes);
    %     currAxes.Visible = 'off';
    %     pause(1/whiskerVideo.FrameRate);
    
    %% get major whisker angle through Fourier transform
    whiskerFT=fft2(vidFrame); %(:,:,frameNum)); %Fourier transform 
    shiftedFreq = abs(fftshift(whiskerFT)); %Shift zero-frequency component to center of spectrum
%     figure; imagesc(shiftedFreq);
    
    [~,maxFreqIdx] = max(shiftedFreq); % Indices of max values
    maxFreqIdx = maxFreqIdx - maxFreqIdx(vWidth/2); %Center indices
    maxFreqIdx([1:100,vWidth-100:vWidth]) = []; % crop outliers
    newXLim=floor(length(maxFreqIdx)/2); xAxis = -newXLim:newXLim; %redefine axis
    slopeCoeff(frameNum) = regress(maxFreqIdx',xAxis'); %get regression slope 
    
    % figure;
    % scatter(xAxis,maxFreqIdx);
    % hold on;
    % plot(xAxis,slopeCoeff(frameNum)*xAxis);
    
    frameNum=frameNum+1;
end
% 
% figure;
% imagesc(vidFrame(:,:,end))






