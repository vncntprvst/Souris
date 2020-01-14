function NBC_Plots_WhiskingPhaseVideoFrame(whiskerTrackingData)

%% get video file
[vFileName,vDirName] = uigetfile({'*.avi';'*.mp4'},'Select Video File');

if ~exist('whiskerTrackingData','var')
    [whiskerTrackingFile,bDirName] = uigetfile({'*.csv'},'Select Behavior File');
    whiskerTrackingData=readtable(fullfile(bDirName,whiskerTrackingFile));
    whiskerTrackingData=ContinuityWhiskerID(whiskerTrackingData);
    whiskerTrackingData=WhiskerAngleSmoothFill(whiskerTrackingData);
end

warning off MATLAB:subscripting:noSubscriptsSpecified
vidObj = VideoReader(fullfile(vDirName,vFileName));
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
vidStruc = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);


%% a few plots to check whisking and spiking 
% figure; plot(whiskerTrackingData(1,:))

% set initial time
if exist('cursor_1','var')
    timePoint=cursor_1.DataIndex;
else
    timePoint=30000; % x axis is in ms already
end

%% if whisker data is in ms  
% vidObj.FrameRate is by default 25. Convert to artificial time
% timePoint=timePoint*20; % this is an average. Would be much better to match index to frame
% minTP=floor(timePoint/60000);secTP=floor((timePoint - minTP*60000)/1000);
% msTP=floor((timePoint - minTP*60000 - secTP*1000));
% vidObj.CurrentTime = (minTP*60)+secTP+msTP*10^-3; %13:34
frameRate = 500;
vidObj.CurrentTime = timePoint/frameRate*vidObj.FrameRate;

% Read one frame at a time using readFrame until the end of 3 sec epoch.
% Append data from each video frame to the structure array.
clipDuration=3; %1000/10^3*20;
% clipDuration=(cursor_2.DataIndex-cursor_1.DataIndex)/10^3*20;
k = 1;
while k<=(clipDuration*frameRate)
    vidStruc(k).cdata = readFrame(vidObj);
    k = k+1;
end
figure;movie(vidStruc,1,50);
% Display single frame stored in vidStruc
figure; image(vidStruc(1).cdata)
figure; image(vidStruc(end).cdata)

figure; plot(whiskerTrackingData(1,timePoint:timePoint+clipDuration*frameRate))
timePoint2=cursor_2.DataIndex;
timePoint3=cursor_3.DataIndex;
figure; hold on 
for frameNum=0:timePoint3-timePoint2
    subplot(4,3,frameNum+1);
    image(vidStruc(timePoint2+frameNum).cdata);
end


% Size the current figure and axes based on the video's width and height.
whiskingEpochh = figure;
set(whiskingEpochh,'position',[1050 150 vidObj.Width vidObj.Height]);
set(gca,'units','pixels');
set(gca,'position',[0 0 vidObj.Width vidObj.Height]);

% Play the movie once at the video's frame rate
figure; 
movie(vidStruc(timePoint2:timePoint3),1,50);
