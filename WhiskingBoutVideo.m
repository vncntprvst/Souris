function [wBoutFrames,frameIndex]=WhiskingBoutVideo(vFileName,vDirName,boutIndex,frameTimes,saveVid)

if saveVid
    recName=regexp(cd,['(?<=\' filesep ')\w+$'],'match','once');
end
%% get video file
if isempty(vFileName)
    [vFileName,vDirName] = uigetfile({'*.avi';'*.mp4'},'Select Video File');
end

warning off MATLAB:subscripting:noSubscriptsSpecified
vidObj = VideoReader(fullfile(vDirName,vFileName));
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
wBoutFrames = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

% set initial time
if ~exist('boutIndex','var') || isempty(boutIndex)
    boutIndex=true(size(frameTimes));
end

frameIndex=find(frameTimes>=boutIndex(1) & frameTimes<=boutIndex(end)); 
% frameDur = mode(diff(frameTimes)); 
% frameRate = 1/frameDur*1000;
% videoTimes=frameTimes/(vidObj.FrameRate*frameDur); %basically, /50

vidObj.CurrentTime = frameIndex(1)/vidObj.FrameRate; %videoTimes(frameIndex(1));% vidObj.CurrentTime = frameIndex(1)/frameRate*vidObj.FrameRate;

% Read one frame at a time using readFrame until the end of epoch.
% Append data from each video frame to the structure array.
clipDuration=numel(frameIndex); %/(vidObj.FrameRate*frameDur); 
% clipDuration=(cursor_2.DataIndex-cursor_1.DataIndex)/10^3*20;
frameNum = 1;
while frameNum<=clipDuration
    wBoutFrames(frameNum).cdata = readFrame(vidObj); %conversion to grayscale doesn't help with video size
    frameNum = frameNum+1;
end

%% Display movie
% Size the  figure and axes based on the video's width and height.
% figure('position',[1050 150 vidObj.Width vidObj.Height]);
% set(gca,'units','pixels');
% movie(vidStruc,1,frameRate);

%% Display single frame stored in vidStruc
% figure; image(vidStruc(end).cdata)
% figure('position',[1050 150 vidObj.Width vidObj.Height]);
% FrameByFrame(vidStruc);

%% save video
if saveVid
    vidOutObj = VideoWriter(['WhiskingBout_' recName '_Slow10x.avi']);
    vidOutObj.FrameRate=50;
    open(vidOutObj);
    writeVideo(vidOutObj, wBoutFrames);
    close(vidOutObj);
end

% vidObj = VideoWriter(path);
%             vidObj.FrameRate = 50;
%             open(vidObj);
%             
%             figure;
%             axes('Parent', gcf, 'Units', 'pixels', 'Position', [ 0 1 640 480 ]);
%             
%             hold on
%             ih = imagesc(this.movie(:,:,1), [0 255]);
%             colormap gray
%             
%             ch = plot(this.ptMasks(1,:,1), this.ptMasks(1,:,2), 'r', 'LineWidth', 2);
%             axis ij off
%             xlim([1 size(this.movie,2)]);
%             ylim([1 size(this.movie,1)]);
%             hold off
%             
%             for i = 1 : size(this.ptMasks, 1)
%                 set(ih, 'CData', this.movie(:,:,i));
%                 set(ch, 'XData', this.ptMasks(i,:,1));
%                 set(ch, 'YData', this.ptMasks(i,:,2));
%                 drawnow;
%                 
%                 frameObj = getframe;
%                 writeVideo(vidObj, frameObj.cdata);
%             end
% 
%             close(vidObj);

