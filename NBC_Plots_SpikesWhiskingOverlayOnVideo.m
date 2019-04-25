function NBC_Plots_SpikesWhiskingOverlayOnVideo(periodBehavData_ms,vidTimes_ms,...
    buSpikeRasters_ms,buSDFs_ms,cursor_1,cursor_2);

%% get video file
videoFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.avi','*.mp4'},'UniformOutput', false);
videoFiles=vertcat(videoFiles{~cellfun('isempty',videoFiles)});

if ~isempty(videoFiles)
    vFileName=videoFiles.name;
    vDirName=videoFiles.folder;
else %ask for it
    [vFileName,vDirName] = uigetfile({'*.mp4;*.avi','video Files';'*.*','All Files'},...
        'Video file');
end
warning off MATLAB:subscripting:noSubscriptsSpecified
vidObj = VideoReader(fullfile(vDirName,vFileName));
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
vidStruc = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

%% whisker and spike data provided as inputs
% [recordingTraces,spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
%     SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
%     BP_periodBehavData_ms,HP_periodBehavData_ms,LP_periodBehavData_ms,...
%     HTBP_periodBehavData_ms, peakWhisking_ms,periodBehavData_ms,...
%     whiskingPhase_ms,instantFreq_ms,sgFreq,sgTime,sgPower,recName,vidTimes_ms] =...
%     NeuronBehaviorCorrelation_GatherData;

%% a few plots to check whisking and spiking 
% figure; plot(periodBehavData_ms(1,:))

% set initial time
if exist('cursor_1','var')
    timePoint=cursor_1.DataIndex;
else
    timePoint=700000; % x axis is in ms already
end
% timePoint=vidTimes_ms(timePoint); %video time in ms
% vidObj.FrameRate is by default 25. Convert to artificial time
timePoint=timePoint*20; % this is an average. Would be much better to match index to frame
minTP=floor(timePoint/60000);secTP=floor((timePoint - minTP*60000)/1000);
msTP=floor((timePoint - minTP*60000 - secTP*1000));
vidObj.CurrentTime = (minTP*60)+secTP+msTP*10^-3; %13:34
% Read one frame at a time using readFrame until the end of 3 sec epoch.
% Append data from each video frame to the structure array.
clipDuration=50/10^3*20;
% clipDuration=(cursor_2.DataIndex-cursor_1.DataIndex)/10^3*20;
k = 1;
while k<=(clipDuration*25)
    vidStruc(k).cdata = readFrame(vidObj);
    k = k+1;
end
numFrames=k;
% Display single frame stored in vidStruc
figure; image(vidStruc(1).cdata)
figure; image(vidStruc(end).cdata)

% Size the current figure and axes based on the video's width and height.
whiskingEpochh = figure;
set(whiskingEpochh,'position',[1050 150 vidObj.Width vidObj.Height]);
set(gca,'units','pixels');
set(gca,'position',[0 0 vidObj.Width vidObj.Height]);

% Play the movie once at the video's frame rate
movie(vidStruc,1,500);

%  rescale
periodBehavData_ms=(periodBehavData_ms-mean(periodBehavData_ms))/...
    max(abs(periodBehavData_ms))*50+50;

% PrV88_125_whiskers=LoadWhiskers('PrV88_125_HSCamClips2.whiskers');
numFrame=size(vidTimes_ms,2);
% M=nan()
figure;    hold on
for frameId = 74896:85079 %numFrame
    
    set video time
    
    hide axes
    
    s(frameId).cdata = readFrame(vidObj); %store image
    %% plot image
    image(flipud(s(frameId).cdata))
    
    %     fields=PrV88_125_whiskers([PrV88_125_whiskers.time]==frames(frameId));
    %     for framePlots=1:size(fields,1)
    %         plot(fields(framePlots).x,fields(framePlots).y);
    %     end
    % get time index
    timeIndices=vidTimes_ms(max([1 frameId-25])):vidTimes_ms(min([numFrame frameId+25]));
    timeIndices=int32(timeIndices(timeIndices>0));
    %% plot angle
    angleTrace=periodBehavData_ms(timeIndices);
    if frameId-25<=1
        angleTrace= [nan(1,101-numel(angleTrace)), angleTrace];
    elseif frameId+25>numFrame
        angleTrace= [angleTrace, nan(1,101-numel(angleTrace))];
    end
    timeAxis=int32(linspace(1,5*numel(angleTrace),numel(angleTrace)));
    plot(timeAxis,angleTrace,'k')
    
    axis([1 640 1 480])
    M(frameId) = getframe;
    cla;
end
figure
movie(M,1) %play 5 times;

%create video
v = VideoWriter('PrV88_125_whiskers_onvideo_171frames.avi');
open(v);
writeVideo(v, M);
close(v);
