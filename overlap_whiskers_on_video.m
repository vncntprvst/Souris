function overlap_whiskers_on_video(fName,numFrames,playVideo,saveVideo)

cmap=lines;

%original video
vidObj = VideoReader([fName '.mp4']);
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

%whisker data
whiskers=Whisker.LoadWhiskers([fName '.whiskers']);
frames=min([whiskers.time]):max([whiskers.time]);

%measurements
measurements = Whisker.LoadMeasurements([fName '.measurements']);

switch nargin
    case 1
        numFrames=max(frames); [playVideo,saveVideo] = deal(false);
    case 2
        [playVideo,saveVideo] = deal(false);
    case 3
        saveVideo = false;
end

figure;
for frameId = 1:numFrames
    hold on
    s(frameId).cdata = readFrame(vidObj); %store image
    image(s(frameId).cdata)
    %plot whiskers
    fields=whiskers([whiskers.time]==frames(frameId));
    for framePlots=1:size(fields,1)
        plot(fields(framePlots).x,fields(framePlots).y,'color',cmap(framePlots,:));
    end
    % then plot follicles from measurments
    fields=measurements([measurements.fid]==frames(frameId));
    for framePlots=1:size(fields,1)
        plot(measurements(framePlots).follicle_x,...
            measurements(framePlots).follicle_y,...
            'marker','O','MarkerEdgeColor','none','MarkerFaceColor',cmap(framePlots,:))
    end
    axis([1 640 1 480])
    if playVideo || saveVideo
        M(frameId) = getframe;
        cla;
    end
end

if playVideo
    figure
    movie(M,5) %play 5 times;
end

%create video
if saveVideo
    v = VideoWriter('whiskers_on_video.avi');
    open(v);
    writeVideo(v, M);
    close(v);
end
