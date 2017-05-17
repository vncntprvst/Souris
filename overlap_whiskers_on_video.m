%original video
vidObj = VideoReader('PrV88_125_HSCamClips2.avi');
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

%whisker data
PrV88_125_whiskers=LoadWhiskers('PrV88_125_HSCamClips2.whiskers');
frames=min([PrV88_125_whiskers.time]):max([PrV88_125_whiskers.time]);

figure;
for frameId = 1:max(frames)
    hold on
    s(frameId).cdata = readFrame(vidObj); %store image
    image(s(frameId).cdata)
    fields=PrV88_125_whiskers([PrV88_125_whiskers.time]==frames(frameId));
    for framePlots=1:size(fields,1)
        plot(fields(framePlots).x,fields(framePlots).y);
    end
	axis([1 640 1 480])
	M(frameId) = getframe;
    cla;
end
figure
movie(M,5) %play 5 times;

%create video
v = VideoWriter('PrV88_125_whiskers_onvideo_171frames.avi');
open(v);
writeVideo(v, M);
close(v);
