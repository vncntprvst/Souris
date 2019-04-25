function NBC_Plots_SpikesWhiskingOverlayOnVideo


%original video
[vFileName,vDirName] = uigetfile({'*.mp4;*.avi','video Files';'*.*','All Files'},...
    'Video file');
vidObj = VideoReader(fullfile(vDirName,vFileName));
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

%whisker data
[recordingTraces,spikeRasters_ms,rasterXInd_ms,rasterYInd_ms,samplingRate,...
    SDFs_ms,spikeTimes,waveForms,unitID,preferredElectrode,keepUnits,...
    BP_periodBehavData_ms,HP_periodBehavData_ms,LP_periodBehavData_ms,...
    HTBP_periodBehavData_ms, peakWhisking_ms,periodBehavData_ms,...
    whiskingPhase_ms,instantFreq_ms,sgFreq,sgTime,sgPower,recName,vidTimes_ms] =...
    NeuronBehaviorCorrelation_GatherData;

%   rescale
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
