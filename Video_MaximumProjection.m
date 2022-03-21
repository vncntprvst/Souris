function Video_MaximumProjection
% use WhiskerDiff.bonsai workflow to create video
%% get video file
videoFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.avi','*.mp4'},'UniformOutput', false);
videoFiles=vertcat(videoFiles{~cellfun('isempty',videoFiles)});

if ~isempty(videoFiles)
    if size(videoFiles,1)>1
        [~,timeIdx]=sort({videoFiles.date});
        videoFiles=videoFiles(timeIdx(end));
    end
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

sumVidIdx=zeros(vidHeight,vidWidth,'uint32');

boutStart=whiskingEpochsList.PixelIdxList{1, 1}(1);
boutEnd=whiskingEpochsList.PixelIdxList{1, 1}(end);

for frNum=1:vidObj.NumFrames %boutEnd
    imCdata = readFrame(vidObj);
%     if frNum<boutStart || frNum>boutEnd
%         continue
%     end
%     sumVidIdx=sumVidIdx+uint32(imcomplement(imCdata(:,:,1)));
    sumVidIdx=sumVidIdx+uint32((imCdata(:,:,1)));
end

sumVidIdx=uint8(255 * mat2gray(sumVidIdx));

% Display single frame stored in vidStruc
figure('color','white','position',[1274         475         613         447]);
image(sumVidIdx)
colormap jet
box off
axis off
