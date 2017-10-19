%%% select video file %%%
[file_name, path_name] = uigetfile('*.avi', 'Select coordinates file', 'MultiSelect', 'on');
if ~iscell(file_name)
    file_name = {file_name};
end

% plot first frame
vid = VideoReader([path_name, file_name{1}]);
firstFrame=readFrame(vid);
figure('position',[1067 217 846 724]); subplot(2,1,1); imagesc(firstFrame);

isRotationGood='No';
while strcmp(isRotationGood,'No')
    %find rotation
    angleRectify=str2double(cell2mat(inputdlg('Enter rotation angle to have head pointing up:',...
        'Rotation', 1,{'15'})));
    rotatedFirstFrame=imrotate(firstFrame, angleRectify);
    subplot(2,1,2); imagesc(rotatedFirstFrame);
    
    isRotationGood = questdlg('Is the angle correct?', ...
        '', 'Yes','No','Yes');
end
close(gcf)

%find whisker region
isCropGood='No';
while strcmp(isCropGood,'No')
    figure('position',[1193 575 560 420]); imagesc(rotatedFirstFrame);
    rect = getrect; %returns [xmin ymin width height]
    whiskerRegion= round([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
    whiskerRegionImage=rotatedFirstFrame(whiskerRegion(3):whiskerRegion(4),whiskerRegion(1): whiskerRegion(2),:);
    figure('position',[1178 63 560 420]); imagesc(whiskerRegionImage);
    isCropGood = questdlg('Is the cropping correct?', ...
            '', 'Yes','No','Yes');
    close all
end

%find nose region
isCropGood='No';
while strcmp(isCropGood,'No')
    figure('position',[1193 575 560 420]); imagesc(whiskerRegionImage);
    rectn = getrect; %returns [xmin ymin width height]
    noseRegion= round([rectn(1) rectn(1)+rectn(3) rectn(2) rectn(2)+rectn(4)]);
    noseRegionImage=whiskerRegionImage(noseRegion(3):noseRegion(4),noseRegion(1): noseRegion(2),:);
    figure('position',[1178 63 560 420]); imagesc(noseRegionImage);
    isCropGood = questdlg('Is the cropping correct?', ...
        '', 'Yes','No','Yes');
    close all
end

[thetas,nosem]=GetWhiskerAngle(vid,angleRectify,whiskerRegion,noseRegion);

%%% save: you can modify the path %%%
cd('C:\Data\Ephys\Behav');
save([ file_name{1}(1:end-3) 'mat'], 'thetas', 'nosem');


