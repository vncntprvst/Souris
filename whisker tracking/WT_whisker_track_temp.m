%%% select video file %%%
[file_name, path_name] = uigetfile('*.avi', 'Select coordinates file', 'MultiSelect', 'on');
if ~iscell(file_name)
    file_name = {file_name};
end

%
vid = VideoReader([path_name, file_name{1}]);
firstFrame=readFrame(vid);
figure; imagesc(firstFrame);

%find rotation 
angleRectify=30;
rotatedFirstFrame=imrotate(firstFrame, angleRectify);
figure; imagesc(rotatedFirstFrame);

%find whisker region
rect = getrect; %returns [xmin ymin width height]
whiskerRegion= round([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);


[thetas,nosem]=GetWhiskerAngle(vid,angleRectify,whiskerRegion,noseRegion);


%%% save: you can modify the path %%%
save([ file_name{1}(1:end-3) 'mat'], 'thetas', 'nosem');







%% demo of Hough transform %%
figure
clf
imagesc(img)
hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

