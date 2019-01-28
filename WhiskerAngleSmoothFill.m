function thetas=WhiskerAngleSmoothFill(varargin)
% convert values to degrees
if nargin==1 %orientation values in radian
    thetas=thetas*180/pi; %equivalent to rad2deg    
%center to zero median
%   thetas=thetas-nanmedian(thetas);
elseif nargin==2 %centroid values
    centroidX=varargin{1}; centroidY=varargin{2};
    %   rebase to whisker referential (mid-point termined from maxima)
    centroidX=centroidX-mean([min(centroidX) max(centroidX)]);
	centroidY=centroidY-mean([min(centroidY) max(centroidY)]);
    centroidY=centroidY-min(centroidY);
    %   convert to radian 
    centroidX=centroidX/max(abs(centroidX))*pi/2;
    centroidY=centroidY/max(abs(centroidY))*pi/2;
    thetas = FindAngle(centroidX,centroidY);
end


% find outliers
outliers=abs(thetas)>30*mad(abs(thetas));
% could use filloutliers
thetas(outliers)=nan;

% plots
% figure; hold on; plot(thetas)
% plot(outliers,zeros(numel(outliers,1)),'rd')

% find missing bits length
% nanBouts=hist(bwlabel(isnan(thetas)),unique(bwlabel(isnan(thetas))));
% nanBouts=sort(nanBouts,'descend'); maxNanBouts=nanBouts(2);

% smooth
thetas=smoothdata(thetas,'rloess',7);
%linear %spline %movmedian
%     [b,a] = butter(3,40/250,'low');
%     foo = filtfilt(b,a,thetas);

%fill missing / NaNs values (if any)
fillDim=find(size(thetas)==max(size(thetas)));
thetas=fillmissing(thetas,'spline',fillDim,'EndValues','nearest');
% plot(thetas)

function wAngle = FindAngle(centroidX,centroidY)
for i = 1:length(centroidX)
    if atan2(centroidY(i),centroidX(i))>=0
        wAngle(i) = (180/pi) * (atan2(centroidY(i),centroidX(i)));
    else
        wAngle(i) = (180/pi) * (atan2(centroidY(i),centroidX(i))+2*pi);
    end
end