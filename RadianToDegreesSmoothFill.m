function thetas=RadianToDegreesSmoothFill(thetas)
% thetas=(multiWhiskerTrackingData(:,7));
% convert to degrees
thetas=thetas*180/pi;

%center to zero median
thetas=thetas-nanmedian(thetas);

% find outliers
outliers=abs(thetas)>30*mad(abs(thetas));
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
thetas=fillmissing(thetas,'spline',1,'EndValues','nearest');
% plot(thetas)