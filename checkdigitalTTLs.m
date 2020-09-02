fname='vIRt47_0803_4651.nev'; %testAMS
dname=cd; %'D:\Data\test\';
% [rec,data,spikes,TTLs] = LoadEphysData(fname,dname);

NEVdata=openNEV(fullfile(dname,fname), 'read');

% unique(NEVdata.Data.SerialDigitalIO.UnparsedData  )
% sum(NEVdata.Data.SerialDigitalIO.UnparsedData==3)

% try change instead of rising edge

digInEvents=NEVdata.Data.SerialDigitalIO.UnparsedData;
digInTimes=NEVdata.Data.SerialDigitalIO.TimeStamp; %TimeStampSec i interval in ms? 

% Given 2 inputs in Port 0 and 1
% No input              => bin2dec('00') = 0
% input on Port 0 only  => bin2dec('10') = 2
% input on Port 1 only  => bin2dec('01') = 1
% input on both ports   => bin2dec('11') = 3

TTL_ID=logical(dec2bin(digInEvents)-'0');
max(bwlabel(TTL_ID(:,end)))
max(bwlabel(TTL_ID(:,end-1)))

% vSyncTTLTimes=digInTimes(TTL_ID(:,1));
% figure;
% plot(diff(vSyncTTLTimes))

videoSyncTTLs_1=regionprops(TTL_ID(:,1),'PixelIdxList');
videoSyncTTLs_2=regionprops(TTL_ID(:,2),'PixelIdxList');

vSyncTTLTimes_1=double(digInTimes(cellfun(@(timeIndex) timeIndex(1),...
    {videoSyncTTLs_1.PixelIdxList})))/30;
vSyncTTLTimes_2=double(digInTimes(cellfun(@(timeIndex) timeIndex(1),...
    {videoSyncTTLs_2.PixelIdxList})))/30;

figure;
plot(diff(vSyncTTLTimes_1))

figure;
plot(vSyncTTLTimes_2-vSyncTTLTimes_1(1:numel(vSyncTTLTimes_2)))

% videoSyncTTLsProperties=regionprops(TTL_ID(:,1),'PixelIdxList');

% Compare with video file
videoFiles = cellfun(@(fileFormat) dir([fileFormat]),...
    {'*.mp4','*.avi'},'UniformOutput', false);
videoFiles=vertcat(videoFiles{~cellfun('isempty',videoFiles)});
videoFileName=videoFiles(1).name;
videoData = py.cv2.VideoCapture(videoFileName);
numFrames=videoData.get(py.cv2.CAP_PROP_FRAME_COUNT);