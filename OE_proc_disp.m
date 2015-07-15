%% Process and display data from Open-Ephys recordings  

cd C:\Data\OpenEphys\PrV25_2015-03-09_20-12-55_test2

%O_E binary format - see also 'load_open_ephys_data.m'
% fid = fopen('102_CH1.continuous'); 
fid = fopen('all_channels.events');
hdr = fread(fid, 1024, 'char*1');
timestamp = fread(fid, 1, 'int64',0,'l');
N = fread(fid, 1, 'uint16',0,'l');
recordingNumber = fread(fid, 1, 'uint16', 0, 'l');
samples = fread(fid, N, 'int16',0,'b');
recordmarker = fread(fid, 10, 'char*1');
fclose(fid);
figure; plot(samples);

%kwik data
h5disp('experiment1.kwik','/event_types/TTL')
data = h5read('experiment1.kwik','/event_types/TTL/events/recording')