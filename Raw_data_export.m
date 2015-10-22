%% Process and display data from Open-Ephys recordings

%% Get file path
[fname,dname] = uigetfile('*.*','C:\Data\OpenEphys\','Open Ephys Reordings');
[data, timestamps, info] = load_open_ephys_data([dname fname]);

%% Create fake 32 ch. data
convertdata=int16(data);
% convertdata=reshape(convertdata,size(info.nsamples,2),size(data,1)/size(info.nsamples,2));
fake32=repmat(convertdata,32,1);

%% Export
cd('C:\Data\fake_data');
save 'fake32.mat' fake32 -V7
fprintf('%d\n',round(a));

fileID = fopen('fake32.dat','w');
fprintf(fileID,'%d\n',fake32);
fclose(fileID);
