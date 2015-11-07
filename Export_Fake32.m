%% Load one channel and export to "fake" 32 channel data 

%% Get file path
[fname,dname] = uigetfile({'*.continuous;*.kwik;*.kwd;*.kwx;*.nex','All Data Formats';...
    '*.*','All Files' },'Open Ephys Recordings','C:\Data\OpenEphys\');
cd(dname);

[data, timestamps, info] = load_open_ephys_data([dname fname]);

%% Create fake 32 ch. data
convertdata=int16(data);
fake32=repmat(convertdata,32,1);
save 'fake32.mat' fake32 -V7

%% Export
cd('C:\Data\export');
expname=regexp(strrep(dname,'-','_'),'\\\w+','match');
uname=getenv('username');
expname=[expname{end}(2:end) '_' expname{end-1}(2:end) '_' uname(1)];
tic;
fileID = fopen([expname '.dat'],'w');
fprintf(fileID,'%d\n',formatdata);
fclose(fileID);
disp(['took ' num2str(toc) ' secondes to export data']);
