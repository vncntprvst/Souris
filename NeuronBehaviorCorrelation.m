%% Correlation between bursts/spike rate and periodical behaviors (whisking, breathing)

% load periodical behavior data
[filename,dname] = uigetfile({'*.csv','.csv Files';...
    '*.*','All Files' },'TTL Onset Data','C:\Data\Behav');
cd(dname)
fileID = fopen(filename,'r');

