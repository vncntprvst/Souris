function spikes=Spike2_Export(dirName,fileName)
% dirName='D:\Vincent\vIRt57\vIRt57_0216';
%'D:\Vincent\vIRt61\vIRt61_0302'

%'D:\Vincent\sc001\sc001_072221';
% fileName='vIRt57_0216_5732.ns6';
%'vIRt61_0302_5631.ns6'
%'vIRt61_0302_5926.ns6'
%'sc001_072221_run001.ns6';

% clearvars

data=LoadEphys_Blackrock(dirName,fileName);

preprocOption={'CAR','all'};
samplingRate=30000;
allTraces=PreProcData(data,samplingRate,preprocOption);

%% save channel of interest
chNum=3;
exportDir=fullfile(dirName,'SpikeSorting',fileName(1:end-4));
cd(exportDir)
if ~exist('Spike2Export','dir'); mkdir('Spike2Export'); end

fileID = fopen(fullfile(exportDir,'Spike2Export',...
    [fileName(1:end-4) '_' num2str(chNum) '.bin']),'w');
fwrite(fileID,allTraces(chNum,:),'int16');
fclose(fileID);

%% Load data
fName=['vIRt57_0216_5732_Ch' num2str(chNum) '.mat'];
spikes=LoadSpikes_Spike2(fName);

if false
    figure; plot([1:size(allTraces,2)]/30000, allTraces(chNum,:)) % 10*10^6:11*10^6)) %1:300000))
    for chNum=1:32
        figure; plot(allTraces(chNum,1:300000))
    end
    save(fullfile(dirName, 'allTraces.mat'), 'allTraces')
end
end