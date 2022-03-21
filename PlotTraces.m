clear all
dirName='D:\Vincent\vIRt57\vIRt57_0216'; 
%'D:\Vincent\vIRt61\vIRt61_0302'
% 
%'D:\Vincent\sc001\sc001_072221';
fileName='vIRt57_0216_5732.ns6'; 
%'vIRt61_0302_5631.ns6'
% 'vIRt61_0302_5926.ns6'
% 
%'sc001_072221_run001.ns6';

[data,rec,spikes,TTLs]=LoadEphys_Blackrock(dirName,fileName);

preprocOption={'CAR','all'};
samplingRate=30000;
allTraces=PreProcData(data,samplingRate,preprocOption);

figure; plot([1:size(allTraces,2)]/30000, allTraces(27,:)) % 10*10^6:11*10^6)) %1:300000))

for chNum=1:32
    figure; plot(allTraces(chNum,1:300000))
end
% save(fullfile(dirName, 'allTraces.mat'), 'allTraces')

