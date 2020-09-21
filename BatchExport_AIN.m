function [dataFiles,allRecInfo]=BatchExport_AIN(exportDir)
% not finalized yet
% if export for SC, must not be run as root (so start Matlab from /bin/matlab as user)
% Vincent Prevosto 10/16/2018
rootDir=cd;
if ~isfolder('SpikeSorting')
    %create export directory
    mkdir('SpikeSorting');
end
dataFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*raw.kwd','*RAW*Ch*.nex','*.ns6'},'UniformOutput', false);
dataFiles=vertcat(dataFiles{~cellfun('isempty',dataFiles)});
% just in case other export / spike sorting has been performed, do not include those files
dataFiles=dataFiles(~cellfun(@(flnm) contains(flnm,{'_export';'_TTLs'; '_trialTTLs'; '_vSyncTTLs';...
    'temp_wh';'_nopp.dat';'_all_sc';'_VideoFrameTimes';'_Wheel'}),...
    {dataFiles.name})); %by filename
dataFiles=dataFiles(~cellfun(@(flnm) contains(flnm,{'_SC';'_JR';'_ML'}),...
    {dataFiles.folder})); % by folder name
if ~exist('exportDir','var')
    exportDir=(fullfile(rootDir,'SpikeSorting'));
end

%% export each file
for fileNum=1:size(dataFiles,1)
    
    %% load other data
    NEVdata=openNEV(fullfile(dataFiles(fileNum).folder,[dataFiles(fileNum).name(1:end-3), 'nev']));
    fsIdx=cellfun(@(x) contains(x','FlowSensor'),{NEVdata.ElectrodesInfo.ElectrodeLabel});
    if any(fsIdx)
        try
            fsData = openNSx(fullfile(dataFiles(fileNum).folder,[dataFiles(fileNum).name(1:end-3), 'ns4']));
            fsIdx = cellfun(@(x) contains(x,'FlowSensor'),{fsData.ElectrodesInfo.Label});
            fsData = fsData.Data(fsIdx,:);
            % downsample to 1kHz
            fsData = decimate(double(fsData),10);
        catch
            fsData = openNSx(fullfile(dataFiles(fileNum).folder,[dataFiles(fileNum).name(1:end-3), 'ns2']));
            fsIdx = cellfun(@(x) contains(x,'FlowSensor'),{fsData.ElectrodesInfo.Label});
            fsData = fsData.Data(fsIdx,:);
        end
    end
    reIdx = cellfun(@(x) contains(x','RotaryEncoder'),{NEVdata.ElectrodesInfo.ElectrodeLabel});
    if any(reIdx)
        reData = openNSx(fullfile(dataFiles(fileNum).folder,[dataFiles(fileNum).name(1:end-3), 'ns2']));
        reIdx = cellfun(@(x) contains(x,'RotaryEncoder'),{reData.ElectrodesInfo.Label});
        reData = reData.Data(reIdx,:);
    end
    
    vSyncTTLDir=cd;
    %% get recording name
    % (in case they're called 'continuous' or some bland thing like this)
    % basically, Open Ephys
    if contains(dataFiles(fileNum).name,'continuous')
        foldersList=regexp(strrep(dataFiles(fileNum).folder,'-','_'),...
            ['(?<=\' filesep ').+?(?=\' filesep ')'],'match');
        expNum=foldersList{cellfun(@(fl) contains(fl,'experiment'),foldersList)}(end);
        recNum=foldersList{cellfun(@(fl) contains(fl,'recording'),foldersList)}(end);
        recordingName=foldersList{find(cellfun(@(fl) contains(fl,'experiment'),foldersList))-1};
        recordingName=[recordingName '_' expNum '_' recNum];
    elseif contains(dataFiles(fileNum).name,'experiment')
        folderIdx=regexp(dataFiles(fileNum).folder,['(?<=\w\' filesep ').+?']);
        if isempty(folderIdx)
            folderIdx=1;
        end
        recordingName=strrep(dataFiles(fileNum).folder(folderIdx(end):end),'-','_');
    else
        recordingName=dataFiles(fileNum).name(1:end-4);
    end
    
    % collect info
%     recInfo.dataPoints=int32(recInfo.dataPoints);
%     recInfo.baseName=recordingName;
%     recNameComp=regexp(strrep(recordingName,'_','-'),'\w+','match');
%     recInfo.subject=recNameComp{1};
%     recInfo.shortDate=recNameComp{2};
%     recInfo.probeDepth=recNameComp{3};
% 
    %% check that recordingName doesn't have special characters
    recordingName=regexprep(recordingName,'\W','');
%     allRecInfo{fileNum}.recordingName=recordingName;
    
    cd(exportDir)
    if ~isfolder(recordingName)
        %create export directory
        mkdir(recordingName);
    end
    cd(recordingName)
%     recInfo.export.directory=fullfile(exportDir,recordingName);
        
    %% save other data
    if exist('fsData','var')
        if ~isfolder(fullfile(dataFiles(fileNum).folder,'FlowSensor'))
            mkdir(fullfile(dataFiles(fileNum).folder,'FlowSensor'));
        end
        fileID = fopen(fullfile(dataFiles(fileNum).folder,'FlowSensor',[recordingName '_fs.bin']),'w');
        fwrite(fileID,fsData,'int16');
        fclose(fileID);
        recInfo.export.binFile=[recordingName '_fs.bin'];
    end
    if exist('reData','var')
        if ~isfolder(fullfile(dataFiles(fileNum).folder,'RotaryEncoder'))
            mkdir(fullfile(dataFiles(fileNum).folder,'RotaryEncoder'));
        end
        fileID = fopen(fullfile(dataFiles(fileNum).folder,'RotaryEncoder',[recordingName '_re.bin']),'w');
        fwrite(fileID,reData,'int16');
        fclose(fileID);
        recInfo.export.binFile=[recordingName '_re.bin'];
    end
        
end
cd ..


