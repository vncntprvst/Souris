%% Convert and export data to "raw" format

extra=0;
filterOption='difffilt'; % 'difffilt lowpass 'movav' 'multifilt';

%% Get file path
[fname,dname] = uigetfile({'*.continuous;*.kwik;*.kwd;*.kwx;*.nex','All Data Formats';...
    '*.*','All Files' },'Open Ephys Recordings','C:\Data\OpenEphys\');
cd(dname);
disp(['loading ' dname fname]);
if strfind(fname,'continuous')
    [data, timestamps, info] = load_open_ephys_data([dname fname]);
elseif strfind(fname,'raw.kwd')
    % The last number in file name from Open-Ephys recording is Node number
    % e.g., experiment1_100.raw.kwd is "raw" recording from Node 100 for
    % experiment #1 in that session.
    % Full recording parameters can be recovered from settings.xml file.
    % -<SIGNALCHAIN>
    %     -<PROCESSOR NodeId="100" insertionPoint="1" name="Sources/Rhythm FPGA">
    %         -<CHANNEL_INFO>
    %             ...
    %     -<CHANNEL name="0" number="0">
    %         <SELECTIONSTATE audio="0" record="1" param="1"/>
    %             ...
    %   ...
    %     -<PROCESSOR NodeId="105" insertionPoint="1" name="Filters/Bandpass Filter">
    %         -<CHANNEL name="0" number="0">
    %            <SELECTIONSTATE audio="0" record="1" param="1"/>
    %                <PARAMETERS shouldFilter="1" lowcut="1" highcut="600"/>
    
    %general info: h5disp(fname)
    rawInfo=h5info(fname);%'/recordings/0/data'
    rawInfo=h5info(fname,rawInfo.Groups.Name);
    %   chanInfo=h5info([regexp(fname,'^[a-z]+1','match','once') '.kwx']);
    %get basic info about recording
    rec.dur=rawInfo.Groups.Datasets.Dataspace.Size;
    rec.samplingRate=h5readatt(fname,rawInfo.Groups.Name,'sample_rate');
    rec.bitDepth=h5readatt(fname,rawInfo.Groups.Name,'bit_depth');
    %   rec.numSpikeChan= size(chanInfo.Groups.Groups,1); %number of channels with recored spikes
    rec.numRecChan=rawInfo.Groups.Datasets.Dataspace.Size-3;  %number of raw data channels.
    % Last 3 are headstage's AUX channels (e.g accelerometer)
    
    %load data (only recording channels)
    tic;
    data=h5read(fname,'/recordings/0/data',[1 1],[rec.numRecChan(1) Inf]);
    disp(['took ' num2str(toc) ' seconds to load data']);
elseif strfind(fname,'kwik')
    disp('Check out OE_proc_disp instead');
    return
elseif strfind(fname,'nex')
    disp('Only TBSI_proc_disp available right now');
    return
end

%% amplitude-dependant bandpass
if strcmp(filterOption,'lowpass')
    %% butterworth low-pass
    tic
    [b,a] = butter(3,500/rec.samplingRate,'high');
    % delay = round(max(grpdelay(b,a)));
    for chNm=1:size(data,1)
        data(chNm,:)= filter(b,a,single(data(chNm,:)));
    end
    disp(['lowpass done in ' num2str(toc) 'seconds']);
elseif strcmp(filterOption,'movav')
    %% substract moving average
    % mvaData=nans(size(data,1),size(data,2));
    avOver=5; %ms
    windowSize=avOver*(rec.samplingRate/1000);
    meanDelay = round(mean(grpdelay((1/windowSize)*ones(1,windowSize),1)));
    for chNm=1:size(data,1)
        movAverage = filter((1/windowSize)*ones(1,windowSize),1,single(data(chNm,:)));
        movAverageTail=filter((1/windowSize)*ones(1,windowSize),1,single(data(chNm,end-meanDelay+1:end)));
        movAverage = [movAverage(meanDelay+1:end) movAverageTail(avOver+2:end)] ;
        data(chNm,:)=data(chNm,:)-int16(movAverage);
    end
elseif  strcmp(filterOption,'difffilt')   
    %% multi-stage filter / amplification
    tic;
    disp('diff started');
    for chNm=1:size(data,1)
        
        disp(['Channel ' num2str(chNm) ' of ' num2str(size(data,1))]);
        
        %         figure; hold on;
        %         plot(data(chNm,:))
        %start with fairly high HP for a basis
        [b,a] = butter(3,500/rec.samplingRate,'high'); %'bandpass' is the default when Wn has two elements.
        HPfiltSample = filter(b,a,single(data(chNm,:)));
        %         plot(HPfiltSample)
        %         patch([1:length(zeros(1,rec.samplingRate)),fliplr(1:length(zeros(1,rec.samplingRate)))],...
        %         [zeros(1,rec.samplingRate)-double(fooMAD),fliplr(zeros(1,rec.samplingRate)+double(fooMAD))],'r','EdgeColor','none','FaceAlpha',0.1);
        
        % get the underlying low-frequency baseline
        [b,a] = butter(3,1000/rec.samplingRate,'low'); %'bandpass' is the default when Wn has two elements.
        filtSampleDelay= round(max(grpdelay(b,a)));
        LPfiltSample = filter(b,a,single([HPfiltSample zeros(1,filtSampleDelay)]));
        LPfiltSample = LPfiltSample(filtSampleDelay+1:length(HPfiltSample)+filtSampleDelay);
        %         plot(LPfiltSample)
        
        % compute median absolute deviation (MAD)
        sampleMAD = 4*mad(HPfiltSample,1); %sampleSTD = std(adjustSample)
        % find high amplitude signal
        %         plot(abs(HPfiltSample-LPfiltSample));
        highAmpSampleIdx=abs(HPfiltSample-LPfiltSample)>sampleMAD;
        highAmpSampleIdx([false highAmpSampleIdx(1:end-1)])=true;
        highAmpSampleIdx([highAmpSampleIdx(2:end) false])=true;
        %         plot(60*highAmpSampleIdx)
        % find high diff
        samplediff=diff(HPfiltSample); %figure; plot(samplediff)
        diffMAD=mad(samplediff);
        %         plot(abs([samplediff(1) samplediff]-LPfiltSample))
        highDiffIdx=abs([samplediff(1) samplediff]-LPfiltSample)>2*diffMAD;
        %         plot(100*highDiffIdx)
        
        % merge the two indices and enlarge detected regions by one sample on each side
        highAmpandDiffSampleIdx=highAmpSampleIdx & highDiffIdx;
        highAmpandDiffSampleIdx([false highAmpandDiffSampleIdx(1:end-1)])=true;
        highAmpandDiffSampleIdx([highAmpandDiffSampleIdx(2:end) false])=true;
        %         plot(70*highAmpandDiffSampleIdx)
        
        % select highAmpandDiffSampleIdx-identified "spikes" based on minimal size
        highAmpSample=bwlabel(highAmpSampleIdx);%identify high amplitude regions
        selectHighAmpSample=unique(highAmpSample(highAmpandDiffSampleIdx));%select id number of those that are also high diff
        selectHighAmpSample=selectHighAmpSample(selectHighAmpSample>0);% remove zero id
        highAmpSampleSize=regionprops(highAmpSampleIdx,'area');% find high amplitude segmengts' size
        highAmpSampleSize=highAmpSampleSize(selectHighAmpSample);% restrict to high diff ones
        selectHighAmpSample=selectHighAmpSample([highAmpSampleSize.Area]>mean([highAmpSampleSize.Area]));% select those that are 5 samples or more
        
        %use amplitude as template
        highAmpSampleIdx(~(highAmpandDiffSampleIdx & ismember(highAmpSample,selectHighAmpSample)))=0;
%         plot(200*highAmpSampleIdx)
        % increase chunk size
        highAmpSampleIdx([false highAmpSampleIdx(1:end-1)])=true;
        highAmpSampleIdx([highAmpSampleIdx(2:end) false])=true;
        
        %and add threshold
%         mean(HPfiltSample(HPfiltSample>1.5*std(HPfiltSample)))
     
        % lowpass data
        [b,a] = butter(3,6000/rec.samplingRate,'low');
        filtSampleDelay = round(mean(grpdelay(b,a)));
%         LPfiltSample= filter(b,a,HPfiltSample);
        LPfiltSample = filter(b,a,single([HPfiltSample zeros(1,filtSampleDelay)]));
        LPfiltSample = LPfiltSample(filtSampleDelay+1:length(HPfiltSample)+filtSampleDelay);
        LPfiltSample(highAmpSampleIdx)=LPfiltSample(highAmpSampleIdx)+HPfiltSample(highAmpSampleIdx);
        

        % save
        data(chNm,:)=LPfiltSample;
        
%         plot(1:length(HPfiltSample),ones(1,length(HPfiltSample))*2*std(HPfiltSample))
%         plot(1:length(HPfiltSample),-ones(1,length(HPfiltSample))*2*std(HPfiltSample))
        
%         plot(LPfiltSample)
        
%         BPfiltSample=BPfiltSample(filtSampleDelay+1:length(adjustSample)+filtSampleDelay);
%         data(chNm,:)=int16(BPfiltSample)+adjAmpTmplt;
        

        
        clearvars -except data extra filterOptio fname dname chNm rawInfo rec

        
        % optional: find their coordinates
%         coordPutSpikes=regionprops(highAmpandDiffSampleIdx & ismember(highAmpSample,selectHighAmpSample),'Centroid');
%         coordPutSpikes=[coordPutSpikes.Centroid];coordPutSpikes=coordPutSpikes(1:2:end);
        %         plot(coordPutSpikes,160*ones(1,length(coordPutSpikes)),'*');

        %delineate window around each "spike"
        % ideal would be 1ms, but too slow and unrealistic. 
%         So, just the time of rise and fall, as above
%         spikeTemplate=false(size(data(chNm,:)));
%         for spknm=1:size(coordPutSpikes,2)
%             tic
%             spkArea=sum(highAmpSample==highAmpSample(round(coordPutSpikes(spknm))));
%             if round(coordPutSpikes(spknm))-spkArea<1
%                 [spikeTemplate(1:round(coordPutSpikes(spknm))+rec.samplingRate/1000/2)]=deal(true);
%             elseif round(coordPutSpikes(spknm))+spkArea>size(spikeTemplate,2)
%                 [spikeTemplate(round(coordPutSpikes(spknm))-spkArea:...
%                     size(spikeTemplate,2))]=deal(true);
%             else
%                 [spikeTemplate(round(coordPutSpikes(spknm))-spkArea:...
%                     round(coordPutSpikes(spknm))+spkArea)]=deal(true);
%             end
%             toc
%         end
        
    end
        disp(['difffilt done in ' num2str(toc) 'seconds']);
elseif strcmp(filterOption,'multifilt')
    %% multi-stage filter / amplification
    tic;
    disp('multifilt started');
    for chNm=1:size(data,1)
        disp(['Channel ' num2str(chNm) ' of ' num2str(size(data,1))]);
        [b,a] = butter(3,500/rec.samplingRate,'low'); %'bandpass' is the default when Wn has two elements.
        filtSampleDelay = round(max(grpdelay(b,a)));
        LPfiltSample = filter(b,a,single([data(chNm,:) zeros(1,filtSampleDelay)]));
        LPfiltSample = LPfiltSample(filtSampleDelay+1:length(data(chNm,:))+filtSampleDelay);
        
        adjustSample=data(chNm,:)-int16(LPfiltSample);
        [b,a] = butter(3,1000/rec.samplingRate,'low'); %'bandpass' is the default when Wn has two elements.
        filtSampleDelay= round(max(grpdelay(b,a)));
        LPfiltSample = filter(b,a,single([adjustSample zeros(1,filtSampleDelay)]));
        LPfiltSample = LPfiltSample(filtSampleDelay+1:length(adjustSample)+filtSampleDelay);
        
        % compute median absolute deviation (MAD), qui vaut 1.5 STD. Donc un seuil a 7 en MAD, ca veut dire 4*STD (car 4*1.5=7).
        sampleMAD = mad(adjustSample,1); %sampleSTD = std(adjustSample)
        % find low amplitude signal
        lowampSampleIdx=abs(adjustSample-(int16(LPfiltSample)))<sampleMAD;
        %mid-amp
        midampSampleIdx=abs(adjustSample-(int16(LPfiltSample)))>=sampleMAD & abs(adjustSample-(int16(LPfiltSample)))<2*sampleMAD;
        %keep high amplitude template
        adjAmpTmplt=adjustSample;adjAmpTmplt(midampSampleIdx)=adjAmpTmplt(midampSampleIdx)/2;adjAmpTmplt(lowampSampleIdx)=0;
        %filter adjustSample
        [b,a] = butter(3,6000/rec.samplingRate,'low');
        filtSampleDelay = round(max(grpdelay(b,a)));
        BPfiltSample= filter(b,a,single([adjustSample zeros(1,filtSampleDelay)]));
        %     actualDelay=find(BPfiltSample(find(adjustSample==max(adjustSample)):end)==max(BPfiltSample(find(adjustSample==max(adjustSample)):end)));
        BPfiltSample=BPfiltSample(filtSampleDelay+1:length(adjustSample)+filtSampleDelay);
        data(chNm,:)=int16(BPfiltSample)+adjAmpTmplt;
        %         foo=int16(BPfiltSample)+adjAmpTmplt;
        %         fooMAD = 9*mad(foo,1);
        %
        %             figure; hold on;
        %             plot(data(chNm,1:rec.samplingRate));
        %             plot(LPfiltSample(1:rec.samplingRate));
        %             plot(adjustSample(1:rec.samplingRate));
        %             plot(foo(1:rec.samplingRate));
        %             plot(600*lowampSampleIdx(1:rec.samplingRate));
        %             plot(800*midampSampleIdx(1:rec.samplingRate));
        %
        %         patch([1:length(zeros(1,rec.samplingRate)),fliplr(1:length(zeros(1,rec.samplingRate)))],...
        %         [zeros(1,rec.samplingRate)-double(fooMAD),fliplr(zeros(1,rec.samplingRate)+double(fooMAD))],'r','EdgeColor','none','FaceAlpha',0.1);
        
        
        %             patch([1:length(foo(1:rec.samplingRate)),fliplr(1:length(foo(1:rec.samplingRate)))],...
        %         [foo(1:rec.samplingRate)-fooMAD,fliplr(foo(1:rec.samplingRate)+fooMAD)],'r','EdgeColor','none','FaceAlpha',0.1);
        
        
        %     [b,a] = butter(3,[600 6000]/rec.samplingRate); %'bandpass' is the default when Wn has two elements.
        %     filtSampleDelay = round(max(grpdelay(b,a)));
        %     filtfoo = filter(b,a,single([foo zeros(1,filtSampleDelay)]));
        %     filtfoo=filtfoo(filtSampleDelay+1:length(foo)+filtSampleDelay);
        %     plot(filtfoo(1:rec.samplingRate));
        
    end
    disp(['multifilt done in ' num2str(toc) 'seconds']);
    %preview
    %     [b,a] = butter(3,500/rec.samplingRate,'high'); %'bandpass' is the default when Wn has two elements.
    %     delay{1} = round(max(grpdelay(b,a)));
    %     filtData{1} = filter(b,a,single([adjustSample zeros(1,delay{1})]));
    %     filtData{1}=filtData{1}( delay{1}+1:length(adjustSample)+delay{1});
    %     figure; plot(filtData{1}(1,4400:5450));
end

%% convert data class and format
if ~isa(data,'int16')
    data=int16(data);
end
% formatdata=reshape(data',[size(data,1)*size(data,2) 1]);
%test on truncated data
% cd('C:\Data\fake_data');
% trFormatdata=data(:,1:chanInfo.Groups.Groups(1).Datasets(2).Dataspace.Size);
% expname='Trunc16ch';
% fileID = fopen([expname '.dat'],'w');
% fwrite(fileID,trFormatdata,'int16');

%% Export
cd('C:\Data\export');
expname=regexp(strrep(dname,'-','_'),'\\\w+','match');
uname=getenv('username');
expname=[expname{end}(2:end) '_' expname{end-1}(2:end) '_' uname(1)];
tic;
fileID = fopen([expname '.dat'],'w');
fwrite(fileID,data,'int16');
% fprintf(fileID,'%d\n',formatdata);
fclose(fileID);
disp(['took ' num2str(toc) ' seconds to export data']);

%% Create fake 32 ch. data
% convertdata=int16(data);
% fake32=repmat(convertdata,32,1);
% save 'fake32.mat' fake32 -V7
%%Export
% fprintf(fileID,'%d\n',formatdata);

if extra
    % Look at some data
    
    timeWindow=1:1*rec.samplingRate;
    iter=floor(size(data,2)/size(timeWindow,2)/60);
    
    period=1;
    % for period=1:2:iter
    sampleStart=(period-1)*600*rec.samplingRate;
    chan5sample=single(data(5,sampleStart+timeWindow));
    
    % substract moving average
    % mvaData=nans(size(data,1),size(data,2));
    % avOver=5; %ms
    % windowSize=avOver*(rec.samplingRate/1000);
    % meanDelay = round(mean(grpdelay((1/windowSize)*ones(1,windowSize),1)));
    % for chNm=1:size(chan5sample,1)
    %     movAverage = filter((1/windowSize)*ones(1,windowSize),1,single(chan5sample(chNm,:)));
    %     movAverageTail=filter((1/windowSize)*ones(1,windowSize),1,single(chan5sample(chNm,end-meanDelay+1:end)));
    %     movAverage = [movAverage(meanDelay+1:end) movAverageTail] ;
    %     chan5sample(chNm,:)=chan5sample(chNm,:)-(movAverage);
    % end
    
    %filters
    % [b,a] = butter(3,40/rec.samplingRate,'high'); %'bandpass' is the default when Wn has two elements.
    % filtData{1} = filter(b,a,chan5sample);
    [b,a] = butter(3,500/rec.samplingRate,'low'); %'bandpass' is the default when Wn has two elements.
    delay{1} = round(max(grpdelay(b,a)));
    filtData{1} = filter(b,a,[chan5sample zeros(1,delay{1})]);
    [b,a] = butter(3,500/rec.samplingRate,'high'); %'bandpass' is the default when Wn has two elements.
    delay{2} = round(max(grpdelay(b,a)));
    filtData{2} = filter(b,a,[chan5sample zeros(1,delay{2})]);
    [b,a] = butter(3,[600 6000]/rec.samplingRate); %'bandpass' is the default when Wn has two elements.
    delay{3} = round(max(grpdelay(b,a)));
    filtData{3} = filter(b,a,[chan5sample zeros(1,delay{3})]);
    % [b,a] = butter(3,[2800 3750]/rec.samplingRate,'stop');
    % filtData{3} = filter(b,a,filtData{3});
    
    % substract moving average
    % mvaData=nans(size(data,1),size(data,2));
    % avOver=5; %ms
    % windowSize=avOver*(rec.samplingRate/1000);
    % meanDelay = round(mean(grpdelay((1/windowSize)*ones(1,windowSize),1)));
    % for chNm=1:size(filtData{2},1)
    %     movAverage = filter((1/windowSize)*ones(1,windowSize),1,single(filtData{2}(chNm,:)));
    %     movAverageTail=filter((1/windowSize)*ones(1,windowSize),1,single(filtData{2}(chNm,end-meanDelay+1:end)));
    %     movAverage = [movAverage(meanDelay+1:end) movAverageTail] ;
    %     filtData{2}(chNm,:)=filtData{2}(chNm,:)-(movAverage);
    % end
    
    figure;
    %plot
    % subplot(floor(iter/2),2,period);
    plot(chan5sample(4400:5450));
    hold on
    plot(filtData{1}(4400+delay{1}-10:5450+delay{1}-10));
    grid on
    axis tight
    % set(gca,'ylim',[-5000 5000]);
    
    sample=chan5sample(4400:5450)-filtData{1}(4400+delay{1}-10:5450+delay{1}-10);
    [b,a] = butter(3,1000/rec.samplingRate,'low'); %'bandpass' is the default when Wn has two elements.
    filtSampleDelay= round(max(grpdelay(b,a)));
    LPfiltSample = filter(b,a,[sample zeros(1,filtSampleDelay)]);
    LPfiltSample = LPfiltSample(filtSampleDelay+1:end);
    
    % compute median absolute deviation (MAD), qui vaut 1.5 STD. Donc un seuil a 7 en MAD, ca veut dire 4*STD (car 4*1.5=7).
    sampleMAD = mad(sample,1); %sampleSTD = std(sample)
    % find low amplitude signal
    lowampSampleIdx=abs(sample-(LPfiltSample))<sampleMAD;
    %keep high amplitude template
    adjAmpTmplt=sample;adjAmpTmplt(lowampSampleIdx)=0;
    %filter sample (should just take filtData{3} eventually)
    [b,a] = butter(3,[600 6000]/rec.samplingRate);
    filtSampleDelay = round(max(grpdelay(b,a)));
    BPfiltSample= filter(b,a,[sample zeros(1,filtSampleDelay)]);
    actualDelay=find(BPfiltSample(find(sample==max(sample)):end)==max(BPfiltSample(find(sample==max(sample)):end)));
    BPfiltSample=BPfiltSample(actualDelay+1:length(sample)+actualDelay);
    adjustSample=BPfiltSample+adjAmpTmplt;
    % adjustSample=BPfiltSample-LPfiltSample+highampTemp;
    
    figure;
    plot(sample);
    hold on
    plot(LPfiltSample)
    patch([1:length(LPfiltSample),fliplr(1:length(LPfiltSample))],...
        [LPfiltSample-sampleMAD,fliplr(LPfiltSample+sampleMAD)],'r','EdgeColor','none','FaceAlpha',0.1);
    % plot(lowampSampleIdx*600)
    plot(BPfiltSample)
    plot(adjustSample)
    % plot(filtData{2}(4400:5450));
    % plot(filtData{3}(4400:5450));
    legend('sample','LPfiltSample','sampleMAD','BPfiltSample','adjustSample')
    
    %power spectrum
    % foo=chan5sample(4400:5450)-filtData{1}(4400+delay{1}-10:5450+delay{1}-10);
    % foopws=spa(foo);
    % foop=etfe(foo);
    % spectrum(foopws,foop);
    
    % fft
    figure;  plot(filtData{3}(5400:5450));
    NFFT = 2^nextpow2(size(filtData{3}(5400:5450),2)); % Next power of 2 from signal length
    Y = fft(filtData{3}(5400:5450),NFFT)/size(filtData{3}(5400:5450),2);
    f = rec.samplingRate/2*linspace(0,1,NFFT/2+1);
    figure; plot(f,2*abs(Y(1:NFFT/2+1)))
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    % subplot(211)
    % subplot(floor(iter/2),2,period+1);
    % plot(filtData{1});
    % subplot(212)
    plot(filtData{2});hold on;
    plot(filtData{3})
    hold off
    grid on
    axis tight
    set(gca,'ylim',[-2000 2000])
    set(gca,'xtick',linspace(0,max(timeWindow),100),...
        'xticklabel',linspace(0,max(timeWindow)/rec.samplingRate*1000,numel(linspace(0,max(timeWindow),100))))
    legend('filt1','filt2','filt3')
    % zoom on
    % end
    
    % bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    %     'HalfPowerFrequency1',300,'HalfPowerFrequency2',double(rec.samplingRate/2), ...
    %     'SampleRate',rec.samplingRate); %4th-order Butterworth bandpass
    % grpdelay(bpFilt,2048,rec.samplingRate)   % plot group delay
    % fvtool(bpFilt)
    % delay = round(mean(grpdelay(df))); % filter delay in samples
    %
    % filtData = filter(df,[chan5sample; zeros(delay,1)]); % Append D zeros to the input data
    % filtData = filtData(delay+1:end);                  % Shift data to compensate for delay
end
