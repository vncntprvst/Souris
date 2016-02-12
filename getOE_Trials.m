function TTL=getOE_Trials(fileName)
%% get Trial structure from TTLs

switch nargin
    case 0
        [fileName,dname] = uigetfile({'*.kwe','.Kwik Files';...
            '*.events','OE format Files'; '*.*','All Files'},'Events Data');
        cd(dname)
    case 1
        % fine
end

if strfind(fileName,'events')
   [~, TTL_times, info] = load_open_ephys_data(fileName);
   TTLevents=info.eventType==3;
   TTL_ID=info.eventId(TTLevents);
   TTL_times=TTL_times(TTLevents).*1000; %convert to ms scale
else
    % h5disp('experiment1.kwe','/event_types/TTL')
    % TTLinfo=h5info('experiment1.kwe','/event_types/TTL');
    TTL_ID = h5read(fileName,'/event_types/TTL/events/user_data/eventID');
    TTL_times = h5read(fileName,'/event_types/TTL/events/time_samples');
    sampleRate = h5readatt(fileName,'/recordings/0/','sample_rate');
    TTL_times = TTL_times./uint64(sampleRate/1000); %convert to ms scale
end

% keep absolute time of TTL onset
TTL_times=TTL_times(diff([0;TTL_ID])>0);
% TTL sequence
TTL_seq=diff(TTL_times);

    % sending 10ms TTLs, with 10ms interval. Two pulses for begining of trial
    %(e.g.,head through front panel). With sampling rate of 30kHz, that
    % interval should be 601 samples (20ms*30+1). Or 602 accounting for jitter.
    % TTL_seq at native sampling rate should thus read as:
    %   601
    %   end of trial time
    %   inter-trial interval
    %   601 ... etc

TTLlength=mode(TTL_seq)/2; %might want to code a better TTL length detector

if TTL_seq(1)>=TTLlength*2+10 %missed first trial initiation, discard times
    TTL_seq(1)=TTLlength+300;
end
if TTL_seq(end)<=TTLlength*2+10 %unfinished last trial
    TTL_seq(end)=TTLlength+300;
end

TTL.start=TTL_times([TTL_seq<=TTLlength*2+10;false]);%Trials.start=Trials.start./uint64(SamplingRate/1000)
TTL.end=TTL_times(find([TTL_seq<=TTLlength*2+10;false])+2);
try
    TTL.interval=TTL_times(find([TTL_seq(1:end-3)<=TTLlength+10;false])+3)-TTL_times(find([TTL_seq(1:end-3)<=TTLlength+10;false])+2);
catch
    TTL.interval=[]; %
end
end