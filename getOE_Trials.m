function TTL=getOE_Trials(fileName)
        %% get Trial structure from TTLs  
        % h5disp('experiment1.kwe','/event_types/TTL')
        % TTLinfo=h5info('experiment1.kwe','/event_types/TTL');
        
        TTL_ID = h5read('experiment1.kwe','/event_types/TTL/events/user_data/eventID');
        TTL_times = h5read('experiment1.kwe','/event_types/TTL/events/time_samples');
        
        if mode(diff(TTL_times)/30000)==1 %assuming 30k sampling rate
            % 1 second TTL pulses
            TTLlength=30000; % actually closer to 29979 but that's fine
        else
        %we send 10ms TTLs, with 10ms interval. Two pulses for begining of trial
        %(e.g.,head through front panel). With sampling rate of 30kHz, that
        % interval should be 601 samples (20ms*30+1). Or 602 accounting for jitter.
        % TTL_seq should thus read as:
        %   601
        %   end of trial time
        %   inter-trial interval
            TTLlength=600; %might want to code a better TTL length detector
        end
        
        
        % keep absolute time of TTL onset
        TTL_times=TTL_times(diff([0;TTL_ID])>0);
        % TTL sequence
        TTL_seq=diff(TTL_times);
        
        if mode(TTL_seq/30000)==1 %assuming 30k sampling rate
            % 1 second TTL pulses
            TTLlength=30000; % actually closer to 29979 but that's fine
        else
        %we send 10ms TTLs, with 10ms interval. Two pulses for begining of trial
        %(e.g.,head through front panel). With sampling rate of 30kHz, that
        % interval should be 601 samples (20ms*30+1). Or 602 accounting for jitter.
        % TTL_seq should thus read as:
        %   601
        %   end of trial time
        %   inter-trial interval
            TTLlength=600; %might want to code a better TTL length detector
        end
        
        if TTL_seq(1)>=TTLlength+10 %missed first trial initiation, discard times
            TTL_seq(1)=TTLlength+300;
        end
        if TTL_seq(end)<=TTLlength+10 %unfinished last trial
            TTL_seq(end)=TTLlength+300;
        end
        
        TTL.start=TTL_times([TTL_seq<=TTLlength+10;false]);%Trials.start=Trials.start./uint64(SamplingRate/1000)
        TTL.end=TTL_times(find([TTL_seq<=TTLlength+10;false])+2);
        try
            TTL.interval=TTL_times(find([TTL_seq(1:end-3)<=TTLlength+10;false])+3)-TTL_times(find([TTL_seq(1:end-3)<=TTLlength+10;false])+2);
        catch
            TTL.interval=[]; %
        end
end