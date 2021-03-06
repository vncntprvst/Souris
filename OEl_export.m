% Export function for analysis in OpenElectrophy

function OEl_export(data)
%% Data structure:
% see   http://pythonhosted.org/neo/core.html
%       http://pythonhosted.org/neo/core.html#grouping-objects
%       http://pythonhosted.org/neo/core.html#neo-diagram
% Block:
%     The top-level container gathering all of the data, discrete and
%     continuous, for a given recording session. Contains Segment and
%     RecordingChannelGroup objects.

% Segment:
%     A container for heterogeneous discrete or continous data sharing a
%     common clock (time basis) but not necessarily the same sampling rate,
%     start time or end time. A Segment can be considered as equivalent to
%     a �trial�, �episode�, �run�, �recording�, etc., depending on the
%     experimental context. May contain any of the data objects.


% Data objects
% 
    % These objects directly represent data as arrays of numerical values
    % with associated metadata (units, sampling frequency, etc.).
    % 
    % AnalogSignal:
    %     A regular sampling of a continuous, analog signal.
    % AnalogSignalArray:
    %     A regular sampling of a multichannel continuous analog signal.
    %     This representation (as a 2D NumPy array) may be more efficient
    %     for subsequent analysis than the equivalent list of individual
    %     AnalogSignal objects.
    % Spike:
    %     One action potential characterized by its time and waveform.
    % SpikeTrain:
    %     A set of action potentials (spikes) emitted by the same unit in a
    %     period of time (with optional waveforms).
    % Event and EventArray:
    %     A time point representng an event in the data, or an array of
    %     such time points.
    % Epoch and EpochArray:
    %     An interval of time representing a period of time in the data, or
    %     an array of such intervals.

%% Preallocate
block = struct();
block.segments = { };
block.name = 'my block with matlab';
for s = 1:3
    seg = struct();
    seg.name = strcat('segment ',num2str(s));
    seg.analogsignals = { };
    for a = 1:5
        anasig = struct();
        anasig.array = rand(100,1);
        anasig.units = 'mV';
        anasig.t_start = 0;
        anasig.t_start_units = 's';
        anasig.sampling_rate = 100;
        anasig.sampling_rate_units = 'Hz';
        seg.analogsignals{a} = anasig;
    end
    seg.spiketrains = { };
    for t = 1:7
        sptr = struct();
        sptr.array = rand(30,1)*10;
        sptr.units = 'ms';
        sptr.t_start = 0;
        sptr.t_start_units = 'ms';
        sptr.t_stop = 10;
        sptr.t_stop_units = 'ms';
        seg.spiketrains{t} = sptr;
    end

    block.segments{s} = seg;
end

%% Export
save 'myblock.mat' block -V7

end