function wBoutAudio=WhiskingBoutAudio(ephys,unitNum,boutIndex,factor)

%% Add spike audio
% FacePro.MakeSpikeAudio()
sampleRate=1000; %factor=20;
spikeTimes = find(ephys.rasters(unitNum,boutIndex))/sampleRate;
% waveforms = sortedSpikeWaveforms;
waveforms = repmat(mean([[gausswin(6); zeros(2,1)],[zeros(2,1);-gausswin(6)]],2), 1, sum(ephys.rasters(unitNum,boutIndex))); %gausswin(8)
spikeTimesScaled = round(spikeTimes * factor * sampleRate);
wBoutAudio = zeros(length(boutIndex)*factor, 1);
%             else
%                 spikeBinary = sweep(:);
%                 waveforms = repmat(-gausswin(8), 1, sum(spikeBinary));
%                 spikeTimesScaled = find(spikeBinary) * factor;
%                 aud = zeros(length(spikeBinary)*factor, 1);
%             end

for j = 1 : size(waveforms,2)
    wBoutAudio(spikeTimesScaled(j):spikeTimesScaled(j)+size(waveforms,1)-1) = waveforms(:,j);
end

wBoutAudio = interp1(1:length(wBoutAudio),wBoutAudio,1:sampleRate*(length(boutIndex)/1000)*factor);
wBoutAudio = reshape(wBoutAudio, [ round(length(wBoutAudio)/length(boutIndex)), length(boutIndex) ]);
% wBoutAudio = wBoutAudio(:, frameRange(1):frameRange(2));
wBoutAudio = int16(wBoutAudio / max(abs(wBoutAudio(:))) * double(intmax('int16')));