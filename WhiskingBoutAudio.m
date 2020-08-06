function wBoutAudio=WhiskingBoutAudio(rasters,sampleRate,factor)

%% Add spike audio
% FacePro.MakeSpikeAudio()
% sampleRate=1000; %factor=20;
spikeTimes = find(rasters)/sampleRate;
% waveforms = sortedSpikeWaveforms;
adjFact=factor/10;
% waveforms = repmat(mean([[gausswin(6); zeros(2,1)],...
%     [zeros(2,1);-gausswin(6)]],2), 1, sum(rasters(:,boutIndex))); %gausswin(8)
% waveforms = repmat(mean([[gausswin(6*adjFact); zeros(2*adjFact,1)],...
%     [zeros(2*adjFact,1);-gausswin(6*adjFact)]],2), 1, sum(rasters(:,boutIndex))); %gausswin(8)
waveforms = repmat(gausswin(8*adjFact), 1, sum(rasters)); 

spikeTimesScaled = round(spikeTimes * factor * sampleRate);
wBoutAudio = zeros(numel(rasters)*factor, 1);
%             else
%                 spikeBinary = sweep(:);
%                 waveforms = repmat(-gausswin(8), 1, sum(spikeBinary));
%                 spikeTimesScaled = find(spikeBinary) * factor;
%                 aud = zeros(length(spikeBinary)*factor, 1);
%             end

for spikeNum = 1 : size(waveforms,2)
    wBoutAudio(spikeTimesScaled(spikeNum):spikeTimesScaled(spikeNum)+size(waveforms,1)-1) = ...
        wBoutAudio(spikeTimesScaled(spikeNum):spikeTimesScaled(spikeNum)+size(waveforms,1)-1) + waveforms(:,spikeNum);
end

% wBoutAudio = interp1(1:length(wBoutAudio),wBoutAudio,1:sampleRate*(numel(rasters)/sampleRate)*factor);

wBoutAudio = reshape(wBoutAudio, [ round(length(wBoutAudio)/numel(rasters)), numel(rasters) ]);

% wBoutAudio = wBoutAudio(:, frameRange(1):frameRange(2));
wBoutAudio = int16(wBoutAudio / max(abs(wBoutAudio(:))) * double(intmax('int16')));
% 
%             aud = interp1(1:length(aud),aud,1:sampleRate*5*factor);
%             aud = reshape(aud, [ round(length(aud)/2500), 2500 ]);
%             aud = aud(:, frameRange(1):frameRange(2));
%             aud = int16(aud / max(abs(aud(:))) * double(intmax('int16')));
% 
% 
%             aud = FacePro.MakeSpikeAudio(rasters(boutIndex), 10, [1 2500], sampleRate);
%             
% %             
% soundsc(wBoutAudio',sampleRate*factor/20)
% soundsc(double(vertcat(wBoutAudio(:))));
% soundsc(waveforms(:,1),sampleRate/2*10) 
% soundsc(gausswin(8),sampleRate/2*10) 
% soundsc(double(vertcat(aud(:))))
% 
% foo=rasters(unitNum,boutIndex);
% figure; hold on;
% imagesc(wBoutAudio)
% % plot(rasters(unitNum,boutIndex))
% plot(spikeTimes*sampleRate+0.5,ones(numel(spikeTimes),1),'kd')
% 
% figure; imagesc(waveforms)
% figure; imagesc(wBoutAudio)
% figure; imagesc(aud)

