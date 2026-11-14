function [onsetTimes, onsetInfo] = onset_detection(audioFile, showPlot)
%ONSET_DETECTION Spectral-based novelty onset detection (Muller, FMP Sec. 6.1.2)
%   onsetTimes = onset_detection(audioFile)
%   [onsetTimes, onsetInfo] = onset_detection(audioFile, showPlot)
%
%   audioFile  - path to a WAV file
%   showPlot   - true to show the diagnostic figure (default false)
%   onsetTimes - column vector of onset times (s)
%   onsetInfo  - struct with fs, frame times, flux curves, and parameters

if nargin < 2
    showPlot = false;
end
if ~isfile(audioFile)
    error('onset_detection:fileNotFound', 'File not found: %s', audioFile);
end

winLen     = 2048;    % N (samples)
hopLen     = 512;     % H (samples)
gammaComp  = 100;     % log compression factor
avgHalfWin = 10;      % M (frames)
peakThresh = 0.2;     % min peak height (normalized flux)
minSepSec  = 0.1;     % min time between onsets (s)

%% Load audio
[audioRaw, fs] = audioread(audioFile);
audioMono = mean(audioRaw, 2);
audioMono = audioMono / max(abs(audioMono));

% Zero-pad so the first note has a silent frame before it
audioPad = [zeros(winLen,1); audioMono; zeros(winLen,1)];

%% STFT
hannWin   = 0.5 - 0.5*cos(2*pi*(0:winLen-1)'/winLen);
numFrames = floor((length(audioPad) - winLen)/hopLen) + 1;
numBins   = winLen/2 + 1;
stftMat   = zeros(numBins, numFrames);
for m = 1:numFrames
    segFFT = fft(audioPad((m-1)*hopLen + (1:winLen)) .* hannWin);
    stftMat(:, m) = segFFT(1:numBins);
end
frameTime = ((0:numFrames-1)*hopLen + winLen/2 - winLen) / fs;

%% Log compression:  Y = log(1 + gamma*|X|)
logSpec = log(1 + gammaComp*abs(stftMat));

%% Spectral flux:  sum_k H( Y(n+1,k) - Y(n,k) ),  assigned to frame n+1
specFlux = [0, sum(max(diff(logSpec, 1, 2), 0), 1)];

%% Local average and enhanced flux
localAvg = movmean(specFlux, 2*avgHalfWin + 1);
enhFlux  = max(specFlux - localAvg, 0);
enhFlux  = enhFlux / max(enhFlux);

%% Peak picking
[pkHeights, peakTimes] = findpeaks(enhFlux, frameTime, ...
    'MinPeakHeight', peakThresh, 'MinPeakDistance', minSepSec);
onsetTimes = max(peakTimes(:), 0);   % frame-center bias can put the first onset slightly < 0

%% Outputs
onsetInfo = struct('fs', fs, 'frameTime', frameTime, 'specFlux', specFlux, ...
    'localAvg', localAvg, 'enhFlux', enhFlux, 'peakHeights', pkHeights(:), ...
    'params', struct('winLen', winLen, 'hopLen', hopLen, 'gammaComp', gammaComp, ...
                     'avgHalfWin', avgHalfWin, 'peakThresh', peakThresh, 'minSepSec', minSepSec));

if showPlot
    plotOnsets(audioFile, audioMono, fs, frameTime, logSpec, specFlux, localAvg, ...
               enhFlux, peakTimes, pkHeights, onsetInfo.params);
end
end

function plotOnsets(audioFile, audioMono, fs, frameTime, logSpec, specFlux, localAvg, ...
                    enhFlux, peakTimes, pkHeights, prm)
audioTime = (0:length(audioMono)-1) / fs;
freqAxis  = (0:size(logSpec,1)-1) * fs / prm.winLen;
timeLims  = [frameTime(1) audioTime(end)];

figure('Position', [100 100 900 800]);

subplot(4,1,1);
plot(audioTime, audioMono);
xlim(timeLims); ylabel('Amplitude');
title('Waveform');

subplot(4,1,2);
imagesc(frameTime, freqAxis/1000, logSpec); axis xy;
xlim(timeLims); ylim([0 5]); ylabel('Freq (kHz)');
title(sprintf('Log-Compressed Spectrogram (\\gamma = %d)', prm.gammaComp));

subplot(4,1,3);
plot(frameTime, specFlux); hold on;
plot(frameTime, localAvg, 'r', 'LineWidth', 1.5);
xlim(timeLims); ylabel('\Delta_{spectral}');
legend('Spectral flux', sprintf('Local avg (M = %d)', prm.avgHalfWin));
title('Spectral Flux');

subplot(4,1,4);
plot(frameTime, enhFlux); hold on;
plot(peakTimes, pkHeights, 'rv', 'MarkerFaceColor', 'r');
yline(prm.peakThresh, 'k--');
xlim(timeLims); xlabel('Time (s)'); ylabel('Normalized');
legend('Enhanced flux', 'Detected onsets', 'Threshold');
title('Enhanced Spectral Flux + Detected Onsets');

sgtitle(fileTitle(audioFile), 'Interpreter', 'none', 'FontWeight', 'bold');
end

function figTitle = fileTitle(audioFile)
% e.g. Cmaj_120_qrtr_asc.wav -> "Cmaj_120_qrtr_asc.wav  -  Cmaj | 120 BPM | quarter, ascending"
[~, baseName, fileExt] = fileparts(audioFile);
nameParts = split(baseName, '_');
if numel(nameParts) < 2
    figTitle = [baseName fileExt];
    return
end
termMap = containers.Map({'qrtr','eighth','asc','dsc'}, ...
                         {'quarter','eighth','ascending','descending'});
descTerms = nameParts(3:end);
for t = 1:numel(descTerms)
    if isKey(termMap, descTerms{t})
        descTerms{t} = termMap(descTerms{t});
    end
end
descStr  = strrep(strjoin(descTerms, ', '), 'ascending, descending', 'ascending/descending');
figTitle = sprintf('%s%s  -  %s | %s BPM | %s', baseName, fileExt, ...
    nameParts{1}, nameParts{2}, descStr);
end
