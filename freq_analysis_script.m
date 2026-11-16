%% run_on_wav.m
% Reads a WAV file, runs freq_analysis, plots the results, and prints a rough
% note list. Put this file, freq_analysis.m, and your WAV in the same folder.

clear; clc; close all;

wavFile = 'Cmaj_60_qrtr_asc.wav';     % <-- change to your file
minNoteDur = 0.10;                 % ignore pitch runs shorter than this (s)

%% Read + analyze (freq_analysis reads the file itself)
out = freq_analysis(wavFile);      % pass a struct as 2nd arg to override defaults
% e.g. out = freq_analysis(wavFile, struct('Nh', 4, 'zpFactor', 1));

[x, fs] = audioread(wavFile);      % only needed here for the waveform plot
x = mean(x, 2);
fprintf('Read %s: %.2f s at %d Hz\n', wavFile, length(x)/fs, fs);

%% Plots
figure('Name', 'Frequency analysis');

subplot(4,1,1);
plot((0:length(x)-1)/fs, x); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Waveform');

subplot(4,1,2);
imagesc(out.tFrames, out.midiAxis, out.Yc.'); axis xy; colorbar;
xlabel('Time (s)'); ylabel('MIDI'); title('Semitone-pooled log spectrogram');

subplot(4,1,3);
imagesc(out.tFrames, out.midiAxis, out.Yh.'); axis xy; colorbar;
xlabel('Time (s)'); ylabel('MIDI'); title('After harmonic summation');

subplot(4,1,4);
plot(out.tFrames, out.midiFrame, '.-'); grid on;
xlabel('Time (s)'); ylabel('MIDI pitch'); title('Detected pitch per frame');

%% Rough note list (runs of constant pitch; replace with onset-based
%  segmentation once the onset detector is merged)
mm = out.midiFrame;
mm(isnan(mm)) = 0;                                  % 0 = unvoiced
starts = find([true; diff(mm) ~= 0]);
ends   = [starts(2:end) - 1; numel(mm)];
names  = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};

fprintf('\n  Start(s)  End(s)   MIDI  Note\n');
for i = 1:numel(starts)
    pitch = mm(starts(i));
    if pitch == 0, continue; end
    t0 = out.tFrames(starts(i));
    t1 = out.tFrames(ends(i));
    if (t1 - t0) < minNoteDur, continue; end
    fprintf('  %7.3f  %7.3f   %3d   %s%d\n', t0, t1, pitch, ...
        names{mod(pitch,12)+1}, floor(pitch/12)-1);
end


%% ------------------------- Local functions (must stay at the end of the script) -------------------------

function out = freq_analysis(x, fs, opts)
%FREQ_ANALYSIS  Long-window STFT pitch analysis (frequency section only).
%
%   out = freq_analysis('song.wav')          % reads the WAV file for you
%   out = freq_analysis('song.wav', opts)
%   out = freq_analysis(x, fs)               % or pass samples directly
%   out = freq_analysis(x, fs, opts)
%
%   Steps (per the "Frequency" and "Segment Notes" sections of the write-up):
%     1. optional highpass filter
%     2. zero-pad N/2 so frame m is centered at t_m = m*H/fs
%     3. long-window STFT
%     4. semitone pooling   Y_LF(n,b) = sum_{k in P(b)} |X(n,k)|^2
%     5. log compression + harmonic summation
%     6. per-frame MIDI pitch = argmax_b Y~(n,b)
%
%   Inputs
%     x     audio vector (mono)
%     fs    sampling rate (Hz)
%     opts  (optional) struct; any field can be omitted:
%           N          window length                    [4096]
%           H          hop in samples                   [round(fs/100) = 441 @ 44.1k]
%                      ** must equal the hop used by the onset STFT **
%           zpFactor   Nfft = zpFactor*N                [4]  (1 = exactly as in report)
%           Nh         harmonics in harmonic summation  [5]
%           gammaLog   log-compression factor (>1)      [1000]
%           hpCutoff   highpass cutoff Hz ([] = off)    [60]
%           minMidi    lowest pitch searched            [36]  (C2)
%           maxMidi    highest pitch searched           [96]  (C7)
%           gateRel    voiced if energy > gateRel*max   [0.02]
%           keepSTFT   return the full complex STFT     [true]
%
%   Output struct "out"
%     tFrames    numFrames x 1   frame times (s), t_m = m*H/fs   <- shared time grid
%     midiAxis   1 x nb          MIDI number for each semitone column
%     Y_LF       numFrames x nb  semitone-pooled power spectrogram
%     Yc         numFrames x nb  log-compressed Y_LF
%     Yh         numFrames x nb  harmonic-summed salience Y~
%     midiFrame  numFrames x 1   detected MIDI pitch (NaN where unvoiced)
%     f0Hz       numFrames x 1   detected f0 in Hz (NaN where unvoiced)
%     voiced     numFrames x 1   logical gate
%     X, fk      full STFT (numFrames x nBins) and bin frequencies (Hz), if keepSTFT
%     params     the parameters actually used (for logging / merging)

    %% 0. Input handling: filename or raw samples
    if ischar(x) || isstring(x)
        if nargin >= 2 && isstruct(fs)
            opts = fs;                       % freq_analysis('file.wav', opts)
        else
            opts = struct();
        end
        [x, fs] = audioread(char(x));        % read the WAV file
    elseif nargin < 3
        opts = struct();
    end
    if size(x, 2) > 1
        x = mean(x, 2);                      % stereo -> mono
    end
    x = x(:);
    p = setDefaults(opts, fs);

    %% 1. Highpass filter (zero-phase so timing is not shifted)
    if ~isempty(p.hpCutoff)
        try
            [bh, ah] = butter(4, p.hpCutoff/(fs/2), 'high');
            x = filtfilt(bh, ah, x);
        catch
            warning('freq_analysis:noToolbox', ...
                'Signal Processing Toolbox not found; skipping highpass filter.');
        end
    end

    %% 2-3. Zero-padded long-window STFT
    N    = p.N;
    H    = p.H;
    Nfft = p.zpFactor * N;
    nBins = Nfft/2 + 1;

    xp        = [zeros(N/2, 1); x; zeros(N, 1)];
    numFrames = floor(length(x)/H) + 1;
    w         = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);      % periodic Hann (no toolbox needed)

    X = zeros(numFrames, nBins);
    for m = 0:numFrames-1
        seg = xp(m*H + (1:N)) .* w;
        Xm  = fft(seg, Nfft);
        X(m+1, :) = Xm(1:nBins).';
    end
    tFrames = (0:numFrames-1)' * H / fs;
    fk      = (0:nBins-1) * fs / Nfft;

    %% 4. Semitone pooling:  P(b) = {k : Bin(f_k) = b},  Bin(f) = floor(12*log2(f/440)+69.5)
    fLo   = 440 * 2.^((p.minMidi - 0.5 - 69)/12);
    fHi   = 440 * 2.^((p.maxMidi + 0.5 - 69)/12);
    kIdx  = find(fk >= fLo & fk < fHi);
    bMidi = floor(12*log2(fk(kIdx)/440) + 69.5);
    nb    = p.maxMidi - p.minMidi + 1;
    bIdx  = min(max(bMidi - p.minMidi + 1, 1), nb);  % clamp for floating-point edge cases

    Pmat = sparse(bIdx, 1:numel(kIdx), 1, nb, numel(kIdx));   % nb x K
    Y_LF = abs(X(:, kIdx)).^2 * Pmat.';                       % numFrames x nb

    %% 5. Log compression + harmonic summation
    Yn = Y_LF / max(max(Y_LF(:)), eps);
    Yc = log(1 + p.gammaLog * Yn);

    shifts = round(12*log2(1:p.Nh));                  % [0 12 19 24 28 ...]
    Ypad   = [Yc, zeros(numFrames, max(shifts))];
    Yh     = zeros(size(Yc));
    for h = 1:p.Nh
        Yh = Yh + Ypad(:, (1:nb) + shifts(h));
    end

    %% 6. Pitch per frame
    [~, iMax] = max(Yh, [], 2);
    midiFrame = p.minMidi + iMax - 1;

    frameEnergy = sum(Y_LF, 2);
    voiced = frameEnergy > p.gateRel * max(frameEnergy);
    midiFrame(~voiced) = NaN;

    f0Hz = 440 * 2.^((midiFrame - 69)/12);            % elementwise power (vector input)

    %% Package output
    out.tFrames   = tFrames;
    out.midiAxis  = p.minMidi:p.maxMidi;
    out.Y_LF      = Y_LF;
    out.Yc        = Yc;
    out.Yh        = Yh;
    out.midiFrame = midiFrame;
    out.f0Hz      = f0Hz;
    out.voiced    = voiced;
    if p.keepSTFT
        out.X  = X;
        out.fk = fk;
    end
    out.params = p;
end

function p = setDefaults(opts, fs)
    d.N        = 4096;
    d.H        = round(fs/100);     % 441 at 44.1 kHz; override to match the onset STFT
    d.zpFactor = 4;
    d.Nh       = 5;
    d.gammaLog = 1000;
    d.hpCutoff = 60;
    d.minMidi  = 36;
    d.maxMidi  = 96;
    d.gateRel  = 0.02;
    d.keepSTFT = true;

    p = d;
    f = fieldnames(opts);
    for i = 1:numel(f)
        p.(f{i}) = opts.(f{i});
    end
end
