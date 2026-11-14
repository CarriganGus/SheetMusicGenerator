function [refNotes, midiInfo] = read_midi_notes(midiFile)
%READ_MIDI_NOTES Parse a Standard MIDI File into a note matrix
%   [refNotes, midiInfo] = read_midi_notes(midiFile)
%
%   refNotes - K x 3 [onset (s), offset (s), MIDI pitch], sorted by onset
%   midiInfo - struct: format, ticksPerQuarter, bpm, timeSig, tempoMap
%   Channel 10 (drums) is ignored.

if ~isfile(midiFile)
    error('read_midi_notes:fileNotFound', 'File not found: %s', midiFile);
end
fid = fopen(midiFile, 'r');
midiBytes = fread(fid, inf, 'uint8=>double')';
fclose(fid);

%% Header chunk
if ~strcmp(char(midiBytes(1:4)), 'MThd')
    error('read_midi_notes:badHeader', 'Not a MIDI file: %s', midiFile);
end
hdrLen     = bigEndian(midiBytes(5:8));
midiFormat = bigEndian(midiBytes(9:10));
numTracks  = bigEndian(midiBytes(11:12));
ticksPerQ  = bigEndian(midiBytes(13:14));
if ticksPerQ >= 32768
    error('read_midi_notes:smpte', 'SMPTE time division not supported');
end
chunkPos = 9 + hdrLen;

noteEvents  = zeros(0, 4);   % [tick, isNoteOn, pitch, channel]
tempoEvents = zeros(0, 2);   % [tick, microseconds per quarter]
timeSig     = [];

%% Track chunks
for trk = 1:numTracks
    if ~strcmp(char(midiBytes(chunkPos:chunkPos+3)), 'MTrk')
        error('read_midi_notes:badTrack', 'Missing MTrk chunk %d', trk);
    end
    trkLen = bigEndian(midiBytes(chunkPos+4:chunkPos+7));
    bytePos = chunkPos + 8;
    trkEnd  = bytePos + trkLen;
    absTick = 0;
    runStatus = 0;

    while bytePos < trkEnd
        [deltaTick, bytePos] = readVLQ(midiBytes, bytePos);
        absTick = absTick + deltaTick;
        statByte = midiBytes(bytePos);

        if statByte == 255                                % meta event
            metaType = midiBytes(bytePos+1);
            [metaLen, bytePos] = readVLQ(midiBytes, bytePos+2);
            metaData = midiBytes(bytePos:bytePos+metaLen-1);
            bytePos  = bytePos + metaLen;
            if metaType == 81                             % set tempo
                tempoEvents(end+1, :) = [absTick, bigEndian(metaData)]; %#ok<AGROW>
            elseif metaType == 88 && isempty(timeSig)     % time signature
                timeSig = [metaData(1), 2^metaData(2)];
            elseif metaType == 47                         % end of track
                break
            end

        elseif statByte == 240 || statByte == 247         % sysex
            [sysLen, bytePos] = readVLQ(midiBytes, bytePos+1);
            bytePos = bytePos + sysLen;

        else                                              % channel event
            if statByte >= 128
                runStatus = statByte;
                bytePos = bytePos + 1;
            elseif runStatus == 0
                error('read_midi_notes:badEvent', 'Data byte with no running status');
            end
            evType  = bitshift(runStatus, -4);
            channel = bitand(runStatus, 15) + 1;
            if evType == 8 || evType == 9
                pitch    = midiBytes(bytePos);
                velocity = midiBytes(bytePos+1);
                bytePos  = bytePos + 2;
                if channel ~= 10
                    noteEvents(end+1, :) = [absTick, evType == 9 && velocity > 0, pitch, channel]; %#ok<AGROW>
                end
            elseif evType == 12 || evType == 13           % program change, channel pressure
                bytePos = bytePos + 1;
            else
                bytePos = bytePos + 2;
            end
        end
    end
    chunkPos = trkEnd;
end

%% Tempo map (default 120 BPM until the first tempo event)
tempoEvents = sortrows(tempoEvents, 1);
if isempty(tempoEvents) || tempoEvents(1,1) > 0
    tempoEvents = [0, 500000; tempoEvents];
end
if isempty(timeSig)
    timeSig = [4 4];
end

%% Pair note-ons with note-offs (offs first at equal ticks for repeated notes)
noteEvents = sortrows(noteEvents, [1 2]);
openOnTick = nan(16, 128);
noteTicks  = zeros(0, 3);    % [onTick, offTick, pitch]
for ev = 1:size(noteEvents, 1)
    evTick = noteEvents(ev,1);
    chIdx  = noteEvents(ev,4);
    pIdx   = noteEvents(ev,3) + 1;
    if ~isnan(openOnTick(chIdx, pIdx))
        noteTicks(end+1, :) = [openOnTick(chIdx, pIdx), evTick, pIdx-1]; %#ok<AGROW>
        openOnTick(chIdx, pIdx) = NaN;
    end
    if noteEvents(ev,2)
        openOnTick(chIdx, pIdx) = evTick;
    end
end
[openCh, openP] = find(~isnan(openOnTick));
if ~isempty(openCh)
    lastTick = max(noteEvents(:,1));
    for k = 1:numel(openCh)
        noteTicks(end+1, :) = [openOnTick(openCh(k), openP(k)), lastTick, openP(k)-1]; %#ok<AGROW>
    end
end

%% Ticks -> seconds
refNotes = [ticksToSec(noteTicks(:,1), tempoEvents, ticksPerQ), ...
            ticksToSec(noteTicks(:,2), tempoEvents, ticksPerQ), ...
            noteTicks(:,3)];
refNotes = sortrows(refNotes, [1 3]);

midiInfo = struct('format', midiFormat, 'ticksPerQuarter', ticksPerQ, ...
    'bpm', 60e6 / tempoEvents(1,2), 'timeSig', timeSig, ...
    'tempoMap', [tempoEvents(:,1), ticksToSec(tempoEvents(:,1), tempoEvents, ticksPerQ), ...
                 60e6 ./ tempoEvents(:,2)]);   % [tick, sec, BPM]
end

function secVals = ticksToSec(tickVals, tempoEvents, ticksPerQ)
tempoTicks = tempoEvents(:,1);
usPerQ     = tempoEvents(:,2);
segStart   = [0; cumsum(diff(tempoTicks) .* usPerQ(1:end-1))] / (ticksPerQ * 1e6);
secVals    = zeros(size(tickVals));
for k = 1:numel(tickVals)
    seg = find(tempoTicks <= tickVals(k), 1, 'last');
    secVals(k) = segStart(seg) + (tickVals(k) - tempoTicks(seg)) * usPerQ(seg) / (ticksPerQ * 1e6);
end
end

function [vlqVal, bytePos] = readVLQ(midiBytes, bytePos)
vlqVal = 0;
while true
    curByte = midiBytes(bytePos);
    bytePos = bytePos + 1;
    vlqVal  = vlqVal * 128 + bitand(curByte, 127);
    if curByte < 128
        break
    end
end
end

function val = bigEndian(byteVec)
val = sum(byteVec .* 256.^(numel(byteVec)-1:-1:0));
end
