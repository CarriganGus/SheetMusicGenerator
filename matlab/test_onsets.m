%% Batch evaluation of onset detection on every WAV/MIDI pair in data/
clear; clc;

scriptDir  = fileparts(mfilename('fullpath'));
dataDir    = fullfile(scriptDir, '..', 'data');
resultsDir = fullfile(scriptDir, '..', 'results');
tolSec     = 0.05;    % match tolerance (s)

wavList = dir(fullfile(dataDir, 'audio', '**', '*.wav'));
if isempty(wavList)
    error('No WAV files found under %s', fullfile(dataDir, 'audio'));
end

rowList = [];
allErr  = [];
for f = 1:numel(wavList)
    wavPath = fullfile(wavList(f).folder, wavList(f).name);
    [~, baseName] = fileparts(wavList(f).name);
    [~, category] = fileparts(wavList(f).folder);
    if strcmp(category, 'audio')
        midiPath = fullfile(dataDir, 'midi', [baseName '.mid']);
    else
        midiPath = fullfile(dataDir, 'midi', category, [baseName '.mid']);
    end
    if ~isfile(midiPath)
        warning('No matching MIDI for %s - skipped', wavList(f).name);
        continue
    end

    onsetTimes = onset_detection(wavPath, false);
    [refNotes, midiInfo] = read_midi_notes(midiPath);
    res = evaluate_onsets(onsetTimes, refNotes(:,1), tolSec, false);

    newRow = struct('Category', string(category), 'File', string(baseName), ...
        'BPM', round(midiInfo.bpm), 'NumRef', res.numRef, 'NumDet', res.numDet, ...
        'TP', res.TP, 'FP', res.FP, 'FN', res.FN, ...
        'Precision', res.precision, 'Recall', res.recall, 'F1', res.f1, ...
        'MeanErr_ms', 1000*res.meanErr, 'MeanAbsErr_ms', 1000*res.meanAbsErr, ...
        'MaxAbsErr_ms', 1000*res.maxAbsErr);
    rowList = [rowList; newRow]; %#ok<AGROW>
    allErr  = [allErr; res.matchErr]; %#ok<AGROW>
end

resultsTbl = struct2table(rowList);

%% Overall row (pooled over all files)
totTP = sum(resultsTbl.TP);
totFP = sum(resultsTbl.FP);
totFN = sum(resultsTbl.FN);
totP  = totTP / max(totTP + totFP, 1);
totR  = totTP / max(totTP + totFN, 1);
totF1 = 2*totP*totR / max(totP + totR, eps);
overallRow = table("ALL", "Overall", NaN, sum(resultsTbl.NumRef), sum(resultsTbl.NumDet), ...
    totTP, totFP, totFN, totP, totR, totF1, ...
    1000*mean(allErr), 1000*mean(abs(allErr)), 1000*max(abs(allErr)), ...
    'VariableNames', resultsTbl.Properties.VariableNames);
resultsTbl = [resultsTbl; overallRow];

%% Display and save
resultsTbl{:, {'Precision','Recall','F1'}} = round(resultsTbl{:, {'Precision','Recall','F1'}}, 3);
resultsTbl{:, {'MeanErr_ms','MeanAbsErr_ms','MaxAbsErr_ms'}} = ...
    round(resultsTbl{:, {'MeanErr_ms','MeanAbsErr_ms','MaxAbsErr_ms'}}, 1);

fprintf('Onset detection results (tolerance = +/-%.0f ms)\n\n', 1000*tolSec);
disp(resultsTbl);

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end
csvPath = fullfile(resultsDir, 'onset_results.csv');
writetable(resultsTbl, csvPath);
fprintf('Saved to %s\n', csvPath);
