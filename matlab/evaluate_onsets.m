function results = evaluate_onsets(detOnsets, refOnsets, tolSec, verbose)
%EVALUATE_ONSETS Score detected onsets against reference (MIDI) onsets
%   results = evaluate_onsets(detOnsets, refOnsets)
%   results = evaluate_onsets(detOnsets, refOnsets, tolSec, verbose)
%
%   detOnsets, refOnsets - onset times (s)
%   tolSec  - match tolerance in seconds (default 0.05)
%   verbose - print a summary (default true)
%
%   One-to-one greedy matching: closest pairs within tolerance are matched
%   first, and each onset can be used at most once.

if nargin < 3 || isempty(tolSec)
    tolSec = 0.05;
end
if nargin < 4
    verbose = true;
end
detOnsets = detOnsets(:);
refOnsets = refOnsets(:);
numDet = numel(detOnsets);
numRef = numel(refOnsets);

%% Greedy one-to-one matching
errMat = detOnsets - refOnsets.';            % numDet x numRef
candIdx = find(abs(errMat) <= tolSec);
[~, sortOrd] = sort(abs(errMat(candIdx)));
candIdx = candIdx(sortOrd);

detUsed = false(numDet, 1);
refUsed = false(numRef, 1);
matchPairs = zeros(0, 2);                    % [detIdx, refIdx]
for c = 1:numel(candIdx)
    [dIdx, rIdx] = ind2sub(size(errMat), candIdx(c));
    if ~detUsed(dIdx) && ~refUsed(rIdx)
        matchPairs(end+1, :) = [dIdx, rIdx]; %#ok<AGROW>
        detUsed(dIdx) = true;
        refUsed(rIdx) = true;
    end
end
matchPairs = sortrows(matchPairs, 2);

%% Metrics
numTP = size(matchPairs, 1);
numFP = numDet - numTP;
numFN = numRef - numTP;
precVal = numTP / max(numDet, 1);
recVal  = numTP / max(numRef, 1);
if precVal + recVal > 0
    f1Val = 2*precVal*recVal / (precVal + recVal);
else
    f1Val = 0;
end

matchErr = detOnsets(matchPairs(:,1)) - refOnsets(matchPairs(:,2));
if isempty(matchErr)
    [meanErr, meanAbsErr, maxAbsErr] = deal(NaN);
else
    meanErr    = mean(matchErr);
    meanAbsErr = mean(abs(matchErr));
    maxAbsErr  = max(abs(matchErr));
end

results = struct('numDet', numDet, 'numRef', numRef, ...
    'TP', numTP, 'FP', numFP, 'FN', numFN, ...
    'precision', precVal, 'recall', recVal, 'f1', f1Val, ...
    'meanErr', meanErr, 'meanAbsErr', meanAbsErr, 'maxAbsErr', maxAbsErr, ...
    'matchPairs', matchPairs, 'matchErr', matchErr, ...
    'falsePosIdx', find(~detUsed), 'missedRefIdx', find(~refUsed), ...
    'tolSec', tolSec);

if verbose
    fprintf('Onsets: %d detected, %d reference | TP %d  FP %d  FN %d\n', ...
        numDet, numRef, numTP, numFP, numFN);
    fprintf('P = %.3f  R = %.3f  F1 = %.3f | mean err %.1f ms, max |err| %.0f ms\n', ...
        precVal, recVal, f1Val, 1000*meanErr, 1000*maxAbsErr);
end
end
