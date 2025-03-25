% processFixationPeriods
% This function processes fixation periods in an experiment by merging short non-fixation gaps
% between consecutive fixation periods, and removing fixation periods that are too short.
% It assumes the input `fixIdx` is a logical matrix indicating fixation (true) and non-fixation (false) periods.
%
% INPUTS:
%   fixIdx       - 3D logical array (time x trial x condition) indicating fixation periods
%   expttime     - Time vector (seconds)
%   max_gap_dur  - Maximum allowable gap duration (seconds) to merge short gaps between fixations
%   min_fix_dur  - Minimum allowable fixation duration (seconds) to keep a fixation period
%
% OUTPUTS:
%   fixIdx       - 3D logical array (same size as input) with processed fixation periods
%
% CREATED: [Date] MC
%
function fixIdx = processFixationPeriods(fixIdx, expttime, max_gap_dur, min_fix_dur)
    % Convert max gap duration and minimum fixation duration to indices
    max_gap_idx = fetchTimeIdx(expttime, max_gap_dur); % Convert the max gap duration to indices
    min_fix_idx = fetchTimeIdx(expttime, min_fix_dur); % Convert the minimum fixation duration to indices
    
    % Convert 0.5 seconds to indices for removing short epochs
    min_epoch_dur = 0.5; % Minimum epoch duration in seconds
    min_epoch_idx = fetchTimeIdx(expttime, min_epoch_dur); % Convert 0.5s duration to indices

    % Loop through trials and conditions (assuming 3D fixIdx: time x trial x condition)
    for trial = 1:size(fixIdx, 2)
        for cond = 1:size(fixIdx, 3)
            % Find start and end of fixation periods
            fix_diff = diff([0; fixIdx(:, trial, cond); 0]); % Edge detection
            fix_starts = find(fix_diff == 1); % Start of fixation
            fix_ends = find(fix_diff == -1) - 1; % End of fixation

            % Remove fixation periods shorter than 0.5 seconds (min_epoch_idx)
            valid_epoch = (fix_ends - fix_starts + 1) >= min_epoch_idx;
            fix_starts = fix_starts(valid_epoch); % Keep only valid start indices
            fix_ends = fix_ends(valid_epoch); % Keep only valid end indices

            % Merge short non-fixation gaps between fixation periods
            i = 1; % Initialize index for loop
            while i < length(fix_starts) % Loop through fixation periods
                gap_duration = fix_starts(i + 1) - fix_ends(i) - 1;
                % If the gap between two fixations is shorter than max_gap_idx, merge them
                if gap_duration <= max_gap_idx
                    fix_ends(i) = fix_ends(i + 1); % Merge the fixation ends
                    fix_starts(i + 1) = []; % Remove the next start
                    fix_ends(i + 1) = []; % Remove the next end
                else
                    i = i + 1; % Only increment if no merge is done
                end
            end

            % Remove fixation periods shorter than the minimum fixation duration (min_fix_dur)
            valid_fix = (fix_ends - fix_starts + 1) >= min_fix_idx;
            fix_starts = fix_starts(valid_fix); % Keep only valid start indices
            fix_ends = fix_ends(valid_fix); % Keep only valid end indices

            % Reset fixation indices for this trial and condition
            fixIdx(:, trial, cond) = false; % Clear current fixation indices
            for i = 1:length(fix_starts)
                fixIdx(fix_starts(i):fix_ends(i), trial, cond) = true; % Set the valid fixation periods
            end
        end
    end
end
