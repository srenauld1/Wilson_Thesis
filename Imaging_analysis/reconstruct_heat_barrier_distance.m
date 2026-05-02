% ── Load your data ────────────────────────────────────────────────────────
% data struct is already in workspace from makeVirmenCSV
% If loading from file:
% load('your_file.mat');   % loads data struct

% ── Run reconstruction using confirmed parameters ─────────────────────────
travel_integration_time = 5;    % seconds — matches vr.tt(1).travel_integration_time
frame_rate = 60;   % = 60.1663 — use exact value from workspace

% Pull the arrays we need directly from the data struct
currTt_full   = data.currTt;                    % 1x53839
hbd_short     = data.heat_barrier_distance;     % 1x36401 — needs fixing
n_full        = length(currTt_full);            % 53839
n_actual      = length(hbd_short);             % 36401

fprintf('Full length:              %d\n', n_full);
fprintf('heat_barrier_distance:    %d\n', n_actual);
fprintf('Missing frames:           %d\n', n_full - n_actual);

% ── Build valid frame mask from trial structure ───────────────────────────
integration_frames = round(travel_integration_time * frame_rate);  % = 300
fprintf('Integration window:       %d frames\n', integration_frames);

should_be_valid = false(1, n_full);

% Find trial boundaries from currTt
transitions = [1, find(diff(currTt_full) ~= 0) + 1, n_full + 1];
n_segs      = length(transitions) - 1;

fprintf('\nTrial segments found: %d\n', n_segs);
fprintf('%-10s %-10s %-12s %-12s %-12s\n', ...
    'Segment','Type','Start','End','Valid from');

for i = 1:n_segs
    seg_start  = transitions(i);
    seg_end    = transitions(i+1) - 1;
    trial_type = currTt_full(seg_start);

    if trial_type == 1 || trial_type == 2
        valid_start = seg_start + integration_frames;
        if valid_start <= seg_end
            should_be_valid(valid_start:seg_end) = true;
            fprintf('%-10d %-10d %-12d %-12d %-12d\n', ...
                i, trial_type, seg_start, seg_end, valid_start);
        else
            fprintf('%-10d %-10d %-12d %-12d %-12s\n', ...
                i, trial_type, seg_start, seg_end, 'too short');
        end
    else
        fprintf('%-10d %-10d %-12d %-12d %-12s\n', ...
            i, trial_type, seg_start, seg_end, 'NaN (accom)');
    end
end

% ── Verify expected vs actual count ──────────────────────────────────────
n_expected = sum(should_be_valid);
fprintf('\nExpected valid frames:    %d\n', n_expected);
fprintf('Actual  valid frames:     %d\n', n_actual);

if n_expected == n_actual
    fprintf('✓ Counts match — reconstruction is exact\n');
else
    fprintf('✗ Count mismatch of %d frames\n', abs(n_expected - n_actual));
    fprintf('  This is likely from trials that ended early via barrier escape\n');
    fprintf('  Check control_trial and currTt to identify short trials\n');
end

% ── Reconstruct full-length array ────────────────────────────────────────
hbd_full = ones(1, n_full) * 100000;   % matches vr.distance_to_barrier init value

if n_expected == n_actual
    % Exact match — place values precisely
    hbd_full(should_be_valid) = hbd_short;
else
    % Mismatch — sequential fill with warning
    warning('Using sequential fill — pre-barrier NaN windows will be approximate');
    valid_idx = find(should_be_valid);
    n_fill    = min(n_actual, length(valid_idx));
    hbd_full(valid_idx(1:n_fill)) = hbd_short(1:n_fill);
end

% ── Write back to data struct ─────────────────────────────────────────────
data.heat_barrier_distance = hbd_full;

% ── Final length check ────────────────────────────────────────────────────
fields_to_check = {'barJump','heat','currTt','control_trial','heat_barrier_distance'};
fprintf('\nFinal length check:\n');
for f = fields_to_check
    L = length(data.(f{1}));
    match = (L == n_full);
    fprintf('  %-30s %d  %s\n', f{1}, L, ...
        char('✓' * match + '✗' * ~match));
end

% ── Sanity check distance values ─────────────────────────────────────────
valid_vals = hbd_full(~isnan(hbd_full));
fprintf('\nDistance value sanity check:\n');
fprintf('  Min:    %.2f\n', min(valid_vals));
fprintf('  Max:    %.2f\n', max(valid_vals));
fprintf('  Mean:   %.2f\n', mean(valid_vals));
fprintf('  Median: %.2f\n', median(valid_vals));
fprintf('  %% NaN:  %.1f%%\n', 100*sum(isnan(hbd_full))/n_full);