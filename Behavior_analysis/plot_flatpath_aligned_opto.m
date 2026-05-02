function exptData = plot_flatpath_aligned_opto(exptData, savepath)
% plot_flatpath_aligned_opto
% Align flat paths around optostim bouts, translating & rotating so that
% the point immediately before opto onset is at (0,0) and the movement
% direction at onset is along +X.
%
% For each opto pulse:
%   - 2 s pre window (before onset): gray
%   - opto period (optoStim > 0): light red
%   - 2 s post window (after offset): blue
% Also plots average aligned trajectory across pulses:
%   - pre segment: black
%   - opto segment: bright red
%   - post segment: blue
%
% Additionally:
%   - Computes net heading displacement from opto start to opto end
%     (unwrapped) for each used pulse, stored in:
%         exptData.opto.headingdisp   [nPulses x 1] (deg)
%   - Prints the mean net displacement to the command window and
%     annotates it on the figure.
%
% Assumes (after process_fictrac_panels):
%   exptData.x, exptData.y, exptData.t, exptData.optoStim are all the
%   same length (downsampled behavioral rate). Handles the general case
%   where t/optoStim differ from x/y length.

%% ------------- PARAMETERS -------------
pre_win_sec   = 2;      % sec before onset
post_win_sec  = 2;      % sec after offset
min_step_dist = 1e-6;   % min displacement to trust direction
min_opto_pts  = 3;      % min samples during opto
%% --------------------------------------

% ---------- Validate required fields ----------
required = {'x', 'y', 't', 'optoStim'};
for f = required
    if ~isfield(exptData, f{1})
        error('plot_flatpath_aligned_opto: missing field exptData.%s', f{1});
    end
end

% Optional but needed for headingdisp
have_heading = isfield(exptData, 'headingPosition');
if ~have_heading
    warning('plot_flatpath_aligned_opto: exptData.headingPosition not found; heading displacement will not be computed.');
end

% Behavioral signals (lower-rate)
x_beh = exptData.x(:);
y_beh = exptData.y(:);
Nb    = numel(x_beh);

% High-rate signals
t_hi    = exptData.t_ds(:);
opto_hi = exptData.optoStim_ds(:);
No      = numel(opto_hi);

if numel(t_hi) ~= No
    error('Length of exptData.t (%d) must match exptData.optoStim (%d).', ...
          numel(t_hi), No);
end

% ---------- Align opto to behavior timebase ----------
if numel(opto_hi) == Nb
    % Already aligned — use directly
    t    = t_hi;
    opto = opto_hi;
else
    % Different lengths — resample opto onto behavior timebase
    warning(['plot_flatpath_aligned_opto: t/optoStim length (%d) does not ' ...
             'match x/y length (%d). Resampling.'], No, Nb);
    t    = linspace(t_hi(1), t_hi(end), Nb).';
    opto = interp1(t_hi, opto_hi, t, 'nearest', 0);
end

% Unified behavior variables
x    = x_beh;
y    = y_beh;

% ---------- Heading (for net displacement) ----------
if have_heading
    heading_raw = exptData.headingPosition(:);  % assume degrees
    if numel(heading_raw) ~= Nb
        warning('plot_flatpath_aligned_opto: exptData.headingPosition length (%d) ~= x/y length (%d); skipping headingdisp.', ...
                numel(heading_raw), Nb);
        have_heading = false;
    else
        % unwrap via radians to handle wrap-around
        h_rad       = deg2rad(heading_raw);
        h_unwrap_r  = unwrap(h_rad);
        heading_unw = rad2deg(h_unwrap_r);   % unwrapped degrees
    end
end

% ---------- Timing info ----------
dt   = median(diff(t));
rate = 1/dt;

pre_pts  = round(pre_win_sec  * rate);
post_pts = round(post_win_sec * rate);

% ---------- Find opto ON segments ----------
is_on      = opto > 0;
change_pts = [1; find(diff(is_on) ~= 0) + 1; Nb];

onsets  = [];
offsets = [];
for k = 1:numel(change_pts)-1
    seg_start = change_pts(k);
    seg_end   = change_pts(k+1) - 1;
    if is_on(seg_start)
        onsets(end+1)  = seg_start; %#ok<AGROW>
        offsets(end+1) = seg_end;   %#ok<AGROW>
    end
end

if isempty(onsets)
    warning('plot_flatpath_aligned_opto: no optostim pulses found after downsampling.');
    return;
end

% ---------- Colors ----------
col_pre      = [0.5 0.5 0.5];   % gray       — individual pre trials
col_post     = [0.5 0.5 1.0];   % blue       — individual post trials
col_opto     = [1.0 0.5 0.5];   % light red  — individual opto trials

% ---------- Set up figure ----------
fig = figure('Name', 'Aligned opto flat paths (2s pre/post)', 'Color', 'w');
hold on;

plotted_any  = false;

% Handles for legend — captured on first instance only
h_pre        = [];
h_post       = [];
h_opto_trial = [];

% Storage for averaging
pre_X_cells  = {};
pre_Y_cells  = {};
opto_X_cells = {};
opto_Y_cells = {};
post_X_cells = {};
post_Y_cells = {};

% Storage for heading displacement (only for pulses actually used)
heading_disp = [];  % will grow as we accept pulses

% ---------- Process each opto pulse ----------
for tr = 1:numel(onsets)
    seg_start = onsets(tr);
    seg_end   = offsets(tr);

    if seg_start > Nb || seg_end > Nb
        continue;
    end

    pre_start  = seg_start - pre_pts;
    pre_end    = seg_start - 1;
    post_start = seg_end   + 1;
    post_end   = seg_end   + post_pts;

    % Skip pulses without a full pre/post window
    if pre_start < 1 || post_end > Nb
        continue;
    end

    idx_pre  = pre_start:pre_end;
    idx_opto = seg_start:seg_end;
    idx_post = post_start:post_end;

    if numel(idx_opto) < min_opto_pts
        continue;
    end

    % Reference: last pre sample → becomes (0,0)
    ref_idx = pre_end;
    x0 = x(ref_idx);
    y0 = y(ref_idx);

    % Direction at onset → aligned to +X
    dx = x(seg_start) - x(ref_idx);
    dy = y(seg_start) - y(ref_idx);

    if hypot(dx, dy) < min_step_dist
        theta = 0;
    else
        theta = -atan2(dy, dx);
    end

    % Translate and rotate the full window
    idx_all = [idx_pre, idx_opto, idx_post];
    X = x(idx_all) - x0;
    Y = y(idx_all) - y0;

    R      = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    XY_rot = R * [X.'; Y.'];
    Xr     = XY_rot(1,:);
    Yr     = XY_rot(2,:);

    % Split back
    n_pre  = numel(idx_pre);
    n_opto = numel(idx_opto);

    Xr_pre  = Xr(1:n_pre);
    Yr_pre  = Yr(1:n_pre);
    Xr_opto = Xr(n_pre+1 : n_pre+n_opto);
    Yr_opto = Yr(n_pre+1 : n_pre+n_opto);
    Xr_post = Xr(n_pre+n_opto+1 : end);
    Yr_post = Yr(n_pre+n_opto+1 : end);

    % Plot individual trial segments; capture first-instance handles
    hp = plot(Xr_pre,  Yr_pre,  '-', 'Color', col_pre,  'LineWidth', 1);
    hq = plot(Xr_post, Yr_post, '-', 'Color', col_post, 'LineWidth', 1);
    ho = plot(Xr_opto, Yr_opto, '-', 'Color', col_opto, 'LineWidth', 1.5);

    if isempty(h_pre),        h_pre        = hp; end
    if isempty(h_post),       h_post       = hq; end
    if isempty(h_opto_trial), h_opto_trial = ho; end

    % Store for averaging
    pre_X_cells{end+1}  = Xr_pre;  %#ok<AGROW>
    pre_Y_cells{end+1}  = Yr_pre;  %#ok<AGROW>
    opto_X_cells{end+1} = Xr_opto; %#ok<AGROW>
    opto_Y_cells{end+1} = Yr_opto; %#ok<AGROW>
    post_X_cells{end+1} = Xr_post; %#ok<AGROW>
    post_Y_cells{end+1} = Yr_post; %#ok<AGROW>

    % --- Net heading displacement (absolute) for this used pulse ---
    if have_heading
        hd_start = heading_unw(seg_start);
        hd_end   = heading_unw(seg_end);
        this_disp = abs(hd_end - hd_start);   % magnitude only
        heading_disp(end+1,1) = this_disp;    %#ok<AGROW>
    end

    plotted_any = true;
end

% ---------- Guard: nothing valid was plotted ----------
if ~plotted_any
    warning('plot_flatpath_aligned_opto: no valid pulses (insufficient pre/post data).');
    close(fig);
    return;
end

% ---------- Attach heading displacement to exptData ----------
if have_heading && ~isempty(heading_disp)
    exptData.opto_sliced_headingdisp = heading_disp;   % [nUsedPulses x 1], in degrees
end

% ---------- Average path ----------
min_pre_L  = min(cellfun(@numel, pre_X_cells));
min_opto_L = min(cellfun(@numel, opto_X_cells));
min_post_L = min(cellfun(@numel, post_X_cells));

nTrials = numel(pre_X_cells);

pre_X_mat  = nan(nTrials, min_pre_L);
pre_Y_mat  = nan(nTrials, min_pre_L);
opto_X_mat = nan(nTrials, min_opto_L);
opto_Y_mat = nan(nTrials, min_opto_L);
post_X_mat = nan(nTrials, min_post_L);
post_Y_mat = nan(nTrials, min_post_L);

for k = 1:nTrials
    % Align pre to its end (reference point)
    pre_X_mat(k,:)  = pre_X_cells{k}(end-min_pre_L+1:end);
    pre_Y_mat(k,:)  = pre_Y_cells{k}(end-min_pre_L+1:end);
    opto_X_mat(k,:) = opto_X_cells{k}(1:min_opto_L);
    opto_Y_mat(k,:) = opto_Y_cells{k}(1:min_opto_L);
    post_X_mat(k,:) = post_X_cells{k}(1:min_post_L);
    post_Y_mat(k,:) = post_Y_cells{k}(1:min_post_L);
end

% this is where the mean is calculated! this is where you may want to
% change things to do displacement or average
mean_pre_X  = mean(pre_X_mat,  1, 'omitnan');
mean_pre_Y  = mean(pre_Y_mat,  1, 'omitnan');
mean_opto_X = mean(opto_X_mat, 1, 'omitnan');
mean_opto_Y = mean(opto_Y_mat, 1, 'omitnan');
mean_post_X = mean(post_X_mat, 1, 'omitnan');
mean_post_Y = mean(post_Y_mat, 1, 'omitnan');



% Plot means: pre = black, opto = bright red, post = blue
h_mean_pre  = plot(mean_pre_X,  mean_pre_Y,  'k-',             'LineWidth', 3);
h_mean_opto = plot(mean_opto_X, mean_opto_Y, 'Color', [1 0 0], 'LineWidth', 3);
h_mean_post = plot(mean_post_X, mean_post_Y, 'Color', [0 0 1], 'LineWidth', 3);

% Mark origin (= reference point at opto onset)
h_origin = plot(0, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);

% ---------- Legend using explicit handles ----------
legend([h_pre, h_opto_trial, h_post, ...
        h_mean_pre, h_mean_opto, h_mean_post, h_origin], ...
       {'Pre 2s (trials)', 'Opto-on (trials)', 'Post 2s (trials)', ...
        'Mean pre', 'Mean opto', 'Mean post', 'Origin'}, ...
       'Location', 'bestoutside');

axis equal;
xlabel('Aligned X (arb.)');
ylabel('Aligned Y (arb.)');
title(sprintf('Flat path aligned to optostim — 2s pre/post  (n = %d pulses)', nTrials));
grid on;

% ---------- Report mean net heading displacement ----------
if have_heading && ~isempty(heading_disp)
    mean_disp = mean(heading_disp, 'omitnan');
    fprintf('Mean net heading displacement during opto (this fly): %.2f deg\n', mean_disp);

    % annotate on the figure (upper-left corner of current axes)
    ax = gca;
    xl = xlim(ax);
    yl = ylim(ax);
    text(xl(1), yl(2), ...
        sprintf('Mean net heading disp: %.2f deg', mean_disp), ...
        'Parent', ax, ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.7]);  % white-ish box
end

hold off;

% ---------- Save ----------
if ~isempty(savepath)
    save_plot_with_title_as_filename('flatpath_aligned_2s', 'opto', savepath);
end
end