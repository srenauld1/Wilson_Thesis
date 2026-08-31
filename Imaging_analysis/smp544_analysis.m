%% smp544_analysis.m
% Dedicated analysis pipeline for smp544 imaging data.
% Loads a processed *_s_.mat file (same way as master_imaging_analysis),
% decides whether heat was on based on the file name, runs the standard
% fictrac/panels processing, and then makes exploratory plots:
%   - dF/F (right + left), forward and rotational velocity time series
%   - flat path colored by dF/F (+ heat on/off if this is a heat expt)
%   - cross-correlations between each dF/F and each kinematic variable
%   - scatter plots of each dF/F against each kinematic variable
%
% dF/F convention for smp544:
%   roi.ts(1,:) = RIGHT side
%   roi.ts(2,:) = LEFT side

%close all
%clear all

%% initialize variables
jump  = 0;
odor  = 0;   % smp544 is not an odor experiment
heat  = 0;   % set below based on the file name

%% load in daq and ts (same pattern as master_imaging_analysis)
baseFolder = '/Users/sophiarenauld/stacks/';

[fileName, pathName] = uigetfile('*_s_.mat', 'Select Processed Data', baseFolder);
if isequal(fileName, 0)
    disp('User canceled the file selection.');
    return;
end
fullFilePath = fullfile(pathName, fileName);
a2p_data = load(fullFilePath);
disp(['User selected "', fullFilePath, '".']);

%% decide whether heat is on based on the file name
% Heat experiments have "heat" somewhere in the path/file name.
if contains(lower(fullFilePath), 'heat')
    heat = 1;
    disp('Heat experiment detected (found "heat" in file name).');
else
    heat = 0;
    disp('No heat detected in file name -> treating as no-heat experiment.');
end

%% create savepath and folder directory
[base_path, ~, ~] = fileparts(fullFilePath);

stim_directory = fullfile(base_path, 'smr_analysis');
if ~exist(stim_directory, 'dir')
    mkdir(stim_directory);
end
savepath = stim_directory;
disp(['Save path: ', savepath]);

%% process and pick out jumps
vel_on  = 0.5;  % Forward speed above which fly is moving (enter movement)
vel_off = 0.5;  % Forward speed below which fly is NOT moving (exit movement)

% process_fictrac_panels_2p expects a second DAQ pathway (dq(2), the 60 Hz
% behavior stream). smp544 files only have dq(1), so run the shared function
% only when dq(2) exists; otherwise do the minimal dq(1) processing inline.
if numel(a2p_data.dq) >= 2
    a2p_data = process_fictrac_panels_2p(a2p_data, vel_on, vel_off, jump);
else
    % --- remove empty fields (same cleanup the shared function does) ---
    fieldNames = fieldnames(a2p_data.dq);
    for i = 1:numel(fieldNames)
        if isempty(a2p_data.dq(1).(fieldNames{i}))
            a2p_data.dq = rmfield(a2p_data.dq, fieldNames{i});
        end
    end

    % --- yaw + side (rad/s) -> deg/s ---
    a2p_data.dq(1).bvy_deg = rad2deg(a2p_data.dq(1).bvy);
    a2p_data.dq(1).bvs_deg = rad2deg(a2p_data.dq(1).bvs);

    % --- total speed + Schmitt-trigger moving/not-moving binary ---
    speed = a2p_data.dq(1).bvf;
    a2p_data.dq(1).totalspeed = abs(a2p_data.dq(1).bvf) + ...
        (4.5*abs(a2p_data.dq(1).bvy)) + (4.5*abs(a2p_data.dq(1).bvs));

    motion = zeros(size(speed));
    motion(1) = speed(1) >= vel_on;
    for i = 2:numel(speed)
        if motion(i-1) == 0 && speed(i) >= vel_on
            motion(i) = 1;                 % start moving
        elseif motion(i-1) == 1 && speed(i) > vel_off
            motion(i) = 1;                 % keep moving
        else
            motion(i) = 0;                 % stopped
        end
    end
    a2p_data.dq(1).motion.moving_or_not = motion;
    a2p_data.dq(1).motion.ftNotMoveInd  = find(motion == 0);
end

%% align imaging and DAQ time bases
% Some acquisitions can contain more imaging volumes than were recorded by
% the DAQ (for example, if ScanImage continued after the FlyG trial ended).
% Analysis is only valid over the interval shared by both streams. Crop the
% in-memory copy used here; this does not alter the processed input file.
n_dq  = numel(a2p_data.dq(1).t);
n_roi = size(a2p_data.roi.ts, 2);
n_common = min(n_dq, n_roi);

if n_common < 2
    error('smp544_analysis:TooFewSharedSamples', ...
        'DAQ and ROI data have fewer than two shared samples (DAQ=%d, ROI=%d).', ...
        n_dq, n_roi);
end

if n_dq ~= n_roi
    warning('smp544_analysis:LengthMismatch', ...
        ['DAQ and ROI lengths differ (DAQ=%d, ROI=%d). Cropping both to ' ...
         'the first %d shared samples for analysis only.'], ...
        n_dq, n_roi, n_common);

    % Crop every top-level dq(1) vector that uses the DAQ time base. Leave
    % scalar metadata, paths, option structs, and unrelated fields intact.
    dq_fields = fieldnames(a2p_data.dq(1));
    for i = 1:numel(dq_fields)
        field_name = dq_fields{i};
        value = a2p_data.dq(1).(field_name);
        if (isnumeric(value) || islogical(value)) && isvector(value) && ...
                numel(value) == n_dq
            if isrow(value)
                a2p_data.dq(1).(field_name) = value(1:n_common);
            else
                a2p_data.dq(1).(field_name) = value(1:n_common, 1);
            end
        end
    end

    % Keep the movement fields internally consistent with the cropped DAQ.
    if isfield(a2p_data.dq(1), 'motion')
        if isfield(a2p_data.dq(1).motion, 'moving_or_not')
            value = a2p_data.dq(1).motion.moving_or_not;
            if numel(value) >= n_common
                a2p_data.dq(1).motion.moving_or_not = value(1:n_common);
            end
        end
        if isfield(a2p_data.dq(1).motion, 'ftNotMoveInd')
            value = a2p_data.dq(1).motion.ftNotMoveInd;
            a2p_data.dq(1).motion.ftNotMoveInd = value(value <= n_common);
        end
    end

    a2p_data.roi.ts = a2p_data.roi.ts(:, 1:n_common);
end

%% define variables
fwd  = a2p_data.dq(1).bvf(:);        % forward velocity (mm/s)
rot  = a2p_data.dq(1).bvy_deg(:);    % rotational velocity (deg/s)
side = a2p_data.dq(1).bvs_deg(:);    % side/lateral velocity (deg/s)
time = a2p_data.dq(1).t(:);          % imaging-volume time base

dff_right = a2p_data.roi.ts(1,:).';  % roi 1 = right side
dff_left  = a2p_data.roi.ts(2,:).';  % roi 2 = left side

% Package kinematic variables and dF/F for looping in later sections
kin_vars   = {fwd, rot, abs(rot), side};
kin_names  = {'Forward velocity (mm/s)', 'Rotational velocity (deg/s)', ...
              'Rotational speed (deg/s)', 'Side velocity (deg/s)'};
kin_tags   = {'fwd', 'rot', 'rotspeed', 'side'};

dff_vars  = {dff_right, dff_left};
dff_names = {'dF/F right', 'dF/F left'};
dff_tags  = {'right', 'left'};

sample_rate = length(time) / (max(time) - min(time));  % samples per second

% Heat state on the imaging grid (only meaningful for heat experiments).
if heat && isfield(a2p_data.dq(1), 'heat_active')
    heat_mask = a2p_data.dq(1).heat_active(:) > 0;   % true where heat is ON
else
    heat_mask = false(size(time));
end

%% ===== EXPLORATORY: one dF/F trace (right ROI), fwd, rot on one figure =====
figure;

% Subplot 1: dF/F from the first/right ROI only
subplot(3, 1, 1);
plot(time, dff_right, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F (right ROI)');
grid on;

% Subplot 2: Forward Velocity
subplot(3, 1, 2);
plot(time, fwd, 'black', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([0 10]);

% Subplot 3: Rotational Velocity
subplot(3, 1, 3);
plot(time, rot, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on;

linkaxes(findall(gcf, 'Type', 'axes'), 'x');
sgtitle('dF/F (right ROI), Forward, and Rotational Velocity');
save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);

%% ===== EXPLORATORY: flat path =====
a2p_data = plot_flat_path_colored(a2p_data, 'dff', 'data_only', savepath, odor);
figure
plot(a2p_data.dq(1).px, a2p_data.dq(1).py)
xlabel('X position (mm)'); ylabel('Y position (mm)');
title('Flat path'); axis equal; grid on;

%% ===== EXPLORATORY: flat path colored by forward velocity =====
% Use the imaging-grid FicTrac path and velocity so all vectors share the
% same time base. The fixed color range matches the forward-velocity panel.
if odor
    flat_x = a2p_data.dq(1).px(:);
    flat_y = a2p_data.dq(1).py(:);
else
    flat_x = a2p_data.dq(1).pxb(:);
    flat_y = a2p_data.dq(1).pyb(:);
end
n_flat = min([numel(flat_x), numel(flat_y), numel(fwd)]);
flat_x = flat_x(1:n_flat);
flat_y = flat_y(1:n_flat);
flat_fwd = fwd(1:n_flat);

% Use the observed range, but prevent extreme outliers from stretching the
% color scale beyond the interpretable range of -5 to 20 mm/s.
flat_fwd_color = min(max(flat_fwd, -5), 20);
fwd_color_min = min(flat_fwd_color, [], 'omitnan');
fwd_color_max = max(flat_fwd_color, [], 'omitnan');
if fwd_color_max <= fwd_color_min
    fwd_color_max = fwd_color_min + eps(fwd_color_min + 1);
end

f_flat_fwd = figure;
scatter(flat_x, flat_y, 10, flat_fwd_color, 'filled', 'MarkerEdgeColor', 'none');
hold on;
plot(flat_x(1), flat_y(1), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
axis equal; grid on;
xlabel('X position (mm)'); ylabel('Y position (mm)');
title('Flat path colored by forward velocity');
colormap(f_flat_fwd, turbo);
c_fwd = colorbar;
ylabel(c_fwd, 'Forward velocity (mm/s)');
clim([fwd_color_min fwd_color_max]);
savefig(f_flat_fwd, fullfile(savepath, 'flat_path_colored_by_fwd.fig'));
exportgraphics(f_flat_fwd, fullfile(savepath, 'flat_path_colored_by_fwd.png'), ...
    'Resolution', 150, 'BackgroundColor', 'white');

%% ===== EXPLORATORY: flat path with heat on/off + dF/F (heat expts only) =====
if heat
    figure(10); clf;

    % Choose which fictive path to plot: 'visual' (dq.px/py) or 'ball' (dq.pxb/pyb)
    usepath = 'visual';
    if strcmp(usepath, 'ball')
        x_vm = a2p_data.dq(1).pxb(:);  y_vm = a2p_data.dq(1).pyb(:);
    else
        x_vm = a2p_data.dq(1).px(:);   y_vm = a2p_data.dq(1).py(:);
    end

    % ---------- Panel 1: heat on/off ----------
    ax1 = subplot(1,2,1); hold on;
    heat_cmd = a2p_data.dq(1).heat_active(:);   % 0/1 command on the imaging grid
    heat_on  = heat_cmd > 0;
    heat_off = ~heat_on;

    h_off = plot(x_vm(heat_off), y_vm(heat_off), 'k.', 'MarkerSize', 10);
    h_on  = plot(x_vm(heat_on),  y_vm(heat_on),  'r.', 'MarkerSize', 10);
    plot(x_vm(1), y_vm(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

    legend([h_off, h_on], {'Heat off', 'Heat on'}, 'Location', 'best');
    xlabel('X position (mm)'); ylabel('Y position (mm)');
    title('Flat path — heat on (red) vs off (black)');
    axis equal; grid on; hold off;

    % ---------- Panel 2: dF/F-colored path (right side) ----------
    ax2 = subplot(1,2,2); hold on;
    dff_on_vm = dff_right;
    valid  = ~isnan(dff_on_vm);
    scatter(x_vm(valid), y_vm(valid), 10, dff_on_vm(valid), 'filled', 'MarkerEdgeColor', 'none');
    axis equal;
    colormap(jet);
    c = colorbar; ylabel(c, 'dF/F (right)');
    plot(x_vm(1), y_vm(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    xlabel('X position (mm)'); ylabel('Y position (mm)');
    title('Flat path colored by dF/F (right)');
    grid on; hold off;

    linkaxes([ax1, ax2], 'xy');
    save_plot_with_title_as_filename('heat_and_dff_colored', 'flatpath', savepath);
end

%% ===== CROSS-CORRELATIONS: each dF/F vs each kinematic variable =====
% Per-lag Pearson correlation of dF/F(t+lag) with behavior(t).
% Positive lag => dF/F lags behavior (behavior precedes dF/F).
% For heat experiments the correlation is computed separately over heat-ON
% and heat-OFF time points so the two relationships can be compared.
maxlag_s  = 5;                              % +/- lag window in seconds
maxlag_N  = round(maxlag_s * sample_rate);
lags_samp = -maxlag_N:maxlag_N;
lags_s    = lags_samp / sample_rate;

for di = 1:numel(dff_vars)
    figure('Name', sprintf('Cross-correlation: %s', dff_names{di}));
    n_kin = numel(kin_vars);
    for ki = 1:n_kin
        d = dff_vars{di};
        k = kin_vars{ki};
        valid = ~isnan(d) & ~isnan(k);

        subplot(n_kin, 1, ki); hold on;
        if heat
            xc_off = lagged_corr(d, k, valid & ~heat_mask, lags_samp);
            xc_on  = lagged_corr(d, k, valid &  heat_mask, lags_samp);
            plot(lags_s, xc_off, 'Color', [0 0 0],   'LineWidth', 1.5);
            plot(lags_s, xc_on,  'Color', [0.85 0 0], 'LineWidth', 1.5);
            legend({'heat off', 'heat on'}, 'Location', 'best');
        else
            xc = lagged_corr(d, k, valid, lags_samp);
            plot(lags_s, xc, 'Color', [0 0 0], 'LineWidth', 1.5);
        end
        line([0 0], ylim, 'Color', 'k', 'LineStyle', ':');
        xlabel('Lag (s)   [+lag: dF/F lags behavior]');
        ylabel('corr coeff');
        title(sprintf('%s  \\times  %s', dff_names{di}, kin_names{ki}));
        grid on; hold off;
    end
    sgtitle(sprintf('Cross-correlation of %s with kinematics', dff_names{di}));
    save_plot_with_title_as_filename('xcorr', sprintf('dff_%s', dff_tags{di}), savepath);
end

%% ===== SCATTER PLOTS: each dF/F vs each kinematic variable =====
for di = 1:numel(dff_vars)
    figure('Name', sprintf('Scatter: %s', dff_names{di}));
    n_kin = numel(kin_vars);
    for ki = 1:n_kin
        d = dff_vars{di};
        k = kin_vars{ki};
        good = ~isnan(d) & ~isnan(k);

        subplot(2, 2, ki); hold on;
        if heat
            % Split points by heat state, fit + annotate each separately
            off = good & ~heat_mask;
            on  = good &  heat_mask;
            [p_off, r_off] = fit_and_r(k(off), d(off));
            [p_on,  r_on ] = fit_and_r(k(on),  d(on));

            scatter(k(off), d(off), 8, [0 0 0],   'filled', 'MarkerFaceAlpha', 0.20);
            scatter(k(on),  d(on),  8, [0.85 0 0], 'filled', 'MarkerFaceAlpha', 0.20);

            xr = [min(k(good)) max(k(good))];
            xfit = linspace(xr(1), xr(2), 100);
            h_off = plot(xfit, polyval(p_off, xfit), 'Color', [0 0 0],   'LineWidth', 1.8);
            h_on  = plot(xfit, polyval(p_on,  xfit), 'Color', [0.85 0 0], 'LineWidth', 1.8);
            legend([h_off, h_on], ...
                {sprintf('heat off (r = %.2f)', r_off), ...
                 sprintf('heat on (r = %.2f)',  r_on)}, 'Location', 'best');
            title(sprintf('%s vs %s', dff_names{di}, kin_tags{ki}));
        else
            [p, r] = fit_and_r(k(good), d(good));
            scatter(k(good), d(good), 8, 'filled', 'MarkerFaceAlpha', 0.25);
            xfit = linspace(min(k(good)), max(k(good)), 100);
            plot(xfit, polyval(p, xfit), 'r-', 'LineWidth', 1.5);
            title(sprintf('%s vs %s (r = %.2f)', dff_names{di}, kin_tags{ki}, r));
        end
        xlabel(kin_names{ki});
        ylabel(dff_names{di});
        grid on; hold off;
    end
    sgtitle(sprintf('%s vs kinematics', dff_names{di}));
    save_plot_with_title_as_filename('scatter', sprintf('dff_%s', dff_tags{di}), savepath);
end

%% ===== local helper functions =====
function xc = lagged_corr(d, k, mask, lags)
% Per-lag Pearson correlation of d(i+lag) with k(i), using only time points
% where mask is true for both i and i+lag. Matches xcorr(d,k) lag sign
% convention; returns NaN for lags with too few valid samples.
    N  = numel(d);
    xc = nan(size(lags));
    for li = 1:numel(lags)
        L  = lags(li);
        i0 = max(1, 1 - L);
        i1 = min(N, N - L);
        if i1 < i0, continue; end
        idx = (i0:i1).';
        a = d(idx + L);
        b = k(idx);
        m = mask(idx) & mask(idx + L);
        a = a(m); b = b(m);
        if numel(a) < 10 || std(a) == 0 || std(b) == 0, continue; end
        R = corrcoef(a, b);
        xc(li) = R(1,2);
    end
end

function [p, r] = fit_and_r(x, y)
% Least-squares line p = [slope intercept] and Pearson r for x,y.
% Returns NaNs if there are too few points or no variance.
    x = x(:); y = y(:);
    if numel(x) < 2 || std(x) == 0 || std(y) == 0
        p = [NaN NaN]; r = NaN; return;
    end
    p = polyfit(x, y, 1);
    R = corrcoef(x, y);
    r = R(1,2);
end
