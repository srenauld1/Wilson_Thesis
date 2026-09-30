%% ANALYZE HEAT BARRIER DATA (new heat / vmData format)
% May 2026

%% SETUP

baseFolder = '/Users/sophiarenauld/stacks/';  % change if needed

%% LOAD DATA

% For batch use, callers may define heatTrialFolder before running this
% script. Interactive use retains the three file pickers.
if exist('heatTrialFolder', 'var') && ~isempty(heatTrialFolder)
    if ~isfolder(heatTrialFolder)
        error('heat_master_analysis:TrialFolderNotFound', ...
            'Heat trial folder does not exist: %s', heatTrialFolder);
    end
    ftFile = find_single_file(heatTrialFolder, '*_ficTracData.mat');
    vmFile = find_single_file(heatTrialFolder, '*virmenData.mat');
    imFile = find_single_file(heatTrialFolder, '*_or_z0z1_s_.mat');
    savepath = fullfile(heatTrialFolder, 'smr_analysis');
else
    ftFile = uigetfile_full('*_ficTracData.mat', 'Select FicTrac Data', baseFolder);
    if isempty(ftFile), return; end
    vmFile = uigetfile_full('*virmenData.mat', 'Select Virmen Data', baseFolder);
    if isempty(vmFile), return; end
    imFile = uigetfile_full('*_or_z0z1_s_.mat', 'Select IMAGING Data', baseFolder);
    if isempty(imFile), return; end
    savepath = fullfile(fileparts(imFile), 'smr_analysis');
end

if ~isfolder(savepath), mkdir(savepath); end
disp(['Heat analysis output: ' savepath]);

load(ftFile, 'ftData');
load(vmFile);   % assumes vmData and daqVolumeTime exist
imaging_data = load(imFile);   % keep as struct

%% COMMON LABELS (dF/F + kinematics on common time bases)

% Imaging time & dF/F
time_dff = daqVolumeTime.daqVolumeTime{1,1};      % imaging time (s, numeric)
dff      = imaging_data.roi.ts;                   % dF/F

% Kinematic time from vmData.Seconds (duration or cell)
t_raw    = vmData.Seconds{1,1};                   % usually duration
time_kin = seconds(t_raw);                        % numeric seconds
time_kin = time_kin(:);                           % column

% Kinematic signals (Virmen)
fwd_vel = smooth(vmData.vel_for{1,1});            % forward velocity
rot_vel = smooth(vmData.vel_yaw{1,1});            % yaw / rotational velocity

% Ensure column vectors
time_dff = time_dff(:);
dff      = dff(:);
fwd_vel  = fwd_vel(:);
rot_vel  = rot_vel(:);

% Clip dF/F to the kinematic time range
tmax_kin         = max(time_kin);
mask_dff         = time_dff <= tmax_kin;
time_dff_clipped = time_dff(mask_dff);
dff_clipped      = dff(mask_dff);

% Interpolate dF/F onto the kinematic time base
fluo = interp1(time_dff_clipped, dff_clipped, time_kin, 'linear', NaN);  % dF/F on time_kin
t    = time_kin;                                                         % unified time
fwd  = fwd_vel;
rot  = rot_vel;

% Motion (binary) from forward speed
vel_thresh = 1;                 % mm/s threshold – adjust as needed
motion     = abs(fwd) > vel_thresh;

% Sample rates
dt_kin        = median(diff(time_kin));      % behavior time step
sample_rate_kin = 1 / dt_kin;                % Hz (behavior)
dt_im         = median(diff(time_dff_clipped));  % imaging dt
sample_rate_im  = 1 / dt_im;                 % Hz (imaging)

%% QUICK CHECK PLOTS (OPTIONAL)

figure(1); clf;
plot(time_kin, fwd_vel);
xlabel('Time (s)');
ylabel('Forward velocity (vmData)');
title('Forward velocity (smoothed)');

figure(2); clf;
plot(ftData.pos_for{1,1}, ftData.pos_side{1,1});
xlabel('Forward position');
ylabel('Side position');
title('FicTrac trajectory');

figure(3); clf; hold on;
plot(time_dff_clipped, dff_clipped*5, 'b');
plot(ftData.seconds{1,1}, ftData.fwVel{1,1}, 'k');
plot(ftData.seconds{1,1}, ftData.yawVel{1,1}, 'r');
xlabel('Time (s)');
ylabel('Signal');
legend({'dF/F (×5)','fwVel (FicTrac)','yawVel (FicTrac)'});
title('Imaging vs FicTrac signals');
grid on; hold off;

%% FLAT PATH: heat on/off AND dF/F, linked axes

figure(10); clf;

% Common Virmen position & time
x_vm = vmData.pos_x{1,1}(:);
y_vm = vmData.pos_y{1,1}(:);
t_vm = seconds(vmData.Seconds{1,1});   % duration -> double (s)
t_vm = t_vm(:);

% ---------- Panel 1: heat on/off ----------
ax1 = subplot(1,2,1); hold on;

% Heat may be stored as cell in table — extract to numeric vector
heat = vmData.heat;
if iscell(heat), heat = cell2mat(heat); end
heat = heat(:);

heat_on  = heat > 0.0001;
heat_off = ~heat_on;

h_off = plot(x_vm(heat_off), y_vm(heat_off), 'k.', 'MarkerSize', 10);
h_on  = plot(x_vm(heat_on),  y_vm(heat_on),  'r.', 'MarkerSize', 10);
plot(x_vm(1), y_vm(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % start point

legend([h_off, h_on], {'Heat off', 'Heat on'}, 'Location', 'best');
xlabel('X position');
ylabel('Y position');
title('Flat path — heat on (red) vs off (black)');
axis equal; grid on;
hold off;

% ---------- Panel 2: dF/F-colored path ----------
ax2 = subplot(1,2,2); hold on;

% Imaging dF/F and time (already clipped to overlap)
t_dff   = time_dff_clipped(:);           % same length as dff_clipped
dff_vec = dff_clipped(:);

% Interpolate dF/F onto Virmen time base
dff_on_vm = interp1(t_dff, dff_vec, t_vm, 'linear', NaN);

valid  = ~isnan(dff_on_vm);
x_plot = x_vm(valid);
y_plot = y_vm(valid);
c_plot = dff_on_vm(valid);

scatter(x_plot, y_plot, 10, c_plot, 'filled', 'MarkerEdgeColor', 'none');
axis equal;
colormap(jet);
c = colorbar;
ylabel(c, 'dF/F');

plot(x_vm(1), y_vm(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % start point
xlabel('X position');
ylabel('Y position');
title('Flat path colored by dF/F');
grid on;
hold off;

% ---------- Link axes so zooming/panning one panel affects the other ----------
linkaxes([ax1, ax2], 'xy');

% ---------- Save figure ----------
save_plot_with_title_as_filename('heat_and_dff_colored', 'flatpath', savepath);

%% CALCULATE MOVING VS NOT
vel_on = 1;
vel_off = 0.5;

% Start state (assume not moving)
if fwd_vel(1) >= vel_on
    motion(1) = 1;
else
    motion(1) = 0;
end

% Apply Schmitt trigger logic
for i = 2:numel(fwd_vel)
    if motion(i-1) == 0 && fwd_vel(i) >= vel_on
        motion(i) = 1; % start moving
    elseif motion(i-1) == 1 && fwd_vel(i) > vel_off
        motion(i) = 1; % keep moving
    else
        motion(i) = 0; % stop moving
    end
end


%% aligned to motion OFFSET, plotting dff and fwd
preN = 10*round(sample_rate_kin);    % how many samples to show before motion offset
min_run_len = 3*round(sample_rate_kin); % 1 second off
offN = 100;     % set >0 if you want more after next offset



%% parameters in SECONDS
pre_run_min_s   = 0.5;    % must be moving at least this long BEFORE offset
off_min_s       = 0.1;   % min time "off" after offset
off_max_s       = 3;    % max time "off" after offset

dt_kin          = median(diff(t));
sample_rate_kin = 1/dt_kin;

pre_run_min_N = round(pre_run_min_s * sample_rate_kin);
off_min_N     = round(off_min_s      * sample_rate_kin);
off_max_N     = round(off_max_s      * sample_rate_kin);

%% find onsets/offsets as before
d = diff([0; motion(:); 0]);
onsets  = find(d == 1);
offsets = find(d == -1) - 1;    % last moving sample

%% for each offset, find next onset (end of off period)
my_onsets = nan(size(offsets));   % index of next onset (or end+1)
for i = 1:numel(offsets)
    next_on = find(onsets > offsets(i), 1, 'first');
    if ~isempty(next_on)
        my_onsets(i) = onsets(next_on);
    else
        my_onsets(i) = numel(motion) + 1;  % off period runs to end
    end
end

off_lens = my_onsets - offsets;   % off-period length (in samples)

%% condition 1: fly has been moving >= pre_run_min_s before offset
pre_run_ok = false(size(offsets));
for i = 1:numel(offsets)
    off_idx = offsets(i);
    start_idx = max(1, off_idx - pre_run_min_N + 1);  % inclusive window
    idx_range = start_idx:off_idx;
    % Require motion == 1 for the entire pre_run window
    pre_run_ok(i) = all(motion(idx_range));
end

%% condition 2: off period between 0.05s and 2s
off_ok = (off_lens >= off_min_N) & (off_lens <= off_max_N);

%% keep only offsets that satisfy BOTH conditions
keep = pre_run_ok & off_ok;

offsets_kept  = offsets(keep);
my_onsets_kept = my_onsets(keep);

%% SEPARATE VARIABLES
fluo_cell = {}; behav_cell = {}; rot_cell = {}; heat_cell = {}; t_cell = {};

for i = 1:numel(offsets_kept)
    t_start = offsets_kept(i) - preN + 1;
    if t_start < 1, continue; end
    t_end = my_onsets_kept(i) - 1 + offN;
    t_end = min(t_end,numel(fluo));
    idx_epoch = t_start:t_end;
    t_offset = t(offsets_kept(i));
    t_rel = t(idx_epoch) - t_offset;

    fluo_cell{end+1}  = fluo(idx_epoch);
    behav_cell{end+1} = fwd(idx_epoch);
    rot_cell{end+1}   = abs(rot_vel(idx_epoch));
    heat_cell{end+1}  = heat(idx_epoch);      % <-- new
    t_cell{end+1}     = t_rel;
end

% Only keep epochs that have enough pre points and motion off period
len_epoch = cellfun(@length, fluo_cell);
minlen    = min(len_epoch);  % For mean, only use shared window

% Truncate/crop for average
fluo_mat  = nan(minlen, numel(fluo_cell));
behav_mat = nan(minlen, numel(behav_cell));
rot_mat   = nan(minlen, numel(rot_cell));
t_mat     = nan(minlen, numel(fluo_cell));
for i = 1:numel(fluo_cell)
    fluo_mat(:,i)  = fluo_cell{i}(1:minlen);
    behav_mat(:,i) = behav_cell{i}(1:minlen);
    rot_mat(:,i)   = rot_cell{i}(1:minlen);
    t_mat(:,i)     = t_cell{i}(1:minlen);
end
t_mean = mean(t_mat,2,'omitnan');



%% plot all

figure;

ax1 = subplot(3,1,1); hold on;
for i = 1:numel(fluo_cell)
    plot(t_cell{i}, fluo_cell{i}, 'Color', [0.2549 0.4118 0.8824], 'HandleVisibility','off');
end
plot(t_mean, nanmean(fluo_mat,2), 'Color', [0 0 0.5451], 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement offset (s)');
ylabel('DF/F');
title('dF/F aligned to movement offset');
hold off;

ax2 = subplot(3,1,2); hold on;
for i = 1:numel(behav_cell)
    plot(t_cell{i}, behav_cell{i}, 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
end
plot(t_mean, nanmean(behav_mat,2), 'k', 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement offset (s)');
ylabel('Forward velocity');
title('Forward velocity aligned to movement offset');
hold off;

ax3 = subplot(3,1,3); hold on;
for i = 1:numel(rot_cell)
    plot(t_cell{i}, rot_cell{i}, 'Color', [0.7 0 0], 'HandleVisibility','off');
end
plot(t_mean, nanmean(rot_mat,2), 'LineWidth', 2, 'Color', [0.3 0 0], 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement offset (s)');
ylabel('Rotational Speed (deg/s)');
title('Rotational speed aligned to movement offset');
hold off;

n_epochs = numel(fluo_cell);
subplot(ax1);   % go back to top subplot
yl = ylim; xl = xlim;
text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'EdgeColor','k');

% Link x-axes of the three subplots
linkaxes([ax1, ax2, ax3], 'x');


save_plot_with_title_as_filename('movement_offset', 'variables', savepath);

%% --- classify epochs by rotational speed around offset (-0.5 to +0.5 s) ---

high_thresh = 1;   % deg/s, example: mean rot speed >= this = "HIGH"
low_thresh  = 0.8;   % deg/s, example: mean rot speed <= this = "LOW"

high_idx = [];   % indices into fluo_cell/... that are HIGH
low_idx  = [];   % indices that are LOW

for i = 1:numel(rot_cell)
    t_rel = t_cell{i};
    r     = rot_cell{i};      % already abs(rot_vel)

    % indices in [-0.5, +0.5] s
    mask_win = (t_rel >= -0.5) & (t_rel <= 0.5);

    if ~any(mask_win)
        continue;
    end

    mean_rot = mean(r(mask_win), 'omitnan');

    if mean_rot >= high_thresh
        high_idx(end+1) = i; 
    elseif mean_rot <= low_thresh
        low_idx(end+1)  = i;
    end
end


%% plot high rot and low rot
plot_epoch_group(high_idx, fluo_cell, behav_cell, rot_cell, t_cell, savepath, ...
                 sprintf('HIGH rot (>= %.2f deg/s)', high_thresh));

plot_epoch_group(low_idx, fluo_cell, behav_cell, rot_cell, t_cell, savepath, ...
                 sprintf('LOW rot (<= %.2f deg/s)', low_thresh));

%% --- classify epochs by HEAT around offset (-0.5 to +0.5 s) ---

high_heat_thresh = 0.05;   % example units: same as vmData.heat
low_heat_thresh  = 0.0001;   % set these as you like

high_idx = [];   % epochs with high heat
low_idx  = [];   % epochs with low heat

for i = 1:numel(heat_cell)
    t_rel = t_cell{i};
    h     = heat_cell{i};

    % indices in [-0.5, +0.5] s
    mask_win = (t_rel >= -0.5) & (t_rel <= 0.5);
    if ~any(mask_win)
        continue;
    end

    mean_heat = mean(h(mask_win), 'omitnan');

    if mean_heat >= high_heat_thresh
        high_idx(end+1) = i; %#ok<AGROW>
    elseif mean_heat <= low_heat_thresh
        low_idx(end+1)  = i; %#ok<AGROW>
    end
end

%% plot heat groups
plot_epoch_group(high_idx, fluo_cell, behav_cell, rot_cell, t_cell, savepath, ...
    sprintf('HIGH heat (mean heat >= %.2f)', high_heat_thresh));

plot_epoch_group(low_idx, fluo_cell, behav_cell, rot_cell, t_cell, savepath, ...
    sprintf('LOW heat (mean heat <= %.2f)', low_heat_thresh));

%% dff peaks
% Match the VES041 peak definition used by plotting_ves041: trial-specific
% 3 sigma prominence, >=0.25 s width, and >=2 s separation. No upper dF/F
% cutoff is applied. Fill only interpolation-edge gaps for peak detection.
fluo_for_peaks = fillmissing(fluo, 'linear', 'EndValues', 'nearest');
baseline_N = max(3, round(30 * sample_rate_kin));
dff_baseline = movmedian(fluo_for_peaks, baseline_N, 'omitmissing');
noise_sigma = 1.4826 * mad(fluo_for_peaks - dff_baseline, 1);
peak_prominence = 3 * noise_sigma;
peak_args = {'MinPeakProminence', peak_prominence, ...
    'MinPeakDistance', round(2 * sample_rate_kin), ...
    'MinPeakWidth', round(0.25 * sample_rate_kin)};
figure
findpeaks(fluo_for_peaks, peak_args{:});
title(sprintf('VES041 peaks: prominence >= 3 sigma (%.3g dF/F)', peak_prominence));
[pks, locs] = findpeaks(fluo_for_peaks, peak_args{:});

window_s = 2; % seconds before/after
window_N = round(window_s * sample_rate_kin);

% Ensure column vectors
dff  = fluo(:);
bvf  = fwd(:);
rots = abs(rot(:));         % rotational speed (deg/s)
t    = t(:);

% heat on same time base t
heat = vmData.heat;
if iscell(heat), heat = cell2mat(heat); end
heat = heat(:);

dff_epochs  = {};
bvf_epochs  = {};
rots_epochs = {};
heat_epochs = {};
t_epochs    = {};

for i = 1:numel(locs)
    idx_c    = locs(i);
    idx_start = idx_c - window_N;
    idx_end   = idx_c + window_N;
    if idx_start < 1 || idx_end > length(dff), continue, end

    idx_win = idx_start:idx_end;
    t0      = t(idx_c);
    t_epoch = t(idx_win) - t0;

    dff_epochs{end+1}   = dff(idx_win);
    bvf_epochs{end+1}   = bvf(idx_win);
    rots_epochs{end+1}  = rots(idx_win);
    heat_epochs{end+1}  = heat(idx_win);   % <- parallel heat
    t_epochs{end+1}     = t_epoch;
end

n_epochs = numel(bvf_epochs);
if n_epochs < 1
    error('No useable peaks found with sufficient window length!');
end

dt_median = median(cellfun(@(tt) median(diff(tt)), t_epochs));
t_mean = -window_N*dt_median : dt_median : window_N*dt_median;

dff_mat  = nan(length(t_mean), n_epochs);
bvf_mat  = nan(length(t_mean), n_epochs);
rots_mat = nan(length(t_mean), n_epochs);

for i = 1:n_epochs
    dff_mat(:,i)  = interp1(t_epochs{i}, dff_epochs{i},  t_mean, 'linear', nan);
    bvf_mat(:,i)  = interp1(t_epochs{i}, bvf_epochs{i}, t_mean, 'linear', nan);
    rots_mat(:,i) = interp1(t_epochs{i}, rots_epochs{i},t_mean, 'linear', nan);
end

%% --- classify peaks by HEAT around the dF/F peak (-0.5 to +0.5 s) ---

heat_thresh = 0.005;   % define "heat on" as mean heat >= this (change as needed)

heat_idx    = [];    % epochs with heat on
noheat_idx  = [];    % epochs with heat off

for i = 1:n_epochs
    t_rel = t_epochs{i};
    h     = heat_epochs{i};

    % indices in [-0.5, +0.5] s around the dF/F peak
    mask_win = (t_rel >= -0.5) & (t_rel <= 0.5);

    if ~any(mask_win)
        continue;
    end

    mean_heat = mean(h(mask_win), 'omitnan');

    if mean_heat >= heat_thresh
        heat_idx(end+1)   = i; %#ok<AGROW>
    else
        noheat_idx(end+1) = i; %#ok<AGROW>
    end
end


%% plot: ALL, HEAT, NO-HEAT

% All peaks (your original plot)
plot_peak_group(1:n_epochs, t_epochs, dff_epochs, bvf_epochs, rots_epochs, ...
                t_mean, dff_mat, bvf_mat, rots_mat, 'ALL peaks', savepath);

% Heat peaks
plot_peak_group(heat_idx, t_epochs, dff_epochs, bvf_epochs, rots_epochs, ...
                t_mean, dff_mat, bvf_mat, rots_mat, ...
                sprintf('HEAT (mean heat >= %.2f)', heat_thresh), savepath);

% Non-heat peaks
plot_peak_group(noheat_idx, t_epochs, dff_epochs, bvf_epochs, rots_epochs, ...
                t_mean, dff_mat, bvf_mat, rots_mat, ...
                sprintf('NO-HEAT (mean heat < %.2f)', heat_thresh), savepath);

%% HELPER FUNCTION: safe uigetfile + full path
function fullFilePath = uigetfile_full(pattern, prompt, baseFolder)
    [fileName, pathName] = uigetfile(pattern, prompt, baseFolder);
    if isequal(fileName, 0)
        disp(['User canceled: ' prompt]);
        fullFilePath = [];
    else
        fullFilePath = fullfile(pathName, fileName);
        disp(['Loaded: ' fullFilePath]);
    end
end

function fullFilePath = find_single_file(folder, pattern)
    matches = dir(fullfile(folder, pattern));
    matches = matches(~[matches.isdir]);
    if numel(matches) ~= 1
        error('heat_master_analysis:InputFileCount', ...
            'Expected one %s file in %s; found %d.', ...
            pattern, folder, numel(matches));
    end
    fullFilePath = fullfile(matches(1).folder, matches(1).name);
    disp(['Loaded: ' fullFilePath]);
end

%% helper function
function plot_epoch_group(idx_set, fluo_cell, behav_cell, rot_cell, t_cell, savepath, group_name)

    if isempty(idx_set)
        warning('No epochs in group: %s', group_name);
        return;
    end

    % Build matrices on a common time base for this group
    len_epoch = cellfun(@length, fluo_cell(idx_set));
    minlen    = min(len_epoch);
    
    fluo_mat  = nan(minlen, numel(idx_set));
    behav_mat = nan(minlen, numel(idx_set));
    rot_mat   = nan(minlen, numel(idx_set));
    t_mat     = nan(minlen, numel(idx_set));

    for k = 1:numel(idx_set)
        i = idx_set(k);
        fluo_mat(:,k)  = fluo_cell{i}(1:minlen);
        behav_mat(:,k) = behav_cell{i}(1:minlen);
        rot_mat(:,k)   = rot_cell{i}(1:minlen);
        t_mat(:,k)     = t_cell{i}(1:minlen);
    end

    t_mean = mean(t_mat, 2, 'omitnan');

    figure;

    ax1 = subplot(3,1,1); hold on;
    for k = 1:numel(idx_set)
        plot(t_cell{idx_set(k)}, fluo_cell{idx_set(k)}, ...
             'Color', [0.2549 0.4118 0.8824], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(fluo_mat,2), 'Color', [0 0 0.5451], 'LineWidth', 2);
    line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
    xlabel('Time from movement offset (s)');
    ylabel('dF/F');
    title(sprintf('dF/F aligned to offset (%s)', group_name));
    hold off;

    ax2 = subplot(3,1,2); hold on;
    for k = 1:numel(idx_set)
        plot(t_cell{idx_set(k)}, behav_cell{idx_set(k)}, ...
             'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(behav_mat,2), 'k', 'LineWidth', 2);
    line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
    xlabel('Time from movement offset (s)');
    ylabel('Forward velocity');
    title(sprintf('Forward velocity (%s)', group_name));
    hold off;

    ax3 = subplot(3,1,3); hold on;
    for k = 1:numel(idx_set)
        plot(t_cell{idx_set(k)}, rot_cell{idx_set(k)}, ...
             'Color', [0.7 0 0], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(rot_mat,2), 'Color', [0.3 0 0], 'LineWidth', 2);
    line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
    xlabel('Time from movement offset (s)');
    ylabel('Rotational Speed (deg/s)');
    title(sprintf('Rotational speed (%s)', group_name));
    hold off;

    linkaxes([ax1, ax2, ax3], 'x');

    n_epochs = numel(idx_set);
    subplot(ax1);
    yl = ylim; xl = xlim;
    text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

    sgtitle(sprintf('Offset-aligned epochs (%s)', group_name));

    if nargin >= 7 && ~isempty(savepath)
        fname = sprintf('movement_offset_%s', regexprep(group_name,'\W+','_'));
        save_plot_with_title_as_filename(fname, 'variables', savepath);
    end
end

%% helper to plot a subset of epochs
function plot_peak_group(idx_set, t_epochs, dff_epochs, bvf_epochs, rots_epochs, ...
                         t_mean, dff_mat, bvf_mat, rots_mat, label_str, savepath)

    if isempty(idx_set)
        warning('No epochs in group: %s', label_str);
        return;
    end

    % Build group-specific mean matrices
    grp_dff  = nan(size(dff_mat,1), numel(idx_set));
    grp_bvf  = nan(size(bvf_mat,1), numel(idx_set));
    grp_rots = nan(size(rots_mat,1), numel(idx_set));
    for k = 1:numel(idx_set)
        grp_dff(:,k)  = dff_mat(:, idx_set(k));
        grp_bvf(:,k)  = bvf_mat(:, idx_set(k));
        grp_rots(:,k) = rots_mat(:, idx_set(k));
    end

    figure;

    % dF/F
    ax1 = subplot(3,1,1); hold on;
    for k = 1:numel(idx_set)
        i = idx_set(k);
        plot(t_epochs{i}, dff_epochs{i}, 'Color',[0.4 0.8 1 0.2],'HandleVisibility','off');
    end
    plot(t_mean, nanmean(grp_dff,2), 'b', 'LineWidth', 2);
    line([0 0], ylim, 'Color','k','LineStyle',':','LineWidth',2);
    xlabel('Time from dF/F peak (s)');
    ylabel('dF/F');
    title(sprintf('dF/F aligned to dF/F peak (%s)', label_str));
    hold off;

    % Forward velocity
    ax2 = subplot(3,1,2); hold on;
    for k = 1:numel(idx_set)
        i = idx_set(k);
        plot(t_epochs{i}, bvf_epochs{i}, 'Color',[0.5 0.5 0.5],'HandleVisibility','off');
    end
    plot(t_mean, nanmean(grp_bvf,2), 'k','LineWidth',2);
    line([0 0], ylim, 'Color','k','LineStyle',':','LineWidth',2);
    xlabel('Time from dF/F peak (s)');
    ylabel('Forward velocity');
    title(sprintf('Forward velocity aligned to dF/F peak (%s)', label_str));
    hold off;

    % Rotational velocity
    ax3 = subplot(3,1,3); hold on;
    for k = 1:numel(idx_set)
        i = idx_set(k);
        plot(t_epochs{i}, rots_epochs{i}, 'Color',[1 0.5 0.5 0.2],'HandleVisibility','off');
    end
    plot(t_mean, nanmean(grp_rots,2), 'r','LineWidth',2);
    line([0 0], ylim, 'Color','k','LineStyle',':','LineWidth',2);
    xlabel('Time from dF/F peak (s)');
    ylabel('Rotational velocity (deg/s)');
    title(sprintf('Rotational velocity aligned to dF/F peak (%s)', label_str));
    hold off;

    linkaxes([ax1, ax2, ax3], 'x');

    n_grp = numel(idx_set);
    subplot(ax1);
    yl = ylim; xl = xlim;
    text(xl(2), yl(2), sprintf('n peaks = %d', n_grp), ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

    sgtitle(sprintf('dF/F and behavior aligned to dF/F peaks (%s)', label_str));

    if ~isempty(savepath)
        fname = sprintf('dff_peaks_%s', regexprep(label_str,'\W+','_'));
        save_plot_with_title_as_filename(fname, 'kinematics', savepath);
    end
end
