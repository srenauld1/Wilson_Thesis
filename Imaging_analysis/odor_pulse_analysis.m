function odor_pulse_analysis(time, fwd_vel, angular_vel, dff, odor, savepath)
% ODOR_PULSE_ANALYSIS
% Chop data into odor-on pulses with pre/post windows, compute averages,
% and plot fwd vel, rot vel, dF/F, and odor on 4 subplots.
%
% Inputs:
%   time        : [N x 1] time vector (s)
%   fwd_vel     : [N x 1] forward velocity
%   angular_vel : [N x 1] angular velocity
%   dff         : [N x 1] dF/F trace
%   odor        : [N x 1] odor trace
%   savepath    : directory to save the figure

%% ---------------- USER PARAMETERS ----------------
odor_on_thr = -0.5;     % odor ON if odor > this
odor_off_thr = -1;      % odor OFF if odor < this (optional use)
pre_win_sec  = 4;      % seconds before odor onset
post_win_sec = 4;      % seconds after odor offset
%% -------------------------------------------------

% Ensure column vectors
time        = time(:);
fwd_vel     = fwd_vel(:);
angular_vel = angular_vel(:);
dff         = dff(:);
odor        = odor(:);

N  = numel(time);
dt = median(diff(time));
rate = 1/dt;

pre_pts  = round(pre_win_sec  * rate);
post_pts = round(post_win_sec * rate);

%% ---------- Binary odor-on vector ----------
odor_on = odor > odor_on_thr;  % 1 during odor ON
% Optional stricter OFF:
% odor_on(odor < odor_off_thr) = 0;

change_points = [1; find(diff(odor_on) ~= 0) + 1; N];

%% ---------- Containers ----------
fwd_ext_cells  = {};
rot_ext_cells  = {};
dff_ext_cells  = {};
odor_ext_cells = {};

odor_on_markers  = [];  % onset index relative to extended segment
odor_off_markers = [];

%% ---------- Loop over segments ----------
for k = 1:length(change_points)-1
    seg_start = change_points(k);
    seg_end   = change_points(k+1) - 1;

    if odor_on(seg_start)  % odor-on pulse
        pulse_idx = seg_start:seg_end;

        % Pre / Post indices
        pre_start = max(1, seg_start - pre_pts);
        post_end  = min(N, seg_end + post_pts);

        ext_idx = pre_start:post_end;

        % Store extended segments
        fwd_ext_cells{end+1}  = fwd_vel(ext_idx);
        rot_ext_cells{end+1}  = angular_vel(ext_idx);
        dff_ext_cells{end+1}  = dff(ext_idx);
        odor_ext_cells{end+1} = odor(ext_idx);

        onset_rel  = seg_start - pre_start + 1;
        offset_rel = seg_end   - pre_start + 1;
        odor_on_markers(end+1)  = onset_rel;
        odor_off_markers(end+1) = offset_rel;
    end
end

if isempty(fwd_ext_cells)
    warning('No odor pulses meeting criteria.');
    return;
end

%% ---------- Length-based filtering ----------
lens      = cellfun(@length, fwd_ext_cells);
mean_len  = mean(lens);
cutoff_len = mean_len - 10;   % 10 samples below mean length

keep_mask = lens >= cutoff_len;

fwd_ext_cells  = fwd_ext_cells(keep_mask);
rot_ext_cells  = rot_ext_cells(keep_mask);
dff_ext_cells  = dff_ext_cells(keep_mask);
odor_ext_cells = odor_ext_cells(keep_mask);
odor_on_markers  = odor_on_markers(keep_mask);
odor_off_markers = odor_off_markers(keep_mask);

if isempty(fwd_ext_cells)
    warning('No odor pulses left after length-based filtering.');
    return;
end

%% ---------- Trim all segments to common length ----------
min_len = min(cellfun(@length, fwd_ext_cells));

trim = @(C) cellfun(@(x) x(1:min_len), C, 'UniformOutput', false);
fwd_trim  = trim(fwd_ext_cells);
rot_trim  = trim(rot_ext_cells);
dff_trim  = trim(dff_ext_cells);
odor_trim = trim(odor_ext_cells);

nPulses = numel(fwd_trim);
L       = min_len;

fwd_mat  = nan(L, nPulses);
rot_mat  = nan(L, nPulses);
dff_mat  = nan(L, nPulses);
odor_mat = nan(L, nPulses);

for p = 1:nPulses
    fwd_mat(:,p)  = fwd_trim{p};
    rot_mat(:,p)  = rot_trim{p};
    dff_mat(:,p)  = dff_trim{p};
    odor_mat(:,p) = odor_trim{p};
end

% Now each row is a time point, each column a pulse
mean_fwd  = mean(fwd_mat,  2, 'omitnan');   % L×1
mean_rot  = mean(rot_mat,  2, 'omitnan');
mean_dff  = mean(dff_mat,  2, 'omitnan');
mean_odor = mean(odor_mat, 2, 'omitnan');


% Relative time axis (centered roughly at 0 for an onset)
% this uses 2nd pulse markers; you may want index 1
t_rel = ((0:min_len-1) - odor_on_markers(2) + 1) / rate;
onset_t_rel  = 0;
offset_t_rel = (odor_off_markers(2) - odor_on_markers(2)) / rate;

%% ---------- Plot 4-subplot figure with average line ----------
figure;
ax = gobjects(4,1);

% 1) Forward velocity
ax(1) = subplot(4,1,1); hold on;
for i = 1:numel(fwd_trim)
    plot(t_rel, fwd_trim{i}, 'Color', [0.7 1 0.7], 'LineWidth', 1);
end
xline(onset_t_rel,  'k--', 'LineWidth', 1.5);
xline(offset_t_rel, 'k--', 'LineWidth', 1.5);
plot(t_rel, mean_fwd, 'Color',[0 0.5 0], 'LineWidth', 2);   % bold mean
ylabel('Fwd vel (mm/s)');
title('Odor pulse: Forward velocity');
grid on;
hold off;

% 2) Rotational velocity
ax(2) = subplot(4,1,2); hold on;
for i = 1:numel(rot_trim)
    plot(t_rel, rot_trim{i}, 'Color', [0.7 0.7 1], 'LineWidth', 1);
end
xline(onset_t_rel,  'k--', 'LineWidth', 1.5);
xline(offset_t_rel, 'k--', 'LineWidth', 1.5);
plot(t_rel, mean_rot, 'Color',[0 0 1], 'LineWidth', 2);     % bold mean
ylabel('Rot vel');
title('Odor pulse: Angular velocity');
grid on;
hold off;

% 3) dF/F
ax(3) = subplot(4,1,3); hold on;
for i = 1:numel(dff_trim)
    plot(t_rel, dff_trim{i}, 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end
xline(onset_t_rel,  'k--', 'LineWidth', 1.5);
xline(offset_t_rel, 'k--', 'LineWidth', 1.5);
plot(t_rel, mean_dff, 'Color',[0 0 0], 'LineWidth', 2);     % bold mean
ylabel('dF/F');
title('Odor pulse: dF/F');
grid on;
hold off;

% 4) Odor trace
ax(4) = subplot(4,1,4); hold on;
for i = 1:numel(odor_trim)
    plot(t_rel, odor_trim{i}, 'Color', [1 0.8 0.8], 'LineWidth', 1);
end
xline(onset_t_rel,  'k--', 'LineWidth', 1.5);
xline(offset_t_rel, 'k--', 'LineWidth', 1.5);
plot(t_rel, mean_odor, 'Color',[0.8 0 0], 'LineWidth', 2);  % bold mean
xlabel('Time (s)')
ylabel('Odor');
title('Odor pulse: odor trace');
grid on;
hold off;

sgtitle('Odor-pulse aligned activity');
linkaxes(ax,'x');

save_plot_with_title_as_filename('odor_pulses', 'fwd_rot_dff_odor', savepath);

%% all 3 (fwd, rot, dff) on one plot
% Create the figure and subplots
figure;

% Subplot 1: dFF
subplot(4, 1, 1);  % First subplot
hold on
plot(time, dff, 'b', 'LineWidth', 1.5);  % Plot dFF in red
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F');
grid on;

% Subplot 2: Forward Velocity
subplot(4, 1, 2);  % Second subplot
plot(time, fwd_vel, 'black', 'LineWidth', 1.5);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([-5 20]);

% Subplot 3: Rotational Velocity
subplot(4, 1, 3);  % Third subplot
plot(time, angular_vel, 'r', 'LineWidth', 1.5);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on;
ylim([-5 5]);

% Subplot 4: odor
subplot(4, 1, 4);  % Third subplot
plot(time, odor, 'LineWidth', 1.5);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('odor pulses');
title('Odor Pulses');
grid on;
ylim([-6 1]);

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);


end