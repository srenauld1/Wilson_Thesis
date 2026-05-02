
%% function
function a2p_data = plotting_ves041(a2p_data, jump, savepath, odor)

close all
%% define variables
time = a2p_data.dq(1).t;                  % Time points
dff = a2p_data.roi.ts(1, :);                    % Example dFF data
time_kin = a2p_data.dq(2).t;        
fwd = a2p_data.dq(1).bvf; % Example forward velocity data
fwd_supp = a2p_data.dq(2).bvf;
rot = a2p_data.dq(1).bvy_deg; % Example rotational velocity data
rot_supp = a2p_data.dq(2).bvy_deg;
side = a2p_data.dq(1).bvs_deg;
sample_rate = 1/(max(time)/length(time)); %samples per second

if odor
    rot = a2p_data.dq(1).vvy; % Example rotational velocity data
    rot_supp = a2p_data.dq(2).vvy;
end

%% all 3 (fwd, rot, dff) on one plot
% Create the figure and subplots
figure;

% Subplot 1: dFF
subplot(3, 1, 1);  % First subplot
hold on
plot(time, dff, 'b', 'LineWidth', 1.5);  % Plot dFF in red
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F');
grid on;

% Subplot 2: Forward Velocity
subplot(3, 1, 2);  % Second subplot
plot(time_kin, fwd_supp, 'black', 'LineWidth', 1.5);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([-5 20]);

% Subplot 3: Rotational Velocity
subplot(3, 1, 3);  % Third subplot
plot(time_kin, rot_supp, 'r', 'LineWidth', 1.5);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on;
%ylim([-500 500]);

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);


%% now interactively plot dff and yaw velocity
%trajectory_region_selector(a2p_data, dff, jump)

%% aligned to motion on plotting dff, forward vel, rot vel (all aligned)
preN = round(sample_rate);
offN = 0;
min_run_len = round(sample_rate);
motion = a2p_data.dq(1).motion.moving_or_not(:);
fluo = dff';
rot_speed = abs(rot); % rotational velocity in degrees/sec
t = a2p_data.dq(1).t(:); % your time vector

d = diff([0; motion; 0]);
onsets = find(d == 1);
offsets = find(d == -1) - 1;

move_lengths = offsets - onsets + 1;
move_keep = move_lengths >= min_run_len;
onsets = onsets(move_keep);
offsets = offsets(move_keep);

prev_offsets = [0; offsets(1:end-1)];
still_lengths = onsets - prev_offsets - 1;
still_keep = still_lengths >= min_run_len;

onsets = onsets(still_keep);
offsets = offsets(still_keep);

fluo_cell = {}; forvel_cell = {}; rotsel_cell = {}; t_cell = {};
for i=1:numel(onsets)
    t_start = onsets(i) - preN;
    t_end   = offsets(i) + offN;
    if t_start < 1 || t_end > numel(fluo), continue; end
    idx_epoch = t_start:t_end;
    t_epoch = t(idx_epoch);
    t_onset = t(onsets(i));
    t_rel = t_epoch - t_onset;
    fluo_cell{end+1}    = fluo(idx_epoch);
    forvel_cell{end+1}  = fwd(idx_epoch);
    rotsel_cell{end+1}  = rot_speed(idx_epoch);
    t_cell{end+1}       = t_rel;
end

% Robust dt and restrict to epochs with enough pre samples
dt_median = median(cellfun(@(t) median(diff(t)), t_cell));
t_pre = -preN*dt_median;
bouts_with_all_pre = cellfun(@(t) any(t <= t_pre), t_cell);
fluo_cell    = fluo_cell(bouts_with_all_pre);
forvel_cell  = forvel_cell(bouts_with_all_pre);
rotsel_cell  = rotsel_cell(bouts_with_all_pre);
t_cell       = t_cell(bouts_with_all_pre);

% Mean axis: overlapping duration only
t_post_max = min(cellfun(@(t) max(t), t_cell));
t_mean = t_pre:dt_median:t_post_max;

% Interpolate all epochs to t_mean for averaging
fluo_mat    = nan(length(t_mean), numel(fluo_cell));
forvel_mat  = nan(length(t_mean), numel(forvel_cell));
rotsel_mat  = nan(length(t_mean), numel(rotsel_cell));
for i = 1:numel(fluo_cell)
    fluo_mat(:,i)   = interp1(t_cell{i}, fluo_cell{i}, t_mean, 'linear', nan);
    forvel_mat(:,i) = interp1(t_cell{i}, forvel_cell{i}, t_mean, 'linear', nan);
    rotsel_mat(:,i) = interp1(t_cell{i}, rotsel_cell{i}, t_mean, 'linear', nan);
end

% ---- Plot all runs (variable) and mean (fixed axis) ----
figure;

% DFF
subplot(3,1,1); hold on;
for i = 1:numel(fluo_cell)
    plot(t_cell{i}, fluo_cell{i}, 'Color', [0.2549 0.4118 0.8824], 'HandleVisibility','off'); % royal blue
end
plot(t_mean, nanmean(fluo_mat,2), 'LineWidth', 2, 'Color', [0 0 0.5451], 'DisplayName', 'Mean'); % dark blue
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement onset (s)');
ylabel('DF/F');
title('dF/F aligned to movement onset');
hold off;

% Forward velocity
subplot(3,1,2); hold on;
for i = 1:numel(forvel_cell)
    plot(t_cell{i}, forvel_cell{i}, 'Color', [0.5 0.5 0.5], 'HandleVisibility','off'); % grey
end
plot(t_mean, nanmean(forvel_mat,2), 'Color', [0 0 0], 'LineWidth', 2, 'DisplayName', 'Mean'); % black
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement onset (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity aligned to movement onset');
hold off;

% Rotational velocity (deg/sec)
subplot(3,1,3); hold on;
for i = 1:numel(rotsel_cell)
    plot(t_cell{i}, rotsel_cell{i}, 'Color', [0.7 0 0], 'HandleVisibility','off'); % dark red
end
plot(t_mean, nanmean(rotsel_mat,2), 'Color', [0.3 0 0], 'LineWidth', 2, 'DisplayName', 'Mean'); % darker red
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement onset (s)');
ylabel('Rotational Speed (deg/s)');
title('Rotational speed aligned to movement onset');
hold off;

n_epochs = numel(fluo_cell);
subplot(3,1,1); % annotate top plot
yl = ylim; xl = xlim;
text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'EdgeColor','k');

save_plot_with_title_as_filename('movement_onset', 'variables', savepath);


%% aligned to motion OFFSET, plotting dff and fwd
preN = 1*round(sample_rate);    % how many samples to show before motion offset
min_run_len = round(sample_rate); % 1 second off
offN = 0;     % set >0 if you want more after next offset
rot_speed  = abs(rot);    % <--- rotational velocity in deg/s

% Find all motion on (onsets) and off (offsets)
d = diff([0; motion; 0]);
onsets  = find(d==1);
offsets = find(d==-1)-1;

% Offset refers to when motion turns off (so these are your alignment points)
% Only keep off periods longer than min_run_len
my_onsets = nan(size(offsets)); % next onset after each offset
for i = 1:length(offsets)
    next_on = find(onsets > offsets(i),1,'first');
    if ~isempty(next_on)
        my_onsets(i) = onsets(next_on);
    else
        my_onsets(i) = numel(motion)+1; % i.e. end of data
    end
end

off_lens = my_onsets - offsets;
keep = off_lens >= min_run_len;
offsets = offsets(keep);
my_onsets = my_onsets(keep);

fluo_cell = {}; behav_cell = {}; rot_cell = {}; t_cell = {};

for i = 1:numel(offsets)
    t_start = offsets(i) - preN + 1; % +1 so offset itself is always at same place in vector
    if t_start < 1, continue; end  % skip if not enough pre points
    t_end = my_onsets(i) - 1 + offN; % epoch ends at just before the next movement bout, or can add post points
    t_end = min(t_end,numel(fluo));  % don't run off end
    idx_epoch = t_start:t_end;
    t_offset = t(offsets(i));
    t_rel = t(idx_epoch) - t_offset; % time relative to offset (offset=0)
    fluo_cell{end+1}  = fluo(idx_epoch);
    behav_cell{end+1} = fwd(idx_epoch);
    rot_cell{end+1}   = rot_speed(idx_epoch);
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

% Plot
figure;

subplot(3,1,1); hold on;
for i = 1:numel(fluo_cell)
    plot(t_cell{i}, fluo_cell{i}, 'Color', [0.2549 0.4118 0.8824], 'HandleVisibility','off');
end
plot(t_mean, nanmean(fluo_mat,2), 'Color', [0 0 0.5451], 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement offset (s)');
ylabel('DF/F');
title('dF/F aligned to movement offset');
hold off;

subplot(3,1,2); hold on;
for i = 1:numel(behav_cell)
    plot(t_cell{i}, behav_cell{i}, 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
end
plot(t_mean, nanmean(behav_mat,2), 'k', 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time from movement offset (s)');
ylabel('Forward velocity');
title('Forward velocity aligned to movement offset');
hold off;

subplot(3,1,3); hold on;
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
subplot(3,1,1);
yl = ylim; xl = xlim;
text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'EdgeColor','k');

save_plot_with_title_as_filename('movement_offset', 'variables', savepath);


%% old, more detailed plotting
%% Plot the dff and fwd
% % Plot the dff and fwd
% figure;
% 
% 
% yyaxis left;
% plot(time_kin, fwd, '-b', 'LineWidth', 1.5)% Plot first time series in blue
% ylabel('Forward Velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-5 10]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, dff, '-r' , 'LineWidth', 1.5);  % Plot second time series in red
% ylabel('dff');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('Forward velocity and DFF');
% 
% % Optional: Add legend
% legend('Forward Velocity', 'dff');
% save_plot_with_title_as_filename('Fwd vel', 'dff', savepath);
% 
% %% plot dff and yaw
% figure;
% yyaxis left;
% plot(time_kin, yaw, '-g', 'LineWidth', 1.5);  % Plot first time series in blue
% ylabel('Yaw velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-500 500]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, dff, '-r', 'LineWidth', 1.5);  % Plot second time series in red
% ylabel('dff');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('Yaw velocity and DFF');
% 
% % Optional: Add legend
% legend('yaw velocity', 'dff');
% save_plot_with_title_as_filename('Yaw vel', 'dff', savepath);

%% plot fwd and yaw
% figure;
% yyaxis left;
% plot(time_kin, fwd, '-b');  % Plot first time series in blue
% ylabel('Forward velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-5 5]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time_kin, yaw, '-g');  % Plot second time series in red
% ylabel('yaw vel');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('forward and yaw velocity');
% 
% % Optional: Add legend
% legend('fwd velocity', 'yaw velocity');
% save_plot_with_title_as_filename('Fwd vel', 'Yaw vel', savepath);
% 
% 
% %% Plot the side and dff
% % Plot the dff and fwd
% figure;
% yyaxis left;
% plot(time_kin, side, '-y');  % Plot first time series in blue
% ylabel('Sideways velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-500 500]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, dff, '-r');  % Plot second time series in red
% ylabel('dff');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('Side velocity vs dff');
% 
% % Optional: Add legend
% legend('Sideways velocity', 'dff');
% save_plot_with_title_as_filename('Side vel', 'dff', savepath);


%% Parameters
preN = 2*round(sample_rate);
offN = 2*round(sample_rate);
min_run_len = 0.5*round(sample_rate);
min_inbetween_len = round(sample_rate);

%% Find bouts
d = diff([0; motion; 0]);
onsets = find(d == 1);
offsets = find(d == -1) - 1;

move_lengths = offsets - onsets + 1;
move_keep = move_lengths >= min_run_len;
onsets = onsets(move_keep);
offsets = offsets(move_keep);

% Compute PRE-still (before each bout)
prev_offsets = [0; offsets(1:end-1)];
pre_still = onsets - prev_offsets - 1;
is_strict_pre_still = false(size(onsets));
for i = 1:numel(onsets)
    idx = max(prev_offsets(i)+1,1):onsets(i)-1;  % interval BEFORE each onset
    if pre_still(i) >= min_inbetween_len && all(motion(idx)==0)
        is_strict_pre_still(i) = true;
    end
end

% Compute POST-still (after each bout)
next_onsets = [onsets(2:end); numel(motion)+1];
post_still = next_onsets - offsets - 1;
is_strict_post_still = false(size(offsets));
for i = 1:numel(offsets)
    idx = offsets(i)+1 : min(offsets(i)+min_inbetween_len,numel(motion));
    if post_still(i) >= min_inbetween_len && all(motion(idx)==0)
        is_strict_post_still(i) = true;
    end
end

% --- GROUP MASKS (STRICT) ---
onset_long   = is_strict_pre_still;
onset_short  = ~is_strict_pre_still;
offset_long  = is_strict_post_still;
offset_short = ~is_strict_post_still;

% --- PLOT ONSET-ALIGNED GROUPS ---
onset_groups = {'Long-pre-still onset','Short-pre-still onset'};
onset_masks = {onset_long, onset_short};
onset_colors = {[0.3 0.5 1],[1 0.4 0.4]};

for g = 1:2
    mask = onset_masks{g};
    curr_onsets = onsets(mask);
    curr_offsets = offsets(mask);

    fluo_cell = {}; forvel_cell = {}; rotsel_cell = {}; t_cell = {};
    for i = 1:numel(curr_onsets)
        t_start = curr_onsets(i) - preN;
        t_end   = curr_onsets(i) + offN;
        if t_start < 1 || t_end > numel(fluo), continue; end
        idx_epoch = t_start:t_end;
        t_epoch = t(idx_epoch);
        t_onset = t(curr_onsets(i));
        t_rel = t_epoch - t_onset;
        fluo_cell{end+1}    = fluo(idx_epoch);
        forvel_cell{end+1}  = fwd(idx_epoch);
        rotsel_cell{end+1}  = rot_speed(idx_epoch);
        t_cell{end+1}       = t_rel;
    end

    if isempty(fluo_cell), continue; end

    dt_median = median(cellfun(@(t) median(diff(t)), t_cell));
    t_pre = -preN*dt_median;
    t_post = offN*dt_median;
    t_mean = t_pre:dt_median:t_post;

    fluo_mat    = nan(length(t_mean), numel(fluo_cell));
    forvel_mat  = nan(length(t_mean), numel(forvel_cell));
    rotsel_mat  = nan(length(t_mean), numel(rotsel_cell));
    for i = 1:numel(fluo_cell)
        if numel(t_cell{i})<2, continue, end
        fluo_mat(:,i)   = interp1(t_cell{i}, fluo_cell{i}, t_mean, 'linear', nan);
        forvel_mat(:,i) = interp1(t_cell{i}, forvel_cell{i}, t_mean, 'linear', nan);
        rotsel_mat(:,i) = interp1(t_cell{i}, rotsel_cell{i}, t_mean, 'linear', nan);
    end

    figure('Name',onset_groups{g});
    subplot(3,1,1); hold on;
    for i = 1:numel(fluo_cell)
        plot(t_cell{i}, fluo_cell{i}, 'Color', [onset_colors{g}, 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(fluo_mat,2), 'LineWidth', 2, 'Color', onset_colors{g}, 'DisplayName', 'Mean');
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    xlabel('Time from onset (s)'); ylabel('DF/F');
    title(onset_groups{g}); hold off;

    subplot(3,1,2); hold on;
    for i = 1:numel(forvel_cell)
        plot(t_cell{i}, forvel_cell{i}, 'Color', [0.5 0.5 0.5 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(forvel_mat,2), 'k', 'LineWidth', 2);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    ylabel('Forward Velocity'); title('Forward Velocity');
    hold off;

    subplot(3,1,3); hold on;
    for i = 1:numel(rotsel_cell)
        plot(t_cell{i}, rotsel_cell{i}, 'Color', [0.7 0 0 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(rotsel_mat,2), 'Color', [0.3 0 0], 'LineWidth', 2);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    ylabel('Rotational Speed'); xlabel('Time from onset (s)');
    title('Rotational Speed');
    hold off;

    n_epochs = numel(fluo_cell); subplot(3,1,1);
    yl = ylim; xl = xlim;
    text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor','k');
    sgtitle(sprintf('%s | preN=%d, offN=%d, min_run_len=%d, min_inbetween_len=%d, n=%d', ...
        onset_groups{g}, preN, offN, min_run_len, min_inbetween_len, n_epochs));
    
    % ---- Save the current figure ----
    % The filename can encode the group and params
    fname = sprintf('%s_onset_preN%d_offN%d_minRun%d_minStill%d_n%d', ...
        regexprep(onset_groups{g},'[\W]+','_'), preN, offN, min_run_len, min_inbetween_len, n_epochs);
    save_plot_with_title_as_filename('Onset', fname, savepath);
end

% --- OFFSET-ALIGNED ---
offset_groups = {'Long-post-still offset','Short-post-still offset'};
offset_masks = {offset_long, offset_short};
offset_colors = {[0.2 0.8 0.2],[1 0.6 0.1]};  % green and orange
for g = 1:2
    mask = offset_masks{g};
    curr_onsets = onsets(mask);
    curr_offsets = offsets(mask);

    fluo_cell = {}; forvel_cell = {}; rotsel_cell = {}; t_cell = {};
    for i = 1:numel(curr_offsets)
        t_start = curr_offsets(i) - preN;
        t_end = curr_offsets(i) + offN;
        if t_start < 1 || t_end > numel(fluo), continue; end
        idx_epoch = t_start:t_end;
        t_epoch = t(idx_epoch);
        t_offset = t(curr_offsets(i));
        t_rel = t_epoch - t_offset;
        fluo_cell{end+1}    = fluo(idx_epoch);
        forvel_cell{end+1}  = fwd(idx_epoch);
        rotsel_cell{end+1}  = rot_speed(idx_epoch);
        t_cell{end+1}       = t_rel;
    end

    if isempty(fluo_cell), continue; end

    dt_median = median(cellfun(@(t) median(diff(t)), t_cell));
    t_pre = -preN*dt_median;
    t_post = offN*dt_median;
    t_mean = t_pre:dt_median:t_post;

    fluo_mat    = nan(length(t_mean), numel(fluo_cell));
    forvel_mat  = nan(length(t_mean), numel(forvel_cell));
    rotsel_mat  = nan(length(t_mean), numel(rotsel_cell));
    for i = 1:numel(fluo_cell)
        if numel(t_cell{i})<2, continue, end
        fluo_mat(:,i)   = interp1(t_cell{i}, fluo_cell{i}, t_mean, 'linear', nan);
        forvel_mat(:,i) = interp1(t_cell{i}, forvel_cell{i}, t_mean, 'linear', nan);
        rotsel_mat(:,i) = interp1(t_cell{i}, rotsel_cell{i}, t_mean, 'linear', nan);
    end

    figure('Name',offset_groups{g});
    subplot(3,1,1); hold on;
    for i = 1:numel(fluo_cell)
        plot(t_cell{i}, fluo_cell{i}, 'Color', [offset_colors{g}, 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(fluo_mat,2), 'LineWidth', 2, 'Color', offset_colors{g}, 'DisplayName', 'Mean');
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    xlabel('Time from offset (s)'); ylabel('DF/F');
    title(offset_groups{g}); hold off;

    subplot(3,1,2); hold on;
    for i = 1:numel(forvel_cell)
        plot(t_cell{i}, forvel_cell{i}, 'Color', [0.5 0.5 0.5 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(forvel_mat,2), 'k', 'LineWidth', 2);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    ylabel('Forward Velocity'); title('Forward Velocity');
    hold off;

    subplot(3,1,3); hold on;
    for i = 1:numel(rotsel_cell)
        plot(t_cell{i}, rotsel_cell{i}, 'Color', [0.7 0 0 0.2], 'HandleVisibility','off');
    end
    plot(t_mean, nanmean(rotsel_mat,2), 'Color', [0.3 0 0], 'LineWidth', 2);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    ylabel('Rotational Speed'); xlabel('Time from offset (s)');
    title('Rotational Speed');
    hold off;

    n_epochs = numel(fluo_cell); subplot(3,1,1);
    yl = ylim; xl = xlim;
    text(xl(2), yl(2), sprintf('n epochs = %d', n_epochs), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor','k');
    sgtitle(sprintf('%s | preN=%d, offN=%d, min_run_len=%d, min_inbetween_len=%d, n=%d', ...
        offset_groups{g}, preN, offN, min_run_len, min_inbetween_len, n_epochs));

    % ---- Save the current figure ----
    % The filename can encode the group and params
    fname = sprintf('%s_offset_preN%d_offN%d_minRun%d_minStill%d_n%d', ...
        regexprep(offset_groups{g},'[\W]+','_'), preN, offN, min_run_len, min_inbetween_len, n_epochs);
    save_plot_with_title_as_filename('Offset', fname, savepath);
end

%% dff peaks
dff_smooth2 = smooth(smooth(dff, 'rlowess'), 50, 'rlowess');
%findpeaks(dff_smooth2, 'MinPeakDistance', 25, 'MinPeakHeight', 2);
figure
findpeaks(dff, 'MinPeakDistance', 30, 'MinPeakHeight', 3);
[pks, locs] = findpeaks(dff, 'MinPeakDistance', 30, 'MinPeakHeight', 2.5);

window_s = 2; % seconds before/after
window_N = round(window_s * sample_rate);

% Ensure column vectors
dff = dff(:);
bvf = fwd(:);
rots = abs(rot(:)); % rotational speed (deg/s)
t = t(:);

dff_epochs = {};
bvf_epochs = {};
rots_epochs = {};
t_epochs   = {};

for i = 1:numel(locs)
    idx_c = locs(i);
    idx_start = idx_c - window_N;
    idx_end   = idx_c + window_N;
    if idx_start < 1 || idx_end > length(dff), continue, end
    idx_win = idx_start:idx_end;
    t0 = t(idx_c);
    t_epoch = t(idx_win) - t0;
    dff_epochs{end+1}  = dff(idx_win);
    bvf_epochs{end+1}  = bvf(idx_win);
    rots_epochs{end+1} = rots(idx_win);
    t_epochs{end+1}    = t_epoch;
end

n_epochs = numel(bvf_epochs);
if n_epochs < 1
    error('No useable peaks found with sufficient window length!');
end

dt_median = median(cellfun(@(t) median(diff(t)), t_epochs));
t_mean = -window_N*dt_median : dt_median : window_N*dt_median;

dff_mat  = nan(length(t_mean), n_epochs);
bvf_mat  = nan(length(t_mean), n_epochs);
rots_mat = nan(length(t_mean), n_epochs);

for i = 1:n_epochs
    dff_mat(:,i)  = interp1(t_epochs{i}, dff_epochs{i}, t_mean, 'linear', nan);
    bvf_mat(:,i)  = interp1(t_epochs{i}, bvf_epochs{i}, t_mean, 'linear', nan);
    rots_mat(:,i) = interp1(t_epochs{i}, rots_epochs{i}, t_mean, 'linear', nan);
end

figure;

% dF/F
subplot(3,1,1); hold on;
for i = 1:n_epochs
    plot(t_epochs{i}, dff_epochs{i}, 'Color', [0.4 0.8 1 0.2],'HandleVisibility','off');
end
plot(t_mean, nanmean(dff_mat,2), 'b', 'LineWidth', 2);
line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth',2);
xlabel('Time from dF/F peak (s)');
ylabel('dF/F');
title('dF/F aligned to dF/F peak');
hold off;

% Forward velocity
subplot(3,1,2); hold on;
for i = 1:n_epochs
    plot(t_epochs{i}, bvf_epochs{i}, 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
end
plot(t_mean, nanmean(bvf_mat,2), 'k', 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth',2);
xlabel('Time from dF/F peak (s)');
ylabel('Forward velocity');
title(sprintf('Forward velocity aligned to dF/F peak (n peaks=%d)', n_epochs));
hold off;

% Rotational velocity
subplot(3,1,3); hold on;
for i = 1:n_epochs
    plot(t_epochs{i}, rots_epochs{i}, 'Color', [1 0.5 0.5 0.2], 'HandleVisibility','off');
end
plot(t_mean, nanmean(rots_mat,2), 'r', 'LineWidth', 2, 'DisplayName', 'Mean');
line([0 0], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth',2);
xlabel('Time from dF/F peak (s)');
ylabel('Rotational velocity (deg/s)');
title('Rotational velocity aligned to dF/F peak');
hold off;

sgtitle('dF/F and behavior aligned to dF/F peaks');
save_plot_with_title_as_filename('dff_peaks', 'kinematics', savepath);


%% now specifically looking at menotactic periods + 2 sec on either side (??)
is_menotaxing = a2p_data.menotaxis.is_menotaxing;


time = time(:)';
is_menotaxing_on_dff = interp1(time_kin(:), is_menotaxing(:), time, 'nearest', 0); % 0 for extrapolation
is_menotaxing_on_dff = is_menotaxing_on_dff > 0; % logical

figure;

% --- Subplot 1: dF/F ---
subplot(3,1,1); hold on
% Plot light yellow boxes where is_menotaxing_on_dff == 1 SEGMENT FWD VEL
% TOGGLE
meno = is_menotaxing_on_dff(:);
t = time(:);
signal = dff(:);

onst = [false; diff(meno)==1];
offst = [false; diff(meno)==-1];
starts = find(onst);
ends   = find(offst);

if meno(1)>0, starts = [1; starts]; end
if meno(end)>0, ends = [ends; length(meno)]; end

yl = [min(signal) max(signal)];
for i = 1:length(starts)
    patch([t(starts(i)) t(ends(i)) t(ends(i)) t(starts(i))], ...
          [yl(1) yl(1) yl(2) yl(2)], [1 .95 .6], ...
          'EdgeColor','none','FaceAlpha',0.35);
end

plot(t, signal, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F (yellow = menotaxing)');
grid on;
hold off

% --- Subplot 2: Forward Velocity ---
subplot(3,1,2); hold on
% Upsample/interpolate meno to time_kin if needed
meno_kin = interp1(t, double(meno), time_kin(:), 'nearest', 0);
signal = fwd_supp(:);
tk = time_kin(:);

yl = [min(signal) max(signal)];
for i = 1:length(starts)
    % Find interpolated window in time_kin
    t1 = t(starts(i)); t2 = t(ends(i));
    i1 = find(tk >= t1, 1, 'first');
    i2 = find(tk <= t2, 1, 'last');
    if ~isempty(i1) && ~isempty(i2) && i2>i1
        patch([tk(i1) tk(i2) tk(i2) tk(i1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 .95 .6], ...
              'EdgeColor','none','FaceAlpha',0.35);
    end
end
plot(tk, signal, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([-5 20]);
hold off

% --- Subplot 3: Rotational Velocity ---
subplot(3,1,3); hold on
signal = rot_supp(:);
tk = time_kin(:);

yl = [min(signal) max(signal)];
for i = 1:length(starts)
    % Find interpolated window in time_kin
    t1 = t(starts(i)); t2 = t(ends(i));
    i1 = find(tk >= t1, 1, 'first');
    i2 = find(tk <= t2, 1, 'last');
    if ~isempty(i1) && ~isempty(i2) && i2>i1
        patch([tk(i1) tk(i2) tk(i2) tk(i1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 .95 .6], ...
              'EdgeColor','none','FaceAlpha',0.35);
    end
end
plot(tk, signal, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on;
%ylim([-500 500]);
hold off

% --- Link the x-axes of all subplots ---
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% --- Add overall supertitle ---
sgtitle('dF/F, Forward, and Rotational Velocity (yellow = menotaxing)');

save_plot_with_title_as_filename('fwd_rot', 'dff_separate_menotaxis', savepath);

end

