function eight_visual_pattern_analysis(a2p_data, savepath, box, split)
%% calculate normalized dff if exists 2 rois
a=size(a2p_data.roi.ts);
if a(1)>1
    dff = a2p_data.roi.ts(2,:)-a2p_data.roi.ts(1,:);
else
    dff = a2p_data.roi.ts(1,:);
end
%% determine the kinematic variables
time = a2p_data.dq(1).t;
time_supp = a2p_data.dq(2).t;
if isfield(a2p_data.dq, 'vho')
    visual_heading = a2p_data.dq(1).vho;
else
    visual_heading = a2p_data.dq(1).vh;
end
if isfield(a2p_data.dq, 'vvyo')
    visual_pattern_velocity = a2p_data.dq(1).vvyo;
else
    visual_pattern_velocity = a2p_data.dq(1).vvy;
end
ball_forward_velocity = a2p_data.dq(1).bvf;
ball_forward_velocity_supp = a2p_data.dq(2).bvf;
ball_yaw_velocity = a2p_data.dq(1).bvy_deg;
ball_yaw_velocity_supp = a2p_data.dq(2).bvy_deg;
ball_side_velocity = a2p_data.dq(1).bvs;
ball_side_velocity_supp = a2p_data.dq(2).bvs;
totalspeed_supp = a2p_data.dq(2).totalspeed;

%% --- LABEL SWITCH ---
if split
    label_cw  = 'Back-to-Front';
    label_ccw = 'Front-to-Back';
else
    % changed bc bright bar was opposite, may need to change when
    % exmamining dark bar
    label_cw  = 'Counterclockwise';
    label_ccw = 'Clockwise';
end

%% Plot the dff and pattern
% Plot the dff and fwd
figure;
yyaxis left;
plot(time, visual_heading, '-k');
ylabel('pattern');  % Label for the left y-axis


% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-b');  % Plot second time series in red
ylabel('dff');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('Pattern and DFF');

% Optional: Add legend
legend('Pattern', 'dff');
%save_plot_with_title_as_filename('Pattern', 'dff', savepath);
%% variables
% frequency
freq = time(end)/length(time);
% Define the number of timepoints to include before and after motion
len_sec = 2; % length in seconds of time to buffer with
stimlen = 2; % how long the pattern moves for
pre_motion_points = round(len_sec/freq);
post_motion_points = round(len_sec/freq);
remove_if_shorter_than = round(1.5/freq);


%value for total speed (moving or no)
still_cutoff = 1;

%% now need to break up trial
heading  = a2p_data.dq(1).vh;
t  = a2p_data.dq(1).t;

edges = [-pi, -pi/2, 0, pi/2, pi];
% Bin numbers:         1      2     3         4

% Find which bin each heading belongs to
which_bin = discretize(heading, edges);

% Now you can extract the indices for each region:
idx_1 = which_bin == 1;          % 0 <= heading < pi/2
idx_2 = which_bin == 2;          % pi/2 <= heading < pi
idx_3 = which_bin == 3;          % pi <= heading < 3*pi/2
idx_4 = which_bin == 4;          % 3*pi/2 <= heading < 2*pi


%% break up into clockwise ccw

% Sub-indices
idx_1_up = idx_1 & (visual_pattern_velocity > 0.1);   % 0-pi/2, positive slope
idx_1_down = idx_1 & (visual_pattern_velocity < -0.1); % 0-pi/2, negative slope

idx_2_up = idx_2 & (visual_pattern_velocity > 0.1);
idx_2_down = idx_2 & (visual_pattern_velocity < -0.10);

idx_3_up = idx_3 & (visual_pattern_velocity > 0.1);
idx_3_down = idx_3 & (visual_pattern_velocity < -0.1);

idx_4_up = idx_4 & (visual_pattern_velocity > 0.1);
idx_4_down = idx_4 & (visual_pattern_velocity < -0.1);

idx_1_up = keep_runs_min_length(idx_1_up, remove_if_shorter_than);
idx_1_down = keep_runs_min_length(idx_1_down, remove_if_shorter_than);
idx_2_down = keep_runs_min_length(idx_2_down, remove_if_shorter_than);
idx_2_up = keep_runs_min_length(idx_2_up, remove_if_shorter_than);
idx_3_up = keep_runs_min_length(idx_3_up, remove_if_shorter_than);
idx_3_down = keep_runs_min_length(idx_3_down, remove_if_shorter_than);
idx_4_up = keep_runs_min_length(idx_4_up, remove_if_shorter_than);
idx_4_down = keep_runs_min_length(idx_4_down, remove_if_shorter_than);

% code
% -pi to -pi/2 = split, -pi/2 to 0 is vgrating, 0-pi/2 is left eye, pi/2 to
% pi is right eye
% up is cw/f2b, down is ccw/b2f

figure; hold on;
plot(t(idx_1_up),   heading(idx_1_up),   'r.', 'DisplayName', 'F2B');
plot(t(idx_1_down), heading(idx_1_down), 'ro', 'DisplayName', 'B2F');
plot(t(idx_2_up),   heading(idx_2_up),   'g.', 'DisplayName', 'CW');
plot(t(idx_2_down), heading(idx_2_down), 'go', 'DisplayName', 'CCW');
plot(t(idx_3_up),   heading(idx_3_up),   'b.', 'DisplayName', 'Left eye F2B');
plot(t(idx_3_down), heading(idx_3_down), 'bo', 'DisplayName', 'Left eye B2F');
plot(t(idx_4_up),   heading(idx_4_up),   'm.', 'DisplayName', 'Right eye F2B');
plot(t(idx_4_down), heading(idx_4_down), 'mo', 'DisplayName', 'Right eye B2F');
legend; hold off;

totalspeed = a2p_data.dq(1).totalspeed;

all_idx = {idx_1_up, idx_1_down, idx_2_up, idx_2_down, ...
           idx_3_up, idx_3_down, idx_4_up, idx_4_down};
all_labels = {'F2B', 'B2F', 'CW', 'CCW', ...
              'Left eye F2B', 'Left eye B2F', 'Right eye F2B', 'Right eye B2F'};

mean_bout_speed = cell(size(all_idx)); % one cell per presentation type

for idx = 1:numel(all_idx)
    mask = all_idx{idx};
    z = [0; mask(:); 0];
    d = diff(z);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    mean_this = nan(1, length(run_starts));
    for i = 1:length(run_starts)
        bout_idx = run_starts(i) : run_ends(i);
        mean_this(i) = mean(totalspeed(bout_idx));
    end
    mean_bout_speed{idx} = mean_this;
    % mean_this is now a vector, each element is the average totalspeed for one bout of this condition
end

% Example: For idx_1_up (first condition), see mean_bout_speed{1}
% disp('Mean bout speed for F2B:')
% disp(mean_bout_speed{1})

%% loop for all
% Assume you have:
% - t: time vector
% - dff: dff vector (same length as t)
% - idx_1_up, idx_1_down, ..., idx_4_down: 8 logical masks (same length)
speed_cutoff = 3;   % Define cutoff
all_idx = {idx_1_up, idx_1_down, idx_2_up, idx_2_down, ...
           idx_3_up, idx_3_down, idx_4_up, idx_4_down};
all_labels = {'F2B', 'B2F', 'CW', 'CCW', ...
              'Left eye F2B', 'Left eye B2F', 'Right eye F2B', 'Right eye B2F'};

win_sec = 2;    % time before/after event (in sec)
dt = median(diff(t));
win_pts = round(win_sec / dt);

figure('Name','DFF Segments by Heading Quadrant & Slope','Position',[100 200 1200 800]);
n_panels = 8;



for ax_idx = 1:n_panels
    mask = all_idx{ax_idx};
    z = [0; mask(:); 0];
    d = diff(z);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    
    % Get mean speeds for this group (computed as shown previously)
    mean_this = mean_bout_speed{ax_idx};

    dff_segs = {};
    seg_colors = {};
    
    for i = 1:length(run_starts)
        center_idx = round((run_starts(i) + run_ends(i))/2);
        seg_idx = (center_idx-win_pts):(center_idx+win_pts);
        if seg_idx(1) < 1 || seg_idx(end) > length(dff)
            continue
        end
        dff_segs{end+1} = dff(seg_idx);
        if mean_this(i) > speed_cutoff
            seg_colors{end+1} = [0 0.7 0];   % green for fast
        else
            seg_colors{end+1} = [0.7 0 0];   % red for slow
        end
    end

    if isempty(dff_segs)
        mean_dff = [];
        t_rel = linspace(-win_sec, win_sec, 2*win_pts+1);
    else
        t_rel = t(1:length(dff_segs{1})) - t(win_pts+1);
        dff_mat = cell2mat(dff_segs');
        mean_dff = mean(dff_mat, 1, 'omitnan');
    end

    subplot(4,2,ax_idx); hold on;
    % Plot each segment with its assigned color
    for i = 1:length(dff_segs)
        plot(t_rel, dff_segs{i}, 'Color', [seg_colors{i} 0.5], 'LineWidth', 1.2);
    end
    % Add mean line from -1 to 1s (all trials pooled)
    if ~isempty(mean_dff)
        idx_line = (t_rel >= -1) & (t_rel <= 1);
        plot(t_rel(idx_line), mean_dff(idx_line), '-', 'Color', [0.2 0.3 0.9], 'LineWidth', 3);
    end
    xline(1, 'k--', 'LineWidth', 1.1);
    xline(-1,   'k--', 'LineWidth', 1.1);
    if ax_idx==1
        % Add dummy plots for legend
        h_fast = plot(nan, nan, '-', 'Color', [0 0.7 0 0.7], 'LineWidth', 2);
        h_slow = plot(nan, nan, '-', 'Color', [0.7 0 0 0.7], 'LineWidth', 2);
        h_mean = plot(nan, nan, '-', 'Color', [0.2 0.3 0.9], 'LineWidth', 3);
        legend([h_fast h_slow h_mean], ...
            {sprintf('Fast (>%.2f)', speed_cutoff),...
             sprintf('Slow (<%.2f)', speed_cutoff),...
             'Mean DFF'},...
             'Location','best')
    end

    xlabel('Time (s)');
    ylabel('DFF');
    title(all_labels{ax_idx});
    hold off;
end

sgtitle('DFF Based on Pattern Motion & Bout Avg Speed');
save_plot_with_title_as_filename('all8', 'optic', savepath)
%% now plotting all at once

% Assumes you have 'time', 'dff', 'ball_forward_velocity_supp', 'ball_yaw_velocity_supp' (vectors)
% and idx_1_up ... idx_4_down (logical masks)
idxs = {idx_1_up, idx_1_down, idx_2_up, idx_2_down, idx_3_up, idx_3_down, idx_4_up, idx_4_down};
colors = [ 0.8 0 0; 0.8 0.6 0; 0 0.8 0; 0 0.6 0.8; 0 0 0.8; 0.6 0 0.8; 0.2 0.2 0.2; 0.8 0.4 0.8];
labels = {'F2B', 'B2F', 'CW', 'CCW', ...
              'Left eye F2B', 'Left eye B2F', 'Right eye F2B', 'Right eye B2F'};

vals = {dff, ball_forward_velocity, ball_yaw_velocity};
val_labels = {'dF/F','Forward Velocity (mm/s)','Rotational Velocity (deg/s)'};
subplot_titles = {'dF/F','Forward Velocity','Rotational Velocity'};

figure('Position',[200 50 1400 950]);
for i = 1:3
    ax = subplot(3,1,i); hold on;
    h_signal = plot(time, vals{i}, 'k', 'LineWidth', 1.3, 'DisplayName','Signal');
    h_patch = gobjects(1,8); % For legend
    for j = 1:8
        mask = idxs{j};
        mask_padded = [0; mask(:); 0];
        dmask = diff(mask_padded);
        seg_start = find(dmask==1);
        seg_end   = find(dmask==-1)-1;
        for k = 1:length(seg_start)
            inds = seg_start(k):seg_end(k);
            xfill = [time(inds) fliplr(time(inds))];
            ybase = min(vals{i});
            yfill = [ones(1,length(inds))*ybase, vals{i}(inds(end:-1:1))];
            if k==1, % Only add to legend for first segment of this type
                h_patch(j) = fill(xfill, yfill, colors(j,:), ...
                    'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName',labels{j});
            else
                fill(xfill, yfill, colors(j,:), ...
                    'FaceAlpha',0.3, 'EdgeColor','none', 'HandleVisibility','off');
            end
        end
    end
    plot(time, vals{i}, 'k', 'LineWidth', 1.3, 'HandleVisibility','off');
    ylim auto; grid on
    ylabel(val_labels{i});
    title(subplot_titles{i});
    if i == 1
        legend([h_signal h_patch], [{'Signal'}, labels], 'Location', 'northeastoutside');
    end
    if i==3
        xlabel('Time (s)');
    end
end
sgtitle('dF/F, Forward, and Rotational Velocity with Labeled Heading+Slope Regions');
ax = findall(gcf, 'Type', 'axes');
linkaxes(ax, 'x');


save_plot_with_title_as_filename('multi_heading', 'patch_optic', savepath)


end