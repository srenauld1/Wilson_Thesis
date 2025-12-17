function a2p_data = process_pause_data(a2p_data, savepath, box, split)

%% calculate normalized dff if exists 2 rois
dff = a2p_data.roi.ts(1,:);
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
% frequency
freq = time(end)/length(time);
remove_if_shorter_than = round(8/freq);
remove_if_longer_than = round(12/freq);
remove_if_shorter_than_ol = round(1.8/freq);
remove_if_longer_than_ol = round(4.5/freq);

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
save_plot_with_title_as_filename('Pattern', 'dff', savepath);

%% segment out pauses
pauses = zeros(size(visual_heading));
pauses(abs(visual_pattern_velocity)<= 0.03) = 1;

diff_array = diff([0, pauses, 0]);
run_starts = find(diff_array == 1);
run_ends = find(diff_array == -1) - 1;

for i = 1:length(run_starts)
    run_length = run_ends(i) - run_starts(i) + 1;
    if run_length < remove_if_shorter_than
        pauses(run_starts(i):run_ends(i)) = 0;
    end
    if run_length > remove_if_longer_than
        pauses(run_starts(i):run_ends(i)) = 0;
    end
end

figure
hold on
plot(a2p_data.dq(1).t, pauses, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Pauses','Visual Pattern Velocity');
title('Pauses (red) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');

disp('Zoom and pan as needed. Drag a rectangle for each pause interval.');
disp('When done, *click* "Select" on the Figure toolbar, then press Enter in the Command Window.');

while true
    h = drawrectangle('Label','Pause','Color','m'); 
    % Wait for user input (double-click or Enter finishes interval)
    wait(h);
    rectpos = h.Position;
    % If rectangle has zero width (user pressed enter without drawing), break
    if rectpos(3) == 0
        delete(h);
        break;
    end
    t_start = rectpos(1);
    t_end = t_start + rectpos(3);
    indices = find(a2p_data.dq(1).t >= t_start & a2p_data.dq(1).t <= t_end);
    pauses(indices) = 1;
    % Optionally, show points
    plot(a2p_data.dq(1).t(indices), visual_pattern_velocity(indices), 'ko','MarkerFaceColor','y','DisplayName','Hand-picked pause');
    delete(h); % remove the drawn rectangle
end

% Show the final result
figure; hold on
plot(a2p_data.dq(1).t, pauses, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Pauses','Visual Pattern Velocity');
title('Pauses (red, incl. hand-picked) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');


%% segment out open loop rotations

% front to back
open_loop_f2b = zeros(size(visual_heading));
open_loop_f2b(abs(visual_pattern_velocity)>=1) = 1;
% front to back
open_loop_b2f = zeros(size(visual_heading));
open_loop_b2f(abs(visual_pattern_velocity)<=-1) = 1;

diff_array = diff([0, open_loop_f2b, 0]);
run_starts = find(diff_array == 1);
run_ends = find(diff_array == -1) - 1;

for i = 1:length(run_starts)
    run_length = run_ends(i) - run_starts(i) + 1;
    if run_length < remove_if_shorter_than_ol
        open_loop_f2b(run_starts(i):run_ends(i)) = 0;
    end
    if run_length > remove_if_longer_than_ol
        open_loop_f2b(run_starts(i):run_ends(i)) = 0;
    end
end

% plot pauses
% figure
% hold on
% plot(pauses)
% plot(open_loop)
% plot(visual_pattern_velocity)

%% segment out backwards gain periods (by hand)

% Initialize reverse_gain to 0, or load if it already exists
if ~isfield(a2p_data, 'reverse_gain') || isempty(a2p_data.reverse_gain)
    reverse_gain = zeros(size(a2p_data.dq(1).t));
else
    reverse_gain = a2p_data.reverse_gain;
end

figure; hold on
plot(a2p_data.dq(1).t, reverse_gain, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Reverse Gain','Visual Pattern Velocity');
title('Reverse Gain (red) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');

disp('Zoom/pan as needed, then click and drag to draw a rectangle for each reverse-gain interval.');
disp('When you are finished, just single-click on the plot and press Enter without drawing a rectangle.');

while true
    h = drawrectangle('Label','Reverse Gain','Color','m');
    wait(h); % Wait for rectangle to be created/finished
    rectpos = h.Position;
    % If rectangle has zero width, break out
    if rectpos(3) == 0
        delete(h);
        break;
    end
    t_start = rectpos(1);
    t_end = t_start + rectpos(3);
    indices = find(a2p_data.dq(1).t >= t_start & a2p_data.dq(1).t <= t_end);
    reverse_gain(indices) = 1;
    % Visually mark selected interval
    plot(a2p_data.dq(1).t(indices), visual_pattern_velocity(indices),...
        'ko', 'MarkerFaceColor','y','DisplayName','Hand-picked Reverse Gain');
    delete(h);
end

% Assume pauses, open_loop, and reverse_gain are all vectors the same size as a2p_data.dq(1).t

bad_indices = (pauses == 1) | (open_loop_f2b == 1);
reverse_gain(bad_indices) = 0;

% Optionally, update struct if storing
a2p_data.reverse_gain = reverse_gain;

% Show the final result for review
figure; hold on
plot(a2p_data.dq(1).t, reverse_gain, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Reverse Gain','Visual Pattern Velocity');
title('Reverse Gain (red, incl. hand-picked) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');


%% segment out replay (by hand)

% Initialize reverse_gain to 0, or load if it already exists
if ~isfield(a2p_data, 'replay') || isempty(a2p_data.replay)
    replay = zeros(size(a2p_data.dq(1).t));
else
    replay = a2p_data.replay;
end

figure; hold on
plot(a2p_data.dq(1).t, replay, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Replay','Visual Pattern Velocity');
title('Replay (red) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');

disp('Zoom/pan as needed, then click and drag to draw a rectangle for each reverse-gain interval.');
disp('When you are finished, just single-click on the plot and press Enter without drawing a rectangle.');

while true
    h = drawrectangle('Label','Replay','Color','m');
    wait(h); % Wait for rectangle to be created/finished
    rectpos = h.Position;
    % If rectangle has zero width, break out
    if rectpos(3) == 0
        delete(h);
        break;
    end
    t_start = rectpos(1);
    t_end = t_start + rectpos(3);
    indices = find(a2p_data.dq(1).t >= t_start & a2p_data.dq(1).t <= t_end);
    replay(indices) = 1;
    % Visually mark selected interval
    plot(a2p_data.dq(1).t(indices), visual_pattern_velocity(indices),...
        'ko', 'MarkerFaceColor','y','DisplayName','Hand-picked Reverse Gain');
    delete(h);
end

% Assume pauses, open_loop, and reverse_gain are all vectors the same size as a2p_data.dq(1).t

bad_indices = (pauses == 1) | (open_loop_f2b == 1);
replay(bad_indices) = 0;

% Optionally, update struct if storing
a2p_data.replay = replay;

% Show the final result for review
figure; hold on
plot(a2p_data.dq(1).t, replay, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Reverse Gain','Visual Pattern Velocity');
title('Reverse Gain (red, incl. hand-picked) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');



%% Plot dff, forward, rotation, and pattern velocity with pattern pause highlighted in red blocks (each trace on its own subplot)

figure('Position',[100,100,900,800]);
varnames = {'dff', 'ball_forward_velocity', 'ball_yaw_velocity', 'visual_pattern_velocity', 'visual_heading'};
labels = {'dF/F', 'Fwd velocity', 'Yaw velocity', 'Pattern velocity', 'Pattern heading'};
colors = {'b', 'g', 'c', 'm', 'k'};
signals = {dff, ball_forward_velocity, ball_yaw_velocity, visual_pattern_velocity, visual_heading};

% Get pause start/stop intervals
pause_on = [false; diff(pauses(:))]==1;
pause_off = [false; diff(pauses(:))]==-1;
pause_starts = find(pause_on);
pause_ends = find(pause_off);

if pauses(1)==1, pause_starts = [1; pause_starts]; end
if pauses(end)==1, pause_ends = [pause_ends; length(pauses)]; end

% Get open loop start/stop intervals
openloop_on = [false; diff(open_loop_f2b(:))]==1;
openloop_off = [false; diff(open_loop_f2b(:))]==-1;
openloop_starts = find(openloop_on);
openloop_ends = find(openloop_off);

if open_loop_f2b(1)==1, openloop_starts = [1; openloop_starts]; end
if open_loop_f2b(end)==1, openloop_ends = [openloop_ends; length(open_loop_f2b)]; end


% Get reverse gain start/stop intervals
reversegain_on = [false; diff(reverse_gain(:))]==1;
reversegain_off = [false; diff(reverse_gain(:))]==-1;
reversegain_starts = find(reversegain_on);
reversegain_ends = find(reversegain_off);

if reverse_gain(1)==1, reversegain_starts = [1; reversegain_starts]; end
if reverse_gain(end)==1, reversegain_ends = [reversegain_ends; length(reverse_gain)]; end

% Get replay start/stop intervals
replay_on = [false; diff(replay(:))]==1;
replay_off = [false; diff(replay(:))]==-1;
replay_starts = find(replay_on);
replay_ends = find(replay_off);

if reverse_gain(1)==1, replay_starts = [1; replay_starts]; end
if reverse_gain(end)==1, replay_ends = [replay_ends; length(replay)]; end

ax = gobjects(numel(signals),1);  % Pre-allocate axis handles
for k = 1:numel(signals)
    ax(k) = subplot(numel(signals),1,k); hold on;
    y = signals{k};
    yl = [min(y)-0.1*range(y), max(y)+0.1*range(y)];
    % Plot red blocks for pauses
    for i = 1:length(pause_starts)
        patch([time(pause_starts(i)), time(pause_ends(i)), time(pause_ends(i)), time(pause_starts(i))], ...
            [yl(1), yl(1), yl(2), yl(2)], [1 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    for i = 1:length(openloop_starts)
        patch([time(openloop_starts(i)), time(openloop_ends(i)), time(openloop_ends(i)), time(openloop_starts(i))], ...
            [yl(1), yl(1), yl(2), yl(2)], [0.8 1 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    for i = 1:length(reversegain_starts)
        patch([time(reversegain_starts(i)), time(reversegain_ends(i)), time(reversegain_ends(i)), time(reversegain_starts(i))], ...
            [yl(1), yl(1), yl(2), yl(2)], [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    for i = 1:length(replay_starts)
        patch([time(replay_starts(i)), time(replay_ends(i)), time(replay_ends(i)), time(replay_starts(i))], ...
            [yl(1), yl(1), yl(2), yl(2)], [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    % Plot the trace
    plot(time, y, 'Color', colors{k}, 'LineWidth', 1.2);
    xlim([time(1), time(end)]);
    ylim(yl);
    if k == numel(signals)
        xlabel('Time (s)');
    end
    ylabel(labels{k});
    if k == 1
        title('dF/F & Kinematics with Pattern Pause (red blocks)');
    end
end
set(gcf,'Color','w');
linkaxes(ax, 'x');  % Only link the X axes! 

% Optional: Add legend only to the first subplot
subplot(numel(signals),1,1);
save_plot_with_title_as_filename('all_variables', 'pauses', savepath);


%% scatter plot version
% Logical mask for speeds > 0
select = ball_forward_velocity > abs(0);

figure; hold on
scatter(ball_forward_velocity(pauses==0 & select & open_loop_f2b ==0 & reverse_gain ==0), dff(pauses==0 & select & open_loop_f2b ==0 & reverse_gain ==0), 20, 'b', 'filled', 'DisplayName', 'Run');
scatter(ball_forward_velocity(pauses==1 & select), dff(pauses==1 & select), 20, 'r', 'filled', 'DisplayName', 'Pause');
scatter(ball_forward_velocity(open_loop_f2b==1 & select), dff(open_loop_f2b==1 & select), 20, 'g', 'filled', 'DisplayName', 'Open Loop');
scatter(ball_forward_velocity(reverse_gain==1 & select), dff(reverse_gain==1 & select), 20, 'm', 'filled', 'DisplayName', 'Reverse gain');
scatter(ball_forward_velocity(replay==1 & select), dff(replay==1 & select), 20, 'yellow', 'filled', 'DisplayName', 'Reverse gain');
xlabel('Forward velocity');
ylabel('dF/F');
legend('Location','best')
title('dF/F vs Forward Velocity');
grid on
set(gcf, 'Color', 'w');
save_plot_with_title_as_filename('scatterplot', 'pauses_openloop_reversegain', savepath);

%% lm for forward speed and dff
% Logical mask for speeds > 1
select = ball_forward_velocity > 1;

% Group masks
idx_closed_loop = (pauses == 0) & (open_loop_f2b == 0) & select & (reverse_gain==0)& (replay==0);
idx_pause = (pauses == 1) & (open_loop_f2b == 0) & select;
idx_open_loop = (pauses == 0) & (open_loop_f2b == 1) & select;
idx_reverse_gain = (pauses == 0) & (reverse_gain == 1) & select;
idx_replay = (pauses == 0) & (replay == 1) & select;
figure; hold on

% Scatter data
scatter(ball_forward_velocity(idx_closed_loop), dff(idx_closed_loop), 20, 'b', 'filled', 'DisplayName', 'Gain On');
scatter(ball_forward_velocity(pauses==1 & select), dff(pauses==1 & select), 20, 'r', 'filled', 'DisplayName', 'Pause');
scatter(ball_forward_velocity(open_loop_f2b==1 & select), dff(open_loop_f2b==1 & select), 20, 'g', 'filled', 'DisplayName', 'Open Loop');
scatter(ball_forward_velocity(reverse_gain==1 & select), dff(reverse_gain==1 & select), 20, 'm', 'filled', 'DisplayName', 'Reverse gain');
scatter(ball_forward_velocity(replay==1 & select), dff(replay==1 & select), 20, 'yellow', 'filled', 'DisplayName', 'Replay');

% Linear fit for 'run'
if sum(idx_closed_loop) > 1
    x_run = ball_forward_velocity(idx_closed_loop);
    y_run = dff(idx_closed_loop);
    p_run = polyfit(x_run, y_run, 1);
    xfit_run = linspace(min(x_run), max(x_run), 100);
    yfit_run = polyval(p_run, xfit_run);
    plot(xfit_run, yfit_run, 'b-', 'LineWidth', 2, 'DisplayName','Run fit');
    
    % R^2 calculation
    y_pred_run = polyval(p_run, x_run);
    ss_res_run = sum((y_run - y_pred_run).^2);
    ss_tot_run = sum((y_run - mean(y_run)).^2);
    R2_run = 1 - ss_res_run/ss_tot_run;
    % Annotate
    txt1 = sprintf('Run R^2 = %.2f', R2_run);
else
    txt1 = 'Run R^2 = N/A';
end

% Linear fit for 'pause'
if sum(idx_pause) > 1
    x_pause = ball_forward_velocity(idx_pause);
    y_pause = dff(idx_pause);
    p_pause = polyfit(x_pause, y_pause, 1);
    xfit_pause = linspace(min(x_pause), max(x_pause), 100);
    yfit_pause = polyval(p_pause, xfit_pause);
    plot(xfit_pause, yfit_pause, 'r-', 'LineWidth', 2, 'DisplayName','Pause fit');

    % R^2 calculation
    y_pred_pause = polyval(p_pause, x_pause);
    ss_res_pause = sum((y_pause - y_pred_pause).^2);
    ss_tot_pause = sum((y_pause - mean(y_pause)).^2);
    R2_pause = 1 - ss_res_pause/ss_tot_pause;
    % Annotate
    txt2 = sprintf('Pause R^2 = %.2f', R2_pause);
else
    txt2 = 'Pause R^2 = N/A';
end


% Linear fit for 'open loop'
if sum(idx_open_loop) > 1
    x_ol = ball_forward_velocity(idx_open_loop);
    y_ol = dff(idx_open_loop);
    p_ol = polyfit(x_ol, y_ol, 1);
    xfit_ol = linspace(min(x_ol), max(x_ol), 100);
    yfit_ol = polyval(p_ol, xfit_ol);
    plot(xfit_ol, yfit_ol, 'g-', 'LineWidth', 2, 'DisplayName','Open Loop fit');

    % R^2 calculation
    y_pred_ol = polyval(p_ol, x_ol);
    ss_res_ol = sum((y_ol - y_pred_ol).^2);
    ss_tot_ol = sum((y_ol - mean(y_ol)).^2);
    R2_ol = 1 - ss_res_ol/ss_tot_ol;
    % Annotate
    txt3 = sprintf('Open Loop R^2 = %.2f', R2_ol);
else
    txt3 = 'Open Loop R^2 = N/A';
end
% lm reverse gain
if sum(idx_reverse_gain) > 1
    x_rg = ball_forward_velocity(idx_reverse_gain);
    y_rg = dff(idx_reverse_gain);
    p_rg = polyfit(x_rg, y_rg, 1);
    xfit_rg = linspace(min(x_rg), max(x_rg), 100);
    yfit_rg = polyval(p_rg, xfit_rg);
    plot(xfit_rg, yfit_rg, 'm-', 'LineWidth', 2, 'DisplayName','RG fit');
    
    % R^2 calculation
    y_pred_rg = polyval(p_rg, x_rg);
    ss_res_rg = sum((y_rg - y_pred_rg).^2);
    ss_tot_rg = sum((y_rg - mean(y_rg)).^2);
    R2_rg = 1 - ss_res_rg/ss_tot_rg;
    % Annotate
    txt4 = sprintf('Reverse gain R^2 = %.2f', R2_rg);
else
    txt4 = 'RG R^2 = N/A';
end
if sum(idx_replay) > 1
    x_replay = ball_forward_velocity(idx_replay);
    y_replay = dff(idx_replay);
    p_replay = polyfit(x_replay, y_replay, 1);
    xfit_replay = linspace(min(x_replay), max(x_replay), 100);
    yfit_replay = polyval(p_replay, xfit_replay);
    plot(xfit_replay, yfit_replay, 'yellow-', 'LineWidth', 2, 'DisplayName','Replay fit');
    
    % R^2 calculation
    y_pred_replay = polyval(p_replay, x_replay);
    ss_res_replay = sum((y_replay - y_pred_replay).^2);
    ss_tot_replay = sum((y_replay - mean(y_replay)).^2);
    R2_replay = 1 - ss_res_replay/ss_tot_replay;
    % Annotate
    txt5 = sprintf('Replay R^2 = %.2f', R2_replay);
else
    txt5 = 'RG R^2 = N/A';
end
xlabel('Forward velocity');
ylabel('dF/F');
title('dF/F vs Forward Velocity (>1, red=pause)');
legend('Location','best')
grid on
set(gcf, 'Color', 'w');

% R^2 annotation
% Place text in top left and top right corners of the plot
yl = ylim; xl = xlim;
text(xl(1) + 0.03*range(xl), yl(2) - 0.07*range(yl), txt1, 'Color','b','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.15*range(yl), txt2, 'Color','r','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.21*range(yl), txt3, 'Color','g','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.28*range(yl), txt4, 'Color','m','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.35*range(yl), txt5, 'Color','yellow','FontWeight','bold','FontSize',12);

save_plot_with_title_as_filename('scatterplot', 'pauses_ol_lm', savepath);

%% lm for visual pattern and dff

% Group masks
idx_closed_loop = (pauses == 0) & (open_loop_f2b == 0) & (reverse_gain == 0) & (replay==0);
idx_pause = (pauses == 1) & (open_loop_f2b == 0);
idx_open_loop = (pauses == 0) & (open_loop_f2b == 1);
idx_reverse_gain = (reverse_gain ==1);
idx_replay = (replay ==1);
figure; hold on

% Scatter data
scatter(visual_pattern_velocity(idx_closed_loop),   dff(idx_closed_loop),   20, 'b', 'filled', 'DisplayName', 'Run');
scatter(visual_pattern_velocity(idx_pause), dff(idx_pause), 20, 'r', 'filled', 'DisplayName', 'Pause');
scatter(visual_pattern_velocity(idx_open_loop), dff(idx_open_loop), 20, 'g', 'filled', 'DisplayName', 'Open Loop');
scatter(visual_pattern_velocity(idx_reverse_gain), dff(idx_reverse_gain), 20, 'm', 'filled', 'DisplayName', 'Reverse gain');
scatter(visual_pattern_velocity(idx_replay), dff(idx_replay), 20, 'yellow', 'filled', 'DisplayName', 'Reverse gain');


xlabel('Pattern Velocity');
ylabel('dF/F');
title('dF/F vs Visual Pattern Velocity');
%legend('Location','best')
grid on
set(gcf, 'Color', 'w');


save_plot_with_title_as_filename('scatterplot', 'visual_pattern_pauses_ol_lm', savepath);

end