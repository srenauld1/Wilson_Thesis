function a2p_data = process_bar_mp_data(a2p_data, savepath, box, split)

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

%% segment out pauses;
bar_on = visual_heading < 0;

figure
hold on
plot(a2p_data.dq(1).t, bar_on, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_heading, 'b')
xlabel('Time (s)');
legend('Bar on (<0)', 'Visual Heading');
title('Bar on (red) and Visual Heading (blue)');
set(gcf, 'Color', 'w');

%% segment out gain changes

% gain = 0.8
% Initialize gain = 0.8 to 0, or load if it already exists
if ~isfield(a2p_data, 'gain8') || isempty(a2p_data.gain8)
    gain8 = zeros(size(visual_heading));
else
    gain8 = a2p_data.gain8;
end

% Initialize gain = 0.4 to 0, or load if it already exists
if ~isfield(a2p_data, 'gain4') || isempty(a2p_data.gain4)
    gain4 = zeros(size(visual_heading));
else
    gain4 = a2p_data.gain4;
end

% Initialize gain = 1.6 to 0, or load if it already exists
if ~isfield(a2p_data, 'gain16') || isempty(a2p_data.gain16)
    gain16 = zeros(size(visual_heading));
else
    gain16 = a2p_data.gain16;
end


figure; hold on
plot(a2p_data.dq(1).t, gain8, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_heading, 'b')
xlabel('Time (s)');
legend('Gain 0.8','Visual Pattern Velocity');
title('Gain 0.8 (red) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');

disp('Zoom/pan as needed, then click and drag to draw a rectangle for each 0.8 gain interval.');
disp('When you are finished, just single-click on the plot and press Enter without drawing a rectangle.');

while true
    h = drawrectangle('Label','0.8 gain','Color','m');
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
    gain8(indices) = 1;
    % Visually mark selected interval
    plot(a2p_data.dq(1).t(indices), visual_heading(indices),...
        'ko', 'MarkerFaceColor','y','DisplayName','Hand-picked Reverse Gain');
    delete(h);
end


% Show the final result for review
figure; hold on
plot(a2p_data.dq(1).t, gain8, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_pattern_velocity, 'b')
xlabel('Time (s)');
legend('Reverse Gain','Visual Pattern Velocity');
title('Reverse Gain (red, incl. hand-picked) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');


%% segment out gain 4 (by hand)
figure; hold on
plot(a2p_data.dq(1).t, gain4, 'r', 'LineWidth',1.5)
plot(a2p_data.dq(1).t, visual_heading, 'b')
xlabel('Time (s)');
legend('Gain 0.4','Visual Pattern Velocity');
title('Gain 0.4 (red) and Pattern Velocity (blue)');
set(gcf, 'Color', 'w');

disp('Zoom/pan as needed, then click and drag to draw a rectangle for each 0.4 gain interval.');
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
    gain4(indices) = 1;
    % Visually mark selected interval
    plot(a2p_data.dq(1).t(indices), visual_heading(indices),...
        'ko', 'MarkerFaceColor','y','DisplayName','Hand-picked Reverse Gain');
    delete(h);
end

%% save data to a2p
a2p_data.no_bar = bar_on;
a2p_data.gain8 = gain8;
a2p_data.gain16 = gain16;
a2p_data.gain4 = gain4;


%% Plot dff, forward, rotation, and pattern velocity with pattern pause highlighted in red blocks (each trace on its own subplot)

figure('Position',[100,100,900,800]);
varnames = {'dff', 'ball_forward_velocity', 'ball_yaw_velocity', 'visual_heading', 'visual_pattern_velocity'};
labels = {'dF/F', 'Fwd velocity', 'Yaw velocity',  'Pattern heading', 'Pattern velocity'};
colors = {'b', 'g', 'c', 'm', 'k'};
signals = {dff, ball_forward_velocity, ball_yaw_velocity, visual_heading, visual_pattern_velocity};

% NOTE: all these variables must be logical vectors (0/1, same length/timebase)

cat_bar_on_8 = (bar_on==1) & (gain8==1);
cat_bar_on_4 = (bar_on==1) & (gain4==1);
cat_bar_off_8 = (bar_on==0) & (gain8==1);
cat_bar_off_4 = (bar_on==0) & (gain4==1);

% Helper function to get interval indices for a logical mask
get_intervals = @(mask) deal( ...
    find([false; diff(mask(:))]==1), ...       % starts
    find([false; diff(mask(:))]==-1) ...       % ends
);

[on8_starts, on8_ends]   = get_intervals(cat_bar_on_8);
[on4_starts, on4_ends]   = get_intervals(cat_bar_on_4);
[off8_starts, off8_ends] = get_intervals(cat_bar_off_8);
[off4_starts, off4_ends] = get_intervals(cat_bar_off_4);

% Edge correction for each
if cat_bar_on_8(1)==1,   on8_starts = [1; on8_starts]; end
if cat_bar_on_8(end)==1, on8_ends   = [on8_ends; length(cat_bar_on_8)]; end

if cat_bar_on_4(1)==1,   on4_starts = [1; on4_starts]; end
if cat_bar_on_4(end)==1, on4_ends   = [on4_ends; length(cat_bar_on_4)]; end

if cat_bar_off_8(1)==1,  off8_starts = [1; off8_starts]; end
if cat_bar_off_8(end)==1,off8_ends   = [off8_ends; length(cat_bar_off_8)]; end

if cat_bar_off_4(1)==1,  off4_starts = [1; off4_starts]; end
if cat_bar_off_4(end)==1,off4_ends   = [off4_ends; length(cat_bar_off_4)]; end


ax = gobjects(numel(signals),1);  % Pre-allocate axis handles
for k = 1:numel(signals)
    ax(k) = subplot(numel(signals),1,k); hold on;
    y = signals{k};
    yl = [min(y)-0.1*range(y), max(y)+0.1*range(y)];
    % Use the actual y-limits for the patch
    %yl = ylim; % or set ymin/ymax if you want identical across subplots
    
    for i = 1:length(on8_starts)
        patch([time(on8_starts(i)), time(on8_ends(i)), time(on8_ends(i)), time(on8_starts(i))], ...
              [yl(1) yl(1) yl(2) yl(2)], [0 0.4 1], 'EdgeColor','none', 'FaceAlpha',0.25);
    end
    for i = 1:length(on4_starts)
        patch([time(on4_starts(i)), time(on4_ends(i)), time(on4_ends(i)), time(on4_starts(i))], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0.4 0.7], 'EdgeColor','none', 'FaceAlpha',0.25);
    end
    for i = 1:length(off8_starts)
        patch([time(off8_starts(i)), time(off8_ends(i)), time(off8_ends(i)), time(off8_starts(i))], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0 0], 'EdgeColor','none', 'FaceAlpha',0.25);
    end
    for i = 1:length(off4_starts)
        patch([time(off4_starts(i)), time(off4_ends(i)), time(off4_ends(i)), time(off4_starts(i))], ...
              [yl(1) yl(1) yl(2) yl(2)], [0 0.9 0], 'EdgeColor','none', 'FaceAlpha',0.25);
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

subplot(numel(signals),1,5); hold on;  % Activate first subplot
h_on8  = patch([nan nan nan nan], [nan nan nan nan], [0 0.4 1],      'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','Bar ON, Gain 8');
h_on4  = patch([nan nan nan nan], [nan nan nan nan], [1 0.4 0.7],    'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','Bar ON, Gain 4');
h_off8 = patch([nan nan nan nan], [nan nan nan nan], [1 0 0],        'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','Bar OFF, Gain 8');
h_off4 = patch([nan nan nan nan], [nan nan nan nan], [0 0.9 0],      'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','Bar OFF, Gain 4');
legend([h_on8 h_on4 h_off8 h_off4], 'Location', 'northeastoutside', 'Box','off');
save_plot_with_title_as_filename('all_variables', 'colorcoded', savepath);



%% scatter plot version
% Logical mask for speeds > 0
select = ball_forward_velocity > 0;

figure; hold on
scatter(ball_forward_velocity(bar_on==1 & select & gain8==1), dff(bar_on==1 & select & gain8==1), 20, 'b', 'filled', 'DisplayName', 'bar, gain 8');
scatter(ball_forward_velocity(gain8==1 & select & bar_on ==0), dff(gain8==1 & select & bar_on ==0), 20, 'r', 'filled', 'DisplayName', 'no bar, gain 8');
scatter(ball_forward_velocity(gain4==1 & select & bar_on ==0), dff(gain4==1 & select & bar_on ==0), 20, 'g', 'filled', 'DisplayName', 'no bar, gain 4');
scatter(ball_forward_velocity(gain4==1 & select & bar_on ==1), dff(gain4==1 & select & bar_on ==1), 20, 'm', 'filled', 'DisplayName', 'bar, gain 0.4');
xlabel('Forward velocity');
ylabel('dF/F');
legend('Location','best')
title('dF/F vs Forward Velocity');
grid on
set(gcf, 'Color', 'w');
save_plot_with_title_as_filename('scatterplot', 'gains_and_bar', savepath);

%% lm for forward speed and dff
% Logical mask for speeds > 1
select = abs(ball_forward_velocity) > 0;

% Group masks
idx_no_bar8 = (bar_on == 0) & (gain8 == 1);
idx_bar_8 = (bar_on == 1) & (gain8 == 1) & select;
idx_bar4 = (bar_on == 1) & (gain4 == 1) & select;
idx_no_bar4 = (bar_on == 0) & (gain4 == 1) & select;
figure; hold on

% Scatter data
scatter(ball_forward_velocity(bar_on==1 & select & gain8==1), dff(bar_on==1 & select & gain8==1), 20, 'b', 'filled', 'DisplayName', 'bar on, gain 8');
scatter(ball_forward_velocity(gain8==1 & select & bar_on ==0), dff(gain8==1 & select & bar_on ==0), 20, 'r', 'filled', 'DisplayName', 'bar off, gain 8');
scatter(ball_forward_velocity(gain4==1 & select & bar_on ==0), dff(gain4==1 & select & bar_on ==0), 20, 'g', 'filled', 'DisplayName', 'bar off, gain 4');
scatter(ball_forward_velocity(gain4==1 & select & bar_on ==1), dff(gain4==1 & select & bar_on ==1), 20, 'm', 'filled', 'DisplayName', 'bar on, gain 0.4');
% Linear fit for no bar 0.8
if sum(idx_no_bar8) > 1
    x_run = ball_forward_velocity(idx_no_bar8);
    y_run = dff(idx_no_bar8);
    p_run = polyfit(x_run, y_run, 1);
    xfit_run = linspace(min(x_run), max(x_run), 100);
    yfit_run = polyval(p_run, xfit_run);
    plot(xfit_run, yfit_run, 'r-', 'LineWidth', 2, 'DisplayName','No Bar, 0.8');
    
    % R^2 calculation
    y_pred_run = polyval(p_run, x_run);
    ss_res_run = sum((y_run - y_pred_run).^2);
    ss_tot_run = sum((y_run - mean(y_run)).^2);
    R2_run = 1 - ss_res_run/ss_tot_run;
    % Annotate
    txt1 = sprintf('No bar, 0.8 R^2 = %.2f', R2_run);
else
    txt1 = 'No bar 0.8 R^2 = N/A';
end

% Linear fit for 0.8 bar
if sum(idx_bar_8) > 1
    x_pause = ball_forward_velocity(idx_bar_8);
    y_pause = dff(idx_bar_8);
    p_pause = polyfit(x_pause, y_pause, 1);
    xfit_pause = linspace(min(x_pause), max(x_pause), 100);
    yfit_pause = polyval(p_pause, xfit_pause);
    plot(xfit_pause, yfit_pause, 'b-', 'LineWidth', 2, 'DisplayName','Bar, 0.8 fit');

    % R^2 calculation
    y_pred_pause = polyval(p_pause, x_pause);
    ss_res_pause = sum((y_pause - y_pred_pause).^2);
    ss_tot_pause = sum((y_pause - mean(y_pause)).^2);
    R2_pause = 1 - ss_res_pause/ss_tot_pause;
    % Annotate
    txt2 = sprintf('Bar, 0.8 R^2 = %.2f', R2_pause);
else
    txt2 = 'Bar 0.8 R^2 = N/A';
end


% Linear fit for 0.4 bar
if sum(idx_bar4) > 1
    x_ol = ball_forward_velocity(idx_bar4);
    y_ol = dff(idx_bar4);
    p_ol = polyfit(x_ol, y_ol, 1);
    xfit_ol = linspace(min(x_ol), max(x_ol), 100);
    yfit_ol = polyval(p_ol, xfit_ol);
    plot(xfit_ol, yfit_ol, 'm-', 'LineWidth', 2, 'DisplayName','Bar, 0.4 fit');

    % R^2 calculation
    y_pred_ol = polyval(p_ol, x_ol);
    ss_res_ol = sum((y_ol - y_pred_ol).^2);
    ss_tot_ol = sum((y_ol - mean(y_ol)).^2);
    R2_ol = 1 - ss_res_ol/ss_tot_ol;
    % Annotate
    txt3 = sprintf('Bar, 0.4 R^2 = %.2f', R2_ol);
else
    txt3 = 'Bar, 0.4 R^2 = N/A';
end
% lm 0.4 no bar
if sum(idx_no_bar4) > 1
    x_rg = ball_forward_velocity(idx_no_bar4);
    y_rg = dff(idx_no_bar4);
    p_rg = polyfit(x_rg, y_rg, 1);
    xfit_rg = linspace(min(x_rg), max(x_rg), 100);
    yfit_rg = polyval(p_rg, xfit_rg);
    plot(xfit_rg, yfit_rg, 'g-', 'LineWidth', 2, 'DisplayName','No bar, 0.4 fit');
    
    % R^2 calculation
    y_pred_rg = polyval(p_rg, x_rg);
    ss_res_rg = sum((y_rg - y_pred_rg).^2);
    ss_tot_rg = sum((y_rg - mean(y_rg)).^2);
    R2_rg = 1 - ss_res_rg/ss_tot_rg;
    % Annotate
    txt4 = sprintf('No bar, 0.4 R^2 = %.2f', R2_rg);
else
    txt4 = 'No bar 0.4 R^2 = N/A';
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
text(xl(1) + 0.03*range(xl), yl(2) - 0.07*range(yl), txt1, 'Color','r','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.15*range(yl), txt2, 'Color','b','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.21*range(yl), txt3, 'Color','m','FontWeight','bold','FontSize',12);
text(xl(1) + 0.03*range(xl), yl(2) - 0.28*range(yl), txt4, 'Color','g','FontWeight','bold','FontSize',12);

save_plot_with_title_as_filename('scatterplot', 'gains_and_bar_lm', savepath);


%% look at time throughout gain period
% select gain period
gain_choice = input('What gain value are you choosing? (0 for all)');
bar_choice = input('Bar on (1) or off (0) or all (3)?');

% Logical mask for included points
select = abs(ball_forward_velocity) > 0;
if gain_choice == 8 || gain_choice==0.8
    gain_mask = gain8==1;
elseif gain_choice == 4 || gain_choice==0.4
    gain_mask = gain4==1;
elseif gain_choice == 16 || gain_choice==1.6
    gain_mask = gain16==1;
elseif gain_choice == 0
    gain_mask = true(size(gain8)); % selects ALL timepoints
end
t_all = a2p_data.dq(1).t;

% Find contiguous nonzero regions of gain8
d = diff([0; gain_mask(:); 0]);
chunk_starts = find(d == 1);
chunk_ends = find(d == -1) - 1;
valid = (chunk_ends - chunk_starts) > 0;
chunk_starts = chunk_starts(valid);
chunk_ends = chunk_ends(valid);

% Color gradients per chunk (expand as needed)
cmap_list = {
    [linspace(1,0,100)' linspace(0,0,100)' linspace(0,1,100)'];     % red→blue
    %[linspace(1,0.6,100)' linspace(0.4,1,100)' linspace(0.7,0,100)']; % pink→green
    [linspace(1,0.8,100)' linspace(1,0.5,100)' linspace(0,0.8,100)']; % yellow→purple
    [linspace(0,1,100)' linspace(1,0,100)' linspace(0,1,100)']        % cyan→magenta
    };

figure; hold on;
chunk_count = 0;
for c = 1:length(chunk_starts)
    idx_range = chunk_starts(c):chunk_ends(c);
    % Within this chunk, only plot where bar_on==1
    if bar_choice ==1
        bar_in_chunk = bar_on(idx_range)==1;
        bar_title='bar_on';
        bartitle='bar on';
    elseif bar_choice == 0
        bar_in_chunk = bar_on(idx_range)==0;
        bar_title = 'bar_off';
        bartitle='bar off';
    elseif bar_choice == 3
        bar_in_chunk = true(size(idx_range));  % include all points in this chunk
        bar_title = 'bar_on_and_off';
        bartitle='bar on and off';
    end
    if ~any(bar_in_chunk), continue; end

    chunk_count = chunk_count+1;
    cmap = cmap_list{mod(chunk_count-1, numel(cmap_list))+1};
    % Pull values
    x = ball_forward_velocity(idx_range);
    y = dff(idx_range);
    tvals = t_all(idx_range);
    % Keep only bar_on==1 points
    x = x(bar_in_chunk);
    y = y(bar_in_chunk);
    tvals = tvals(bar_in_chunk);

    % order by time for smooth color transition
    [t_order, sort_idx] = sort(tvals);
    x = x(sort_idx);
    y = y(sort_idx);

    npoints = length(x);
    chunk_cmap = cmap(round(linspace(1, size(cmap,1), npoints)), :);

    scatter(x, y, 20, chunk_cmap, 'filled');
end

xlabel('Forward velocity');
ylabel('dF/F');
plot_title = [bar_title '_gain_choice_' num2str(gain_choice)];
figure_title = [bartitle ' gain choice:' num2str(gain_choice)];
title(figure_title);
set(gcf, 'Color', 'w');
save_plot_with_title_as_filename('cl_scatterplot', plot_title, savepath)

end