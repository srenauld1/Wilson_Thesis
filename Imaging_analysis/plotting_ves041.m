
%% function
function plotting_ves041(ts, savepath)

arguments
    ts;
    savepath string = "/Users/sophiarenauld/stacks/20250304-4_ves041-1_fly1_clbarjumps_1";  % Default value
end

close all

% reformat image data
if isfield(ts, 'resp') && isfield(ts.resp, 'i2') && iscell(ts.resp.i2)
    if ~isempty(ts.resp.i2) && numel(ts.resp.i2) >= 1
        % Extract the data from the first cell of ts.resp.i2 and assign it to ts.dff
        ts.dff = ts.resp.i2{1};
    else
        error('ts.resp.i2 is empty or does not have enough elements.');
    end
else
    disp('ts.resp.i2 is not correctly defined.');
end

% Remove the ts.resp field
if isfield(ts, 'resp')
    ts = rmfield(ts, 'resp');
else
    disp('The ts.resp field does not exist and cannot be removed.');
end
[ts, triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(ts, savepath, 10, 2,0.88,0.88);
%% smooth fictrac
fwd = ts.ball.forvel;
yaw = ts.ball.yawvel;
side = ts.ball.sidevel;
dff = ts.dff;
time=ts.t;
%dff = ts.resp.i1.imf_f_dff005000_n;

%% Plot the dff and fwd
% Plot the dff and fwd
figure;
yyaxis left;
plot(time, fwd, '-b');
plot(time, fwd, '-b', 'LineWidth', 3)% Plot first time series in blue
ylabel('Forward Velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 10]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-r' , 'LineWidth', 3);  % Plot second time series in red
ylabel('dff');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('Forward velocity and DFF');

% Optional: Add legend
legend('Forward Velocity', 'dff');
save_plot_with_title_as_filename('Fwd vel', 'dff', savepath);

%% Plot the dff and fwd CONTROL
% Plot the dff and fwd
% figure;
% yyaxis left;
% plot(fwd, '-b');  % Plot first time series in blue
% ylabel('forward velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-5 10]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(ts.resp.control.moex_0_0_0.rawf_f_dff008000_n, '-r');  % Plot second time series in red
% ylabel('dff');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('Forward velocity and dff - control');
% 
% % Optional: Add legend
% legend('Forward velocity', 'dff - control');
% save_plot_with_title_as_filename('Fwd vel', 'control dff', savepath);


%% plot dff and yaw
figure;
yyaxis left;
plot(time, yaw, '-g', 'LineWidth', 3);  % Plot first time series in blue
ylabel('Yaw velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 5]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-r', 'LineWidth', 3);  % Plot second time series in red
ylabel('dff');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('Yaw velocity and DFF');

% Optional: Add legend
legend('yaw velocity', 'dff');
save_plot_with_title_as_filename('Yaw vel', 'dff', savepath);

%% plot fwd and yaw
figure;
yyaxis left;
plot(time, fwd, '-b');  % Plot first time series in blue
ylabel('Forward velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 5]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, yaw, '-g');  % Plot second time series in red
ylabel('yaw vel');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('forward and yaw velocity');

% Optional: Add legend
legend('fwd velocity', 'yaw velocity');
save_plot_with_title_as_filename('Fwd vel', 'Yaw vel', savepath);


%% Plot the side and dff
% Plot the dff and fwd
figure;
yyaxis left;
plot(time, side, '-y');  % Plot first time series in blue
ylabel('Sideways velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 10]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-r');  % Plot second time series in red
ylabel('dff');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('Side velocity vs dff');

% Optional: Add legend
legend('Sideways velocity', 'dff');
save_plot_with_title_as_filename('Side vel', 'dff', savepath);

%% flat path
figure
plot(ts.flypos.x, ts.flypos.y)
hold on
plot(ts.flypos.x(1), ts.flypos.y(1), 'r.', 'MarkerSize', 20)  % The MarkerSize value controls dot size
save_plot_with_title_as_filename('x', 'y', savepath);


%% flat path colored by dff
cmap=jet(100);
dff = dff;
dff_norm = round((dff-min(dff))/(max(dff)-min(dff))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(ts.flypos.x)-1
    plot(ts.flypos.x(i:i+1), ts.flypos.y(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
end

% Create colorbar with actual DFF values
c = colorbar;
caxis([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual DFF values
actual_values = yticks/100 * (max(dff)-min(dff)) + min(dff);
% Format tick labels to 3 decimal places
ylabel(c, 'DFF')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))

hold on
% Add jump points
% plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
%      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
%      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
plot(ts.flypos.x(1), ts.flypos.y(1), 'r.', 'MarkerSize', 20)  % Start point

title("DFF colored, cue flat path")
save_plot_with_title_as_filename('x_color', 'y_color', savepath);
%% flat path colored by dff - ball
figure
hold on
for i = 1:length(ts.flypos.x)-1
    plot(ts.flypos.x_ball(i:i+1), ts.flypos.y_ball(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
end



% Create colorbar with actual DFF values
c = colorbar;
caxis([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual DFF values
actual_values = yticks/100 * (max(dff)-min(dff)) + min(dff);
% Format tick labels to 3 decimal places
ylabel(c, 'DFF')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))

hold on
% Add jump points
% plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
%      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
%      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
plot(ts.flypos.x_ball(1), ts.flypos.y_ball(1), 'r.', 'MarkerSize', 20)  % Start point

title("DFF colored, ball flat path")
save_plot_with_title_as_filename('x_ball_color', 'y_ball_color', savepath);

%% revamp variables

total_mov_mm = (abs(fwd) + abs(ts.ball.sidevel) + abs(yaw));        
not_moving = total_mov_mm <= 3; 
speed = abs(yaw);
% Identify outliers using the default method (usually interquartile range)
outliers1 = isoutlier(fwd, 'mean', 'ThresholdFactor', 6);
outliers2 = isoutlier(speed, 'mean', 'ThresholdFactor', 6);
outliers = outliers1 | outliers2 | not_moving;

% Remove outliers from forward_velocity
cleaned_forward_velocity = fwd(~outliers);
cleaned_rot_speed = speed(~outliers);
cleaned_dff = dff(~outliers);
cleaned_xpos_ball = ts.flypos.x_ball(~outliers);
cleaned_ypos_ball = ts.flypos.y_ball(~outliers);

%% flat path colored by dff - ball
figure
hold on

dffclean_norm = round((cleaned_dff-min(cleaned_dff))/(max(cleaned_dff)-min(cleaned_dff))*(length(cmap)-1)) +1;
for i = 1:length(cleaned_xpos_ball)-1
    plot(cleaned_xpos_ball(i:i+1), cleaned_ypos_ball(i:i+1), 'Color', cmap(dffclean_norm(i), :), 'LineWidth', 2)
end



% Create colorbar with actual DFF values
c = colorbar;
caxis([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual DFF values
actual_values = yticks/100 * (max(dff)-min(dff)) + min(dff);
% Format tick labels to 3 decimal places
ylabel(c, 'DFF')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))

hold on
% Add jump points
% plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
%      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
%      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
plot(cleaned_xpos_ball(1), cleaned_ypos_ball(1), 'r.', 'MarkerSize', 20)  % Start point

title("DFF colored, ball flat path, total vel >3")
save_plot_with_title_as_filename('x_ball_color_moving', 'y_ball_color_moving', savepath);
%% flat path colored by fwd vel
cmap=jet(100);
fwd_norm = round((cleaned_forward_velocity-min(cleaned_forward_velocity))/(max(cleaned_forward_velocity)-min(cleaned_forward_velocity))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(cleaned_xpos_ball)-1
    plot(cleaned_xpos_ball(i:i+1), cleaned_ypos_ball(i:i+1), 'Color', cmap(fwd_norm(i), :), 'LineWidth', 2)
end


% Create colorbar with actual forward velocity values
c = colorbar;
caxis([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual forward velocity values
actual_values = yticks/100 * (max(cleaned_forward_velocity)-min(cleaned_forward_velocity)) + min(cleaned_forward_velocity);
% Format tick labels to 1 decimal place
ylabel(c, 'Forward Velocity (mm/s)')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), actual_values, 'UniformOutput', false))

hold on
% Add jump points
% plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
%      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
%      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
plot(cleaned_xpos_ball(1), cleaned_ypos_ball(1), 'r.', 'MarkerSize', 20)  % Start point

title("Forward Velocity colored")

save_plot_with_title_as_filename('x_fwd', 'y_fwd', savepath);%% flat path colored by fwd vel


title("Forward Velocity colored")
save_plot_with_title_as_filename('x_fwd', 'y_fwd', savepath);

%% flat path colored by rot vel
cmap=jet(100);
speed_norm = round((cleaned_rot_speed-min(cleaned_rot_speed))/(max(cleaned_rot_speed)-min(cleaned_rot_speed))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(cleaned_xpos_ball)-1
    plot(cleaned_xpos_ball(i:i+1), cleaned_ypos_ball(i:i+1), 'Color', cmap(speed_norm(i), :), 'LineWidth', 2)
end

% Create colorbar with actual rotation speed values
c = colorbar;
clim([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual rotation speed values
actual_values = yticks/100 * (max(cleaned_rot_speed)-min(cleaned_rot_speed)) + min(cleaned_rot_speed);
% Format tick labels to 1 decimal place
ylabel(c, 'Rotation Speed (mm/s)')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), actual_values, 'UniformOutput', false))

hold on
% Add jump points
% plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
%      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
%      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
plot(cleaned_xpos_ball(1), cleaned_ypos_ball(1), 'r.', 'MarkerSize', 20)  % Start point

title("Rotation Speed colored")

save_plot_with_title_as_filename('x_rot', 'y_rot', savepath);


%% all 3 on one plot
time = ts.t;                  % Time points
dff = dff;           % Example dFF data
fwd = fwd; % Example forward velocity data
rot = yaw; % Example rotational velocity data

% Create the figure and subplots
figure;

% Subplot 1: dFF
subplot(3, 1, 1);  % First subplot
hold on
plot(time, dff, 'r', 'LineWidth', 2);  % Plot dFF in red
xline(150.363);
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F');
grid on;

% Subplot 2: Forward Velocity
subplot(3, 1, 2);  % Second subplot
plot(time, fwd, 'g', 'LineWidth', 2);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;

% Subplot 3: Rotational Velocity
subplot(3, 1, 3);  % Third subplot
plot(time, rot, 'b', 'LineWidth', 2);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (mm/s)');
title('Rotational Velocity');
grid on;

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);


