
%% function
function plotting_an4(daq, dat, jump, savepath)

close all


%[daq, triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(daq, savepath, 10, 2,0.88,0.88);
%% smooth fictrac
fwd = daq.bvf;
yaw = daq.bvy_deg;
side = daq.bvs_deg;
dff = dat.ts(1,:);
time=daq.t;

%% old granular plots
% %% Plot the dff and fwd
% % Plot the dff and fwd
% figure;
% yyaxis left;
% plot(time, fwd, '-b');
% plot(time, fwd, '-b', 'LineWidth', 3)% Plot first time series in blue
% ylabel('Forward Velocity');  % Label for the left y-axis
% 
% % Set the y-axis limits for the left axis to [-5, 10]
% ylim([-5 10]);
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, dff, '-r' , 'LineWidth', 3);  % Plot second time series in red
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
% plot(time, yaw, '-g', 'LineWidth', 3);  % Plot first time series in blue
% ylabel('Yaw velocity');  % Label for the left y-axis
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, dff, '-r', 'LineWidth', 3);  % Plot second time series in red
% ylabel('dff');  % Label for the right y-axis
% 
% % Add title and x-axis label
% xlabel('Time');
% title('Yaw velocity and DFF');
% 
% % Optional: Add legend
% legend('yaw velocity', 'dff');
% save_plot_with_title_as_filename('Yaw vel', 'dff', savepath);
% 
% %% plot fwd and yaw
% figure;
% yyaxis left;
% plot(time, fwd, '-b');  % Plot first time series in blue
% ylabel('Forward velocity');  % Label for the left y-axis
% 
% % Plot the second time series on the right y-axis
% yyaxis right;
% plot(time, yaw, '-g');  % Plot second time series in red
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
% plot(time, side, '-y');  % Plot first time series in blue
% ylabel('Sideways velocity');  % Label for the left y-axis
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

%% all 3 on one plot
time = daq.t;                  % Time points
dff = dff;           % Example dFF data
fwd = fwd; % Example forward velocity data
rot = yaw; % Example rotational velocity data


% Create the figure and subplots
figure;

% Subplot 1: dFF
subplot(3, 1, 1);  % First subplot
hold on
plot(time, dff, 'b', 'LineWidth', 2);  % Plot dFF in red
xline(150.363);
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F');
grid on;

% Subplot 2: Forward Velocity
subplot(3, 1, 2);  % Second subplot
plot(time, fwd, 'black', 'LineWidth', 2);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;

% Subplot 3: Rotational Velocity
subplot(3, 1, 3);  % Third subplot
plot(time, rot, 'r', 'LineWidth', 2);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (rad/s)');
title('Rotational Velocity');
grid on;

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);


%% flat path
figure
plot(daq.px, daq.py)
hold on
plot(daq.px(1), daq.py(1), 'r.', 'MarkerSize', 20)  % The MarkerSize value controls dot size
save_plot_with_title_as_filename('x', 'y', savepath);


%% flat path colored by dff
cmap=jet(100);
dff_norm = round((dff-min(dff))/(max(dff)-min(dff))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(daq.px)-1
    plot(daq.px(i:i+1), daq.py(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
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
%Add jump points
if jump
    plot(daq.px(logical(daq.jump_detected)), ...
         daq.py(logical(daq.jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(daq.px(1), daq.py(1), 'r.', 'MarkerSize', 20)  % Start point

title("DFF colored, cue flat path")
save_plot_with_title_as_filename('x_color', 'y_color', savepath);
%% flat path colored by dff - ball
% figure
% hold on
% for i = 1:length(daq.px)-1
%     plot(daq.px_ball(i:i+1), daq.py_ball(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
% end
% 
% 
% 
% % Create colorbar with actual DFF values
% c = colorbar;
% caxis([1 100])  % Set colorbar limits to match normalized range
% yticks = get(c, 'YTick');  % Get current tick positions
% % Convert normalized values back to actual DFF values
% actual_values = yticks/100 * (max(dff)-min(dff)) + min(dff);
% % Format tick labels to 3 decimal places
% ylabel(c, 'DFF')
% set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))
% 
% hold on
% % Add jump points
% % plot(cleaned_xpos_ball(logical(ts.vis.jump_detected)), ...
% %      cleaned_ypos_ball(logical(ts.vis.jump_detected)), ...
% %      'm.', 'MarkerSize', 15)  % Adjust color and size as needed
% plot(daq.px_ball(1), daq.py_ball(1), 'r.', 'MarkerSize', 20)  % Start point
% 
% title("DFF colored, ball flat path")
% save_plot_with_title_as_filename('x_ball_color', 'y_ball_color', savepath);

%% revamp variables

total_mov_mm = daq.totalspeed;    
not_moving = daq.motion.moving_or_not;
speed = abs(yaw);
% Identify outliers using the default method (usually interquartile range)
outliers1 = isoutlier(fwd, 'mean', 'ThresholdFactor', 6);
outliers2 = isoutlier(speed, 'mean', 'ThresholdFactor', 6);
outliers = outliers1 | outliers2;

% Remove outliers from forward_velocity
cleaned_forward_velocity = fwd(~outliers);
cleaned_rot_speed = speed(~outliers);
cleaned_dff = dff(~outliers);
cleaned_xpos = daq.px(~outliers);
cleaned_ypos = daq.py(~outliers);

%% flat path colored by fwd vel
cmap=jet(100);
fwd_norm = round((cleaned_forward_velocity-min(cleaned_forward_velocity))/(max(cleaned_forward_velocity)-min(cleaned_forward_velocity))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(cleaned_xpos)-1
    plot(cleaned_xpos(i:i+1), cleaned_ypos(i:i+1), 'Color', cmap(fwd_norm(i), :), 'LineWidth', 2)
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
%Add jump points
if jump
    plot(daq.px(logical(daq.jump_detected)), ...
         daq.py(logical(daq.jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(daq.px(1), daq.py(1), 'r.', 'MarkerSize', 20)  % Start point


title("Forward Velocity colored")
save_plot_with_title_as_filename('x_fwd', 'y_fwd', savepath);

%% flat path colored by rot vel
cmap=jet(100);
speed_norm = round((cleaned_rot_speed-min(cleaned_rot_speed))/(max(cleaned_rot_speed)-min(cleaned_rot_speed))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(cleaned_xpos)-1
    plot(cleaned_xpos(i:i+1), cleaned_ypos(i:i+1), 'Color', cmap(speed_norm(i), :), 'LineWidth', 2)
end

% Create colorbar with actual rotation speed values
c = colorbar;
clim([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual rotation speed values
actual_values = yticks/100 * (max(cleaned_rot_speed)-min(cleaned_rot_speed)) + min(cleaned_rot_speed);
% Format tick labels to 1 decimal place
ylabel(c, 'Rotation Speed (deg/s)')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), actual_values, 'UniformOutput', false))

hold on
%Add jump points
if jump
    plot(daq.px(logical(daq.jump_detected)), ...
         daq.py(logical(daq.jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(daq.px(1), daq.py(1), 'r.', 'MarkerSize', 20)  % Start point
title("Rotation Speed colored")

save_plot_with_title_as_filename('x_rot', 'y_rot', savepath);


end

