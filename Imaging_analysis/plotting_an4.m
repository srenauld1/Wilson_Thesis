
%% function
function plotting_an4(a2p_data, jump, savepath)

close all


%[a2p_data.dq(1), triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(a2p_data.dq(1), savepath, 10, 2,0.88,0.88);
%% smooth fictrac
fwd = a2p_data.dq(1).bvf;
yaw = a2p_data.dq(1).bvy_deg;
side = a2p_data.dq(1).bvs_deg;
dff = a2p_data.roi.ts(1, :);
time= a2p_data.dq(1).t;


%% all 3 one plot
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
plot(time, yaw, 'r', 'LineWidth', 2);  % Plot yawational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (rad/s)');
title('Rotational Velocity');
grid on;

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_yaw', 'dff_separate', savepath);


%% flat path
figure
plot(a2p_data.dq(1).px, a2p_data.dq(1).py)
hold on
plot(a2p_data.dq(1).px(1), a2p_data.dq(1).py(1), 'r.', 'MarkerSize', 20)  % The MarkerSize value controls dot size
save_plot_with_title_as_filename('x', 'y', savepath);


%% flat path colored by dff
cmap=jet(100);
dff_norm = round((dff-min(dff))/(max(dff)-min(dff))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(a2p_data.dq(1).px)-1
    plot(a2p_data.dq(1).px(i:i+1), a2p_data.dq(1).py(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
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
    plot(a2p_data.dq(1).px(logical(a2p_data.dq(1).jump_detected)), ...
         a2p_data.dq(1).py(logical(a2p_data.dq(1).jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(a2p_data.dq(1).px(1), a2p_data.dq(1).py(1), 'r.', 'MarkerSize', 20)  % Start point

title("DFF colored, cue flat path")
save_plot_with_title_as_filename('x_color', 'y_color', savepath);
%% flat path colored by dff - ball
% figure
% hold on
% for i = 1:length(a2p_data.dq(1).px)-1
%     plot(a2p_data.dq(1).px_ball(i:i+1), a2p_data.dq(1).py_ball(i:i+1), 'Color', cmap(dff_norm(i), :), 'LineWidth', 2)
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
% plot(a2p_data.dq(1).px_ball(1), a2p_data.dq(1).py_ball(1), 'r.', 'MarkerSize', 20)  % Start point
% 
% title("DFF colored, ball flat path")
% save_plot_with_title_as_filename('x_ball_color', 'y_ball_color', savepath);

%% revamp variables

total_mov_mm = a2p_data.dq(1).totalspeed;    
not_moving = a2p_data.dq(1).motion.moving_or_not;
speed = abs(yaw);
% Identify outliers using the default method (usually interquartile range)
outliers1 = isoutlier(fwd, 'mean', 'ThresholdFactor', 6);
outliers2 = isoutlier(speed, 'mean', 'ThresholdFactor', 6);
outliers = outliers1 | outliers2;

% Remove outliers from forward_velocity
cleaned_forward_velocity = fwd(~outliers);
cleaned_yaw_speed = speed(~outliers);
cleaned_dff = dff(~outliers);
cleaned_xpos = a2p_data.dq(1).px(~outliers);
cleaned_ypos = a2p_data.dq(1).py(~outliers);

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
    plot(a2p_data.dq(1).px(logical(a2p_data.dq(1).jump_detected)), ...
         a2p_data.dq(1).py(logical(a2p_data.dq(1).jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(a2p_data.dq(1).px(1), a2p_data.dq(1).py(1), 'r.', 'MarkerSize', 20)  % Start point


title("Forward Velocity colored")
save_plot_with_title_as_filename('x_fwd', 'y_fwd', savepath);

%% flat path colored by yaw vel
cmap=jet(100);
speed_norm = round((cleaned_yaw_speed-min(cleaned_yaw_speed))/(max(cleaned_yaw_speed)-min(cleaned_yaw_speed))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(cleaned_xpos)-1
    plot(cleaned_xpos(i:i+1), cleaned_ypos(i:i+1), 'Color', cmap(speed_norm(i), :), 'LineWidth', 2)
end

% Create colorbar with actual yawation speed values
c = colorbar;
clim([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions
% Convert normalized values back to actual yawation speed values
actual_values = yticks/100 * (max(cleaned_yaw_speed)-min(cleaned_yaw_speed)) + min(cleaned_yaw_speed);
% Format tick labels to 1 decimal place
ylabel(c, 'Rotation Speed (deg/s)')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), actual_values, 'UniformOutput', false))

hold on
%Add jump points
if jump
    plot(a2p_data.dq(1).px(logical(a2p_data.dq(1).jump_detected)), ...
         a2p_data.dq(1).py(logical(a2p_data.dq(1).jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end
plot(a2p_data.dq(1).px(1), a2p_data.dq(1).py(1), 'r.', 'MarkerSize', 20)  % Start point
title("Rotation Speed colored")

save_plot_with_title_as_filename('x_yaw', 'y_yaw', savepath);


end

