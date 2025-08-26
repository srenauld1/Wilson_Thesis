
%% function
function plotting_ves041(daq, ts, jump, savepath)

close all


[daq, triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(daq, savepath, 10, 5,0.88,0.88);
%% smooth fictrac
fwd = daq.bfv_supp;
yaw = daq.byv_deg_supp;
side = daq.bsv_deg_supp;
dff = ts{1};
time=daq.t;
time_kin = daq.t_supp;


%% Plot the dff and fwd
% Plot the dff and fwd
figure;


yyaxis left;
plot(time_kin, fwd, '-b', 'LineWidth', 1.5)% Plot first time series in blue
ylabel('Forward Velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 10]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-r' , 'LineWidth', 1.5);  % Plot second time series in red
ylabel('dff');  % Label for the right y-axis

% Add title and x-axis label
xlabel('Time');
title('Forward velocity and DFF');

% Optional: Add legend
legend('Forward Velocity', 'dff');
save_plot_with_title_as_filename('Fwd vel', 'dff', savepath);

%% plot dff and yaw
figure;
yyaxis left;
plot(time_kin, yaw, '-g', 'LineWidth', 1.5);  % Plot first time series in blue
ylabel('Yaw velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-500 500]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time, dff, '-r', 'LineWidth', 1.5);  % Plot second time series in red
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
plot(time_kin, fwd, '-b');  % Plot first time series in blue
ylabel('Forward velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-5 5]);

% Plot the second time series on the right y-axis
yyaxis right;
plot(time_kin, yaw, '-g');  % Plot second time series in red
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
plot(time_kin, side, '-y');  % Plot first time series in blue
ylabel('Sideways velocity');  % Label for the left y-axis

% Set the y-axis limits for the left axis to [-5, 10]
ylim([-500 500]);

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
plot(time, dff, 'r', 'LineWidth', 1.5);  % Plot dFF in red
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F');
grid on;

% Subplot 2: Forward Velocity
subplot(3, 1, 2);  % Second subplot
plot(time_kin, fwd, 'black', 'LineWidth', 1.5);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([-5 20]);

% Subplot 3: Rotational Velocity
subplot(3, 1, 3);  % Third subplot
plot(time_kin, rot, 'b', 'LineWidth', 1.5);  % Plot rotational velocity in blue
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on;
ylim([-500 500]);

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall labels
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title

save_plot_with_title_as_filename('fwd_rot', 'dff_separate', savepath);

%% flat path colored by dff (dots with non-linear scaling)
cmap = jet(100);

% Non-linear scaling options for dff (choose one):
% Option 1: Square root scaling (compresses high values)
dff_scaled = sqrt((dff - min(dff)) / (max(dff) - min(dff)));

% Option 2: Logarithmic scaling (for positive values)
% dff_scaled = log((dff - min(dff) + 1)) / log(max(dff) - min(dff) + 1);

% Option 3: Power scaling (enhances small differences)
% dff_scaled = ((dff - min(dff)) / (max(dff) - min(dff))).^0.5;

% Option 4: Sigmoid scaling (S-curve)
% dff_norm_temp = (dff - min(dff)) / (max(dff) - min(dff));
% dff_scaled = 1 ./ (1 + exp(-10 * (dff_norm_temp - 0.5)));

% Convert scaled values to colormap indices
dff_norm = round(dff_scaled * (length(cmap) - 1)) + 1;

figure
hold on

% Plot as scatter points instead of lines
scatter(daq.px, daq.py, 10, cmap(dff_norm, :), 'filled', 'MarkerEdgeColor', 'none');

% Alternative: Plot as individual dots with plot()
% for i = 1:length(daq.px)
%     plot(daq.px(i), daq.py(i), '.', 'Color', cmap(dff_norm(i), :), 'MarkerSize', 15)
% end

% Create colorbar with actual DFF values
c = colorbar;
caxis([1 100])  % Set colorbar limits to match normalized range
yticks = get(c, 'YTick');  % Get current tick positions

% Convert normalized values back to actual DFF values using inverse transformation
% For square root scaling:
actual_scaled = (yticks - 1) / 99;  % Convert back to 0-1 range
actual_values = (actual_scaled.^2) * (max(dff) - min(dff)) + min(dff);

% For logarithmic scaling (uncomment if using log scaling):
% actual_scaled = (yticks - 1) / 99;
% actual_values = (exp(actual_scaled * log(max(dff) - min(dff) + 1)) - 1) + min(dff);

ylabel(c, 'DFF')
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))

% Add jump points
if jump
    plot(daq.px(logical(daq.jump_detected)), ...
         daq.py(logical(daq.jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
end

% Start point
plot(daq.px(1), daq.py(1), 'r.', 'MarkerSize', 20)

title("DFF colored (non-linear), cue flat path - dots")
save_plot_with_title_as_filename('x_color', 'y_color', savepath);
%% now interactively plot dff and yaw velocity
trajectory_region_selector(daq, ts{1}, jump)

%% can add in specifics with the plotting
end


