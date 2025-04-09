function ts = plotting_neck(daq, ts, savepath)

%% all 3 on one plot
time = daq.t;                  % Time points
dff = ts;           % Example dFF data
fwd = daq.bfv; % Example forward velocity data
rot = daq.byv; % Example rotational velocity data
side = daq.bsv;

% Subplot 1: dFF
numRows = size(dff, 1);
% Initialize a matrix to hold the normalized values
normalized_data = zeros(size(dff));

% Normalize each row
for i = 1:size(dff, 1)  % Iterate through each row
    row_min = min(dff(i, :));  % Minimum value of the current row
    row_max = max(dff(i, :));  % Maximum value of the current row

    % Normalize and assign to the normalized_data matrix
    normalized_data(i, :) = (dff(i, :) - row_min) / (row_max - row_min);
end
figure;

% Subplot 1: Normalized dF/F
subplot(4, 1, 1);  % First subplot
hold on;
for i = 1:numRows  % Assuming you want to plot all rows
    plot(time, normalized_data(i,:), 'LineWidth', 1);  % Plot normalized dFF
end
xlabel('Time (s)');
ylabel('Normalized dF/F');
title('Normalized dF/F');
grid on;
legend

% Subplot 2: Forward Velocity
subplot(4, 1, 2);  % Second subplot
plot(time, fwd, 'g', 'LineWidth', 2);  % Plot forward velocity in green
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;

% Subplot 3: Sideways Velocity
subplot(4, 1, 3);  % Third subplot
plot(time, side, 'b', 'LineWidth', 2);  % Plot sideways velocity in blue
xlabel('Time (s)');
ylabel('Side Velocity (mm/s)');
title('Sideways Velocity');
grid on;

% Subplot 4: Rotational Velocity
subplot(4, 1, 4);  % Fourth subplot
plot(time, rot, 'r', 'LineWidth', 2);  % Plot rotational velocity
xlabel('Time (s)');
ylabel('Rotational Velocity (mm/s)');
title('Rotational Velocity');
grid on;

% Link the x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

% Add overall title for the figure
sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title
end