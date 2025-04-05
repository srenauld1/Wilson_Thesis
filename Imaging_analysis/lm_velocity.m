function lm_velocity(daq, ts, savepath)


close all

%% def variables
forward_velocity = daq.bfv; 
rot_velocity = daq.byv;
side_velocity = daq.bsv;
rot_speed = abs(rot_velocity);
dff = ts{1};
% Identify outliers using the default method (usually interquartile range)
outliers1 = isoutlier(forward_velocity, 'mean', 'ThresholdFactor', 4);
outliers2 = isoutlier(rot_velocity, 'mean', 'ThresholdFactor', 4);
outliers = outliers1 | outliers2;

% Remove outliers from forward_velocity
cleaned_forward_velocity = forward_velocity(~outliers);
cleaned_dff = dff(~outliers);
full_dff = cleaned_dff;
cleaned_rot = rot_speed(~outliers);
cleaned_side = side_velocity(~outliers);
cleaned_rot_vel = rot_velocity(~outliers);
ttime = daq.t(~outliers);
forward_speed = abs(cleaned_forward_velocity);

cleaned_forward_velocity = smooth(cleaned_forward_velocity);
cleaned_rot = smooth(cleaned_rot);
forward_speed = abs(cleaned_forward_velocity);


% Display the number of outliers removed
disp(['Number of outliers removed: ', num2str(sum(outliers))]);


%% Remove stop start transitions
% Set flag to exclude start/stop transitions using transition window settings
ex_startstop = 1;
postStartWin = 0.1; % Time window after start (in seconds)
preStopWin = 0.2;   % Time window before stop (in seconds)
% Number of trials (columns in the cell activity array)
nTrials = 1;
if ex_startstop
    % Convert post-start and pre-stop windows to indices based on time array
    postStartIdx = fetchTimeIdx(ttime, postStartWin);
    preStopIdx = fetchTimeIdx(ttime, preStopWin);
    % Loop over each trial
    for trial = 1:nTrials
        % Calculate run index using Schmitt Trigger
        runIdx = schmittTrigger(forward_speed(:, trial), 5, 5);
        %runIdx2 = schmittTrigger(cleaned_rot(:, trial), 0.1, 0.1);
        %runIdx = runIdx1 | runIdx2;
        % Identify start and stop transitions in runIdx for the current trial
        runTransitions = diff(runIdx);    % Calculate transitions in run state
        startTrans = find(runTransitions == 1); % 0 to 1 (start running)
        stopTrans = find(runTransitions == -1); % 1 to 0 (stop running)
        % Loop over each start transition to set post-start period as NaN
        for st = 1:length(startTrans)
            tStart = startTrans(st); % Start index
            tEnd = min(size(cleaned_dff, 2), tStart + postStartIdx); % End index, within bounds
            cleaned_dff(tStart:tEnd) = nan; % Set post-start window to NaN in cell activity data
        end
        % Loop over each stop transition to set pre-stop period as NaN
        for sp = 1:length(stopTrans)
            tStop = stopTrans(sp); % Stop index
            tStart = max(1, tStop - preStopIdx); % Start index, within bounds
            cleaned_dff(tStart:tStop) = nan; % Set pre-stop window to NaN in cell activity data
        end
        % Set cellactivity to NaN where runIdx is 0 (not running)
        cleaned_dff(runIdx == 0) = nan;
    end
end


%% scatterplot colored by rv
% Example data (replace with your actual data)
forward_velocity = cleaned_forward_velocity;  % Forward velocity data (x-axis)
dff = cleaned_dff;  % dff data (y-axis)
rotational_velocity = cleaned_rot_vel;  % Rotational velocity data (color coding)



% Define the bin edges for rotational velocity
rotational_bins = [-3 -2 -1 0 1 2 3];

% Bin the rotational velocity data
[~, rotational_bin_idx] = histc(rotational_velocity, rotational_bins);

% Initialize a cell array to store dff values and forward velocity for each rotational velocity bin
binned_forward_velocity = cell(length(rotational_bins)-1, 1);
binned_dff = cell(length(rotational_bins)-1, 1);

% Separate forward velocity and dff into bins based on rotational velocity
for j = 1:length(rotational_bins)-1
    % Find the indices where rotational velocity falls into the current bin
    bin_indices = rotational_bin_idx == j;
    
    % Store forward velocity and dff for the current bin
    binned_forward_velocity{j} = forward_velocity(bin_indices);
    binned_dff{j} = dff(bin_indices);
end

% Plot the points for each rotational velocity bin using scatter
figure;
hold on;
colors = lines(length(rotational_bins)-1);  % Get distinct colors for each bin
for j = 1:length(rotational_bins)-1
    scatter(binned_forward_velocity{j}, binned_dff{j}, 50, 'MarkerEdgeColor', colors(j, :), 'MarkerFaceColor', colors(j, :));
end

% Plot average lines for each bin
for i = 1:length(rotational_bins)-1
    % Find the indices of points within the current rotational velocity bin
    bin_indices = find(rotational_velocity >= rotational_bins(i) & rotational_velocity < rotational_bins(i+1));
    
    % Check if there are enough points in the bin
    if ~isempty(bin_indices)
        % Get the forward velocity and dff values for the current bin
        fwd_vel_bin = forward_velocity(bin_indices);
        dff_bin = dff(bin_indices);
        
        % Sort the forward velocity and dff for better plotting
        [fwd_vel_sorted, sort_idx] = sort(fwd_vel_bin);
        dff_sorted = dff_bin(sort_idx);
        
        % Compute a moving average (optional) to smooth the curve, or just plot the raw sorted data
        window_size = 5;  % You can adjust this size
        smoothed_dff = movmean(dff_sorted, window_size);  % Smoothing dff
        
        % Plot the average relationship for this bin
        plot(fwd_vel_sorted, smoothed_dff, 'LineWidth', 2, 'Color', colors(i, :), ...
            'DisplayName', sprintf('RotVel: [%d, %d] deg/s', rotational_bins(i), rotational_bins(i+1)));
    end
end

% Customize plot
xlabel('Forward Velocity (mm/s)');
ylabel('dff');
title('dff vs Forward Velocity for Different Binned Rotational Velocities');
legend(arrayfun(@(x,y) sprintf('%.1f to %.1f mm/s', x, y), rotational_bins(1:end-1), rotational_bins(2:end), 'UniformOutput', false));
grid on;
hold off;

save_plot_with_title_as_filename('scatter', 'average', savepath);

%% line colored by rv, binned forward and rv

% Define the bin edges for forward and rotational velocity
forward_bins = [-3 -1 0 1 3 5 7 9]; %[-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9 9.5, 10];%
rotational_bins = [-3 -2 -1 0 1 2 3];
tabulate(rotational_bin_idx)

% Bin the forward velocity data
[~, forward_bin_idx] = histc(forward_velocity, forward_bins);

% Bin the rotational velocity data
[~, rotational_bin_idx] = histc(rotational_velocity, rotational_bins);

% Initialize arrays to store binned dff values
binned_dff = NaN(length(forward_bins)-1, length(rotational_bins)-1);

% Calculate the mean dff for each combination of forward and rotational bins
for i = 1:length(forward_bins)-1
    for j = 1:length(rotational_bins)-1
        % Find the indices where both forward and rotational velocities fall into the current bins
        bin_indices = forward_bin_idx == i & rotational_bin_idx == j;
        
        % Compute the mean dff for the current bin
        if any(bin_indices)
            binned_dff(i, j) = mean(dff(bin_indices));
        end
    end
end

% Plot the lines for each rotational bin
figure;
hold on;
colors = lines(length(rotational_bins)-1);  % Get distinct colors for each bin
for j = 1:length(rotational_bins)-1
    plot(forward_bins(1:end-1) + diff(forward_bins)/2, binned_dff(:, j), 'LineWidth', 2, 'Color', colors(j, :));
end

% Customize plot
xlabel('Forward Velocity (binned)');
ylabel('Mean dff');
title('dff vs Binned Forward Velocity for Different Binned Rotational Velocities');
legend(arrayfun(@(x,y) sprintf('%.1f to %.1f mm/s', x, y), rotational_bins(1:end-1), rotational_bins(2:end), 'UniformOutput', false));
grid on;
hold off;

%save_plot_with_title_as_filename('binned fwd', 'dff_rv', savepath);

%% heatmap

dff_v_velocityHeat(cleaned_forward_velocity,cleaned_rot_vel,cleaned_side,cleaned_dff,1);
heatmap_save = fullfile(savepath, "heatmap");
%savefig(gcf, heatmap_save)



%% cleaned forward speed

% Create a linear model to predict dff from forward velocity
lm_fwd_speed= fitlm(forward_speed, cleaned_dff);

% Display the linear model summary
disp(lm_fwd_speed);

% Extract and display the R^2 value
r2 = lm_fwd_speed.Rsquared.Ordinary;  % Ordinary R^2
disp(['R^2 for predicting dff from forward speed: ', num2str(r2)]);
figure;
plot(lm_fwd_speed);
xlabel('Forward Speed');
ylabel('dF/F');
title(['Linear Fit: R^2 = ', num2str(r2)]);
grid on;

save_plot_with_title_as_filename('r2_Fwd vel', 'dff', savepath);



%% rot 
% Create a linear model to predict dff from rotaitonal velocity
lm_rot = fitlm(cleaned_rot, cleaned_dff);

% Display the linear model summary
disp(lm_rot);

% Extract and display the R^2 value
r2 = lm_rot.Rsquared.Ordinary;  % Ordinary R^2
disp(['R^2 for predicting dff from rot velocity: ', num2str(r2)]);

% Plot the linear fit and data points
figure;
plot(lm_rot);
xlabel('Rotational Velocity');
ylabel('dF/F');
title(['Linear Fit: R^2 = ', num2str(r2)]);
grid on;

save_plot_with_title_as_filename('r2_rot vel', 'dff', savepath);

%% both

% Combine predictors into a table
data_table = table(forward_speed, cleaned_rot, cleaned_dff', ...
    'VariableNames', {'ForwardVelocity', 'RotationalVelocity', 'dFF'});

% Create a linear model to predict dff from forward and rotational velocity
lm = fitlm(data_table);

% Display the linear model summary
lm

% Extract and display the R^2 value
r2 = lm.Rsquared.Ordinary;  % Ordinary R^2
disp(['R^2 for predicting dff: ', num2str(r2)]);

% Plot the model fit for each predictor
figure;
subplot(2, 1, 1);
scatter(cleaned_forward_velocity, cleaned_dff, 'b.');
xlabel('Forward Velocity');
ylabel('dF/F');
title('dFF vs Forward Velocity');

subplot(2, 1, 2);
scatter(cleaned_rot, cleaned_dff, 'r.')
xlabel('Rotational Velocity');
ylabel('dF/F');
title('dFF vs Rotational Velocity');

% Optional: Plot predicted vs. actual dFF
% figure;
% predicted_dff = predict(lm, predictors);
% scatter(cleaned_dff, predicted_dff, 'k.');
% xlabel('Actual dF/F');
% ylabel('Predicted dF/F');
% title(['Predicted vs. Actual dF/F: R^2 = ', num2str(r2)]);
% grid on;

%% cross correlations
% Specify the maximum lag
max_lag = 40;

% Compute cross-correlation for forward velocity and dFF
[cc_forward, lags_forward] = xcorr(forward_speed, full_dff, max_lag, 'coeff');
[max_cc_fwd, idx_fwd] = max(cc_forward);
max_cc_lag_fwd = lags_forward(idx_fwd);

% Compute cross-correlation for rotational velocity and dFF
[cc_rotational, lags_rotational] = xcorr(cleaned_rot, full_dff, max_lag, 'coeff');
[max_cc_rot, idx_rot] = max(cc_rotational);
max_cc_lag_rot = lags_rotational(idx_rot);

% Display the results
disp(['Maximum cross-correlation (forward speed): ', num2str(max_cc_fwd)]);
disp(['Lag at max cross-correlation (forward speed): ', num2str(max_cc_lag_fwd)]);
disp(['Maximum cross-correlation (rotational): ', num2str(max_cc_rot)]);
disp(['Lag at max cross-correlation (rotational): ', num2str(max_cc_lag_rot)]);

% Plot cross-correlations
figure;

% Forward velocity vs. dFF
subplot(2, 1, 1);
plot(lags_forward, cc_forward, 'b', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title('Cross-Correlation: Forward Velocity vs. dFF');
grid on;

% Rotational velocity vs. dFF
subplot(2, 1, 2);
plot(lags_rotational, cc_rotational, 'g', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title('Cross-Correlation: Rotational Velocity vs. dFF');
grid on;

save_plot_with_title_as_filename('crosscorr', 'velocities', savepath);



% % Shift the data based on the lag
% if max_cc_lag_fwd > 0
%     shifted_forward_velocity = forward_speed(1+max_cc_lag_fwd:end);
%     shifted_dff_fwd = cleaned_dff(1:end-max_cc_lag_fwd);
% elseif max_cc_lag_fwd < 0
%     shifted_forward_velocity = forward_speed(1:end+max_cc_lag_fwd);
%     shifted_dff_fwd = cleaned_dff(1-max_cc_lag_fwd:end);
% else
%     shifted_forward_velocity = forward_speed;
%     shifted_dff_fwd = cleaned_dff;
% end
% 
% if max_cc_lag_rot > 0
%     shifted_rotational_velocity = cleaned_rot(1+max_cc_lag_rot:end);
%     shifted_dff_rot = cleaned_dff(1:end-max_cc_lag_rot);
% elseif max_cc_lag_rot < 0
%     shifted_rotational_velocity = cleaned_rot(1:end+max_cc_lag_rot);
%     shifted_dff_rot = cleaned_dff(1-max_cc_lag_rot:end);
% else
%     shifted_rotational_velocity = cleaned_rot;
%     shifted_dff_rot = cleaned_dff;
% end
% 
% shifted_rotational_velocity_fwd = cleaned_rot(1:end+max_cc_lag_fwd);
% 
% 
% 
% % Fit linear models at the lag with maximum cross-correlation
% lm_fwd = fitlm(shifted_forward_velocity, shifted_dff_fwd);
% lm_rot = fitlm(shifted_rotational_velocity, shifted_dff_rot);
% lm_rot_fwd_shift =  fitlm(shifted_rotational_velocity_fwd, shifted_dff_fwd);
% 
% % Display the model summaries
% disp('Linear model for Forward Velocity:');
% disp(lm_fwd);
% 
% disp('Linear model for Rotational Velocity:');
% disp(lm_rot);
% 
% disp('Linear model for Both Velocity:');
% disp(lm_rot_fwd_shift);
% 
% 
% 
% 
