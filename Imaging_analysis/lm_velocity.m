function lm_velocity(daq, dat, savepath)

%% Define variables and downsample kinematics to match dff
dff = dat.ts(1,:);
time_dff = daq.t;  % dff time base
time_kinematics = daq.t_supp;  % kinematic time base

% Original kinematic variables (high resolution, unaltered)
forward_velocity_hr = daq.bvf_supp;  
rot_velocity_hr = daq.bvy_deg_supp;
side_velocity_hr = daq.bvs_deg_supp;

% Downsample kinematic variables to match dff time base
forward_velocity = interp1(time_kinematics, forward_velocity_hr, time_dff, 'linear', 'extrap');
rot_velocity = interp1(time_kinematics, rot_velocity_hr, time_dff, 'linear', 'extrap');
side_velocity = interp1(time_kinematics, side_velocity_hr, time_dff, 'linear', 'extrap');

% Calculate derived variables
rot_speed = abs(rot_velocity);
total_speed = abs(forward_velocity) + abs(rot_velocity) + abs(side_velocity);

fprintf('Original data points - DFF: %d, Kinematics: %d\n', length(dff), length(forward_velocity_hr));
fprintf('After downsampling - DFF: %d, Kinematics: %d\n', length(dff), length(forward_velocity));

%% Remove start/stop transitions
ex_startstop = 1;
postStartWin = 0.5; % Time window after start (in seconds)
preStopWin = 0.5;   % Time window before stop (in seconds)

if ex_startstop
    % Convert time windows to sample indices
    dt = mean(diff(time_dff)); % Average sampling interval
    postStartIdx = round(postStartWin / dt);
    preStopIdx = round(preStopWin / dt);
    
    % Calculate run index using Schmitt Trigger on total speed
    runIdx = schmittTrigger(total_speed, 1, 1);
    runIdx = runIdx(:);
    
    % Identify start and stop transitions
    runTransitions = diff([0; runIdx; 0]); % Pad to catch edge cases
    startTrans = find(runTransitions == 1); % Start running
    stopTrans = find(runTransitions == -1); % Stop running
    
    % Set post-start periods to NaN
    for st = 1:length(startTrans)
        tStart = startTrans(st);
        tEnd = min(length(dff), tStart + postStartIdx);
        dff(tStart:tEnd) = NaN;
        forward_velocity(tStart:tEnd) = NaN;
        rot_speed(tStart:tEnd) = NaN;
        side_velocity(tStart:tEnd) = NaN;
        total_speed(tStart:tEnd) = NaN;
    end
    
    % Set pre-stop periods to NaN
    for sp = 1:length(stopTrans)
        tStop = stopTrans(sp);
        tStart = max(1, tStop - preStopIdx);
        dff(tStart:tStop) = NaN;
        forward_velocity(tStart:tStop) = NaN;
        rot_speed(tStart:tStop) = NaN;
        side_velocity(tStart:tStop) = NaN;
        total_speed(tStart:tStop) = NaN;
    end
    
    % Set non-running periods to NaN
    dff(runIdx == 0) = NaN;
    forward_velocity(runIdx == 0) = NaN;
    rot_speed(runIdx == 0) = NaN;
    side_velocity(runIdx == 0) = NaN;
    total_speed(runIdx == 0) = NaN;
end

% Remove NaN values for modeling
valid_idx = ~isnan(dff) & ~isnan(forward_velocity) & ~isnan(rot_speed) & ~isnan(side_velocity);
model_dff = dff(valid_idx);
model_forward = forward_velocity(valid_idx);
model_rot = rot_speed(valid_idx);
model_side = side_velocity(valid_idx);
model_total_speed = total_speed(valid_idx);

fprintf('Valid data points for modeling: %d\n', sum(valid_idx));

%% Model 1: DFF vs Forward Velocity
lm_forward = fitlm(model_forward, model_dff);
r2_forward = lm_forward.Rsquared.Ordinary;

disp('=== DFF vs Forward Velocity ===');
disp(lm_forward);
disp(['R^2 = ', num2str(r2_forward)]);

figure;
plot(lm_forward);
xlabel('Forward Velocity (mm/s)');
ylabel('dF/F');
title(['DFF vs Forward Velocity: R^2 = ', num2str(r2_forward, '%.3f')]);
grid on;
save_plot_with_title_as_filename('dff_vs_forward_velocity', 'linear_model', savepath);

%% Model 2: DFF vs Rotational Velocity
lm_rot = fitlm(model_rot, model_dff);
r2_rot = lm_rot.Rsquared.Ordinary;

disp('=== DFF vs Rotational Velocity ===');
disp(lm_rot);
disp(['R^2 = ', num2str(r2_rot)]);

figure;
plot(lm_rot);
xlabel('Rotational Velocity (deg/s)');
ylabel('dF/F');
title(['DFF vs Rotational Velocity: R^2 = ', num2str(r2_rot, '%.3f')]);
grid on;
save_plot_with_title_as_filename('dff_vs_rotational_velocity', 'linear_model', savepath);

%% Model 3: DFF vs Side Velocity
lm_side = fitlm(model_side, model_dff);
r2_side = lm_side.Rsquared.Ordinary;

disp('=== DFF vs Side Velocity ===');
disp(lm_side);
disp(['R^2 = ', num2str(r2_side)]);

figure;
plot(lm_side);
xlabel('Side Velocity (mm/s)');
ylabel('dF/F');
title(['DFF vs Side Velocity: R^2 = ', num2str(r2_side, '%.3f')]);
grid on;
save_plot_with_title_as_filename('dff_vs_side_velocity', 'linear_model', savepath);

%% Model 4: DFF vs Total Speed
lm_total = fitlm(model_total_speed, model_dff);
r2_total = lm_total.Rsquared.Ordinary;

disp('=== DFF vs Total Speed ===');
disp(lm_total);
disp(['R^2 = ', num2str(r2_total)]);

figure;
plot(lm_total);
xlabel('Total Speed (mm/s + deg/s)');
ylabel('dF/F');
title(['DFF vs Total Speed: R^2 = ', num2str(r2_total, '%.3f')]);
grid on;
save_plot_with_title_as_filename('dff_vs_total_speed', 'linear_model', savepath);

%% Combined Model: DFF vs All Velocities
data_table = table(model_forward, model_rot, model_side, model_total_speed, model_dff, ...
    'VariableNames', {'ForwardVelocity', 'RotationalVelocity', 'SideVelocity', 'TotalSpeed', 'DFF'});

lm_combined = fitlm(data_table, 'DFF ~ ForwardVelocity + RotationalVelocity + SideVelocity');
r2_combined = lm_combined.Rsquared.Ordinary;

disp('=== Combined Model: DFF vs All Velocities ===');
disp(lm_combined);
disp(['Combined R^2 = ', num2str(r2_combined)]);

% Plot individual relationships
figure;
subplot(2, 2, 1);
scatter(model_forward, model_dff, 20, 'k.', 'MarkerEdgeAlpha', 0.6);
xlabel('Forward Velocity (mm/s)');
ylabel('dF/F');
title(['Forward: R^2 = ', num2str(r2_forward, '%.3f')]);
grid on;

subplot(2, 2, 2);
scatter(model_rot, model_dff, 20, 'r.', 'MarkerEdgeAlpha', 0.6);
xlabel('Rotational Velocity (deg/s)');
ylabel('dF/F');
title(['Rotational: R^2 = ', num2str(r2_rot, '%.3f')]);
grid on;

subplot(2, 2, 3);
scatter(model_side, model_dff, 20, 'g.', 'MarkerEdgeAlpha', 0.6);
xlabel('Side Velocity (mm/s)');
ylabel('dF/F');
title(['Side: R^2 = ', num2str(r2_side, '%.3f')]);
grid on;

subplot(2, 2, 4);
scatter(model_total_speed, model_dff, 20, 'm.', 'MarkerEdgeAlpha', 0.6);
xlabel('Total Speed');
ylabel('dF/F');
title(['Total Speed: R^2 = ', num2str(r2_total, '%.3f')]);
grid on;

sgtitle(['All Velocity Models (Combined R^2 = ', num2str(r2_combined, '%.3f'), ')']);
save_plot_with_title_as_filename('dff_vs_all_velocities', 'summary', savepath);

%% Cross-correlations
max_lag = 25; % Adjust based on your data

% Forward velocity
[cc_forward, lags] = xcorr(model_forward, model_dff, max_lag, 'coeff');
[max_cc_fwd, idx_fwd] = max(cc_forward);
max_lag_fwd = lags(idx_fwd);

% Rotational velocity  
[cc_rot, lags] = xcorr(model_rot, model_dff, max_lag, 'coeff');
[max_cc_rot, idx_rot] = max(cc_rot);
max_lag_rot = lags(idx_rot);

% Side velocity
[cc_side, lags] = xcorr(model_side, model_dff, max_lag, 'coeff');
[max_cc_side, idx_side] = max(cc_side);
max_lag_side = lags(idx_side);

% Total speed
[cc_total, lags] = xcorr(model_total_speed, model_dff, max_lag, 'coeff');
[max_cc_total, idx_total] = max(cc_total);
max_lag_total = lags(idx_total);

% Display results
fprintf('\n=== Cross-Correlation Results ===\n');
fprintf('Forward velocity - Max CC: %.3f at lag %d\n', max_cc_fwd, max_lag_fwd);
fprintf('Rotational velocity - Max CC: %.3f at lag %d\n', max_cc_rot, max_lag_rot);
fprintf('Side velocity - Max CC: %.3f at lag %d\n', max_cc_side, max_lag_side);
fprintf('Total speed - Max CC: %.3f at lag %d\n', max_cc_total, max_lag_total);

% Plot cross-correlations
figure;
subplot(2, 2, 1);
plot(lags, cc_forward, 'k', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title(['Forward Vel: Max = ', num2str(max_cc_fwd, '%.3f'), ' at lag ', num2str(max_lag_fwd)]);
grid on;

subplot(2, 2, 2);
plot(lags, cc_rot, 'r', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title(['Rotational Vel: Max = ', num2str(max_cc_rot, '%.3f'), ' at lag ', num2str(max_lag_rot)]);
grid on;

subplot(2, 2, 3);
plot(lags, cc_side, 'g', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title(['Side Vel: Max = ', num2str(max_cc_side, '%.3f'), ' at lag ', num2str(max_lag_side)]);
grid on;

subplot(2, 2, 4);
plot(lags, cc_total, 'm', 'LineWidth', 1.5);
xlabel('Lag (samples)');
ylabel('Cross-correlation');
title(['Total Speed: Max = ', num2str(max_cc_total, '%.3f'), ' at lag ', num2str(max_lag_total)]);
grid on;

sgtitle('Cross-Correlations: Velocities vs DFF');
save_plot_with_title_as_filename('cross_correlations', 'velocities_dff', savepath);

%% Summary Statistics
fprintf('\n=== Summary of R^2 Values ===\n');
fprintf('Forward velocity: %.3f\n', r2_forward);
fprintf('Rotational velocity: %.3f\n', r2_rot);
fprintf('Side velocity: %.3f\n', r2_side);
fprintf('Total speed: %.3f\n', r2_total);
fprintf('Combined model: %.3f\n', r2_combined);

% Save results to workspace
results.models.forward = lm_forward;
results.models.rotational = lm_rot;
results.models.side = lm_side;
results.models.total_speed = lm_total;
results.models.combined = lm_combined;

results.r2.forward = r2_forward;
results.r2.rotational = r2_rot;
results.r2.side = r2_side;
results.r2.total_speed = r2_total;
results.r2.combined = r2_combined;

results.cross_correlation.forward = struct('max_cc', max_cc_fwd, 'lag', max_lag_fwd);
results.cross_correlation.rotational = struct('max_cc', max_cc_rot, 'lag', max_lag_rot);
results.cross_correlation.side = struct('max_cc', max_cc_side, 'lag', max_lag_side);
results.cross_correlation.total_speed = struct('max_cc', max_cc_total, 'lag', max_lag_total);

assignin('base', 'velocity_dff_analysis_results', results);
fprintf('\nResults saved to workspace as "velocity_dff_analysis_results"\n');

end


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
% % %% scatterplot colored by rv
% % Example data (replace with your actual data)
% forward_velocity = cleaned_forward_velocity;  % Forward velocity data (x-axis)
% dff = cleaned_dff;  % dff data (y-axis)
% rotational_velocity = cleaned_rot_vel;  % Rotational velocity data (color coding)
% 
% 
% 
% % Define the bin edges for rotational velocity
% rotational_bins = [-3 -2 -1 0 1 2 3];
% 
% % Bin the rotational velocity data
% [~, rotational_bin_idx] = histc(rotational_velocity, rotational_bins);
% 
% % Initialize a cell array to store dff values and forward velocity for each rotational velocity bin
% binned_forward_velocity = cell(length(rotational_bins)-1, 1);
% binned_dff = cell(length(rotational_bins)-1, 1);
% 
% % Separate forward velocity and dff into bins based on rotational velocity
% for j = 1:length(rotational_bins)-1
%     % Find the indices where rotational velocity falls into the current bin
%     bin_indices = rotational_bin_idx == j;
% 
%     % Store forward velocity and dff for the current bin
%     binned_forward_velocity{j} = forward_velocity(bin_indices);
%     binned_dff{j} = dff(bin_indices);
% end
% 
% % Plot the points for each rotational velocity bin using scatter
% figure;
% hold on;
% colors = lines(length(rotational_bins)-1);  % Get distinct colors for each bin
% for j = 1:length(rotational_bins)-1
%     scatter(binned_forward_velocity{j}, binned_dff{j}, 50, 'MarkerEdgeColor', colors(j, :), 'MarkerFaceColor', colors(j, :));
% end
% 
% % Plot average lines for each bin
% for i = 1:length(rotational_bins)-1
%     % Find the indices of points within the current rotational velocity bin
%     bin_indices = find(rotational_velocity >= rotational_bins(i) & rotational_velocity < rotational_bins(i+1));
% 
%     % Check if there are enough points in the bin
%     if ~isempty(bin_indices)
%         % Get the forward velocity and dff values for the current bin
%         fwd_vel_bin = forward_velocity(bin_indices);
%         dff_bin = dff(bin_indices);
% 
%         % Sort the forward velocity and dff for better plotting
%         [fwd_vel_sorted, sort_idx] = sort(fwd_vel_bin);
%         dff_sorted = dff_bin(sort_idx);
% 
%         % Compute a moving average (optional) to smooth the curve, or just plot the raw sorted data
%         window_size = 5;  % You can adjust this size
%         smoothed_dff = movmean(dff_sorted, window_size);  % Smoothing dff
% 
%         % Plot the average relationship for this bin
%         plot(fwd_vel_sorted, smoothed_dff, 'LineWidth', 2, 'Color', colors(i, :), ...
%             'DisplayName', sprintf('RotVel: [%d, %d] deg/s', rotational_bins(i), rotational_bins(i+1)));
%     end
% end
% 
% % Customize plot
% xlabel('Forward Velocity (mm/s)');
% ylabel('dff');
% title('dff vs Forward Velocity for Different Binned Rotational Velocities');
% legend(arrayfun(@(x,y) sprintf('%.1f to %.1f mm/s', x, y), rotational_bins(1:end-1), rotational_bins(2:end), 'UniformOutput', false));
% grid on;
% hold off;
% 
% save_plot_with_title_as_filename('scatter', 'average', savepath);
% 
% %% line colored by rv, binned forward and rv
% 
% % Define the bin edges for forward and rotational velocity
% forward_bins = [-3 -1 0 1 3 5 7 9]; %[-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9 9.5, 10];%
% rotational_bins = [-3 -2 -1 0 1 2 3];
% tabulate(rotational_bin_idx)
% 
% % Bin the forward velocity data
% [~, forward_bin_idx] = histc(forward_velocity, forward_bins);
% 
% % Bin the rotational velocity data
% [~, rotational_bin_idx] = histc(rotational_velocity, rotational_bins);
% 
% % Initialize arrays to store binned dff values
% binned_dff = NaN(length(forward_bins)-1, length(rotational_bins)-1);
% 
% % Calculate the mean dff for each combination of forward and rotational bins
% for i = 1:length(forward_bins)-1
%     for j = 1:length(rotational_bins)-1
%         % Find the indices where both forward and rotational velocities fall into the current bins
%         bin_indices = forward_bin_idx == i & rotational_bin_idx == j;
% 
%         % Compute the mean dff for the current bin
%         if any(bin_indices)
%             binned_dff(i, j) = mean(dff(bin_indices));
%         end
%     end
% end
% 
% % Plot the lines for each rotational bin
% figure;
% hold on;
% colors = lines(length(rotational_bins)-1);  % Get distinct colors for each bin
% for j = 1:length(rotational_bins)-1
%     plot(forward_bins(1:end-1) + diff(forward_bins)/2, binned_dff(:, j), 'LineWidth', 2, 'Color', colors(j, :));
% end
% 
% % Customize plot
% xlabel('Forward Velocity (binned)');
% ylabel('Mean dff');
% title('dff vs Binned Forward Velocity for Different Binned Rotational Velocities');
% legend(arrayfun(@(x,y) sprintf('%.1f to %.1f mm/s', x, y), rotational_bins(1:end-1), rotational_bins(2:end), 'UniformOutput', false));
% grid on;
% hold off;

%save_plot_with_title_as_filename('binned fwd', 'dff_rv', savepath);

%% heatmap

%dff_v_velocityHeat(cleaned_forward_velocity,cleaned_rot_vel,cleaned_side,cleaned_dff,1);
%heatmap_save = fullfile(savepath, "heatmap");
%savefig(gcf, heatmap_save)

