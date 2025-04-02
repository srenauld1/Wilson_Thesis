function plot_raster(exptData,saccade_both, savepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% determine whether stim is 5s or 0.5s
if sum(exptData.optoStim) == 319968
    num_trials = 48;
    trial_duration = 15;
elseif sum(exptData.optoStim)== 31968
    num_trials = 48;
    trial_duration = 15;
else
    num_trials = 48;
    trial_duration = 15;

%% downsample data for plotting
newSampRate = 60; 
currentSampRate =length(exptData.t)/max(exptData.t);

% Assume exptData is defined and contains fields with time series data
downsampling_factor = round(currentSampRate/ newSampRate);  % Factor to downsample

% Initialize a new structure to hold the downsampled data
downsampled_exptData = struct();

% Loop through each field in exptData and downsample
fieldNames = fieldnames(exptData);  % Get all field names in the structure

for i = 1:length(fieldNames)
    field = fieldNames{i};
    
    % Check if the field contains numeric data and downsample if it does
    if isnumeric(exptData.(field))  % Check if the field is numeric
        % Downsample the data
        downsampled_exptData.(field) = exptData.(field)(1:downsampling_factor:end);  % Keep every nth sample
    else
        % If it's not numeric, copy the field directly (or handle differently if needed)
        downsampled_exptData.(field) = exptData.(field);
    end
end

true_sample_rate = length(downsampled_exptData.t)/max(downsampled_exptData.t);

%% now i need to separate by trial
% Create a figure
figure;
hold on;

% Loop through trials
for trial = 1:num_trials
    % Calculate the start and end time for this trial
    trial_start_idx = round((trial - 1) * trial_duration);
    trial_end_idx = round(trial * trial_duration);

    % Determine indices for this trial
    trial_indices = find(downsampled_exptData.t >= trial_start_idx & downsampled_exptData.t < trial_end_idx);

    % Plot opto on periods
    opto_on_times = downsampled_exptData.t(downsampled_exptData.optoStim(trial_indices) ~= 0);
    for i = 1:length(opto_on_times)
        plot([opto_on_times(i) opto_on_times(i)], [trial - 0.5, trial + 0.5], 'g', 'LineWidth', 2);  % Green for opto on
    end
    
    % Plot opto off periods
    opto_off_times = downsampled_exptData.t(downsampled_exptData.optoStim(trial_indices) == 0);
    for i = 1:length(opto_off_times)
        plot([opto_off_times(i) opto_off_times(i)], [trial - 0.5, trial + 0.5], 'r', 'LineWidth', 2);  % Red for opto off
    end

    % Mark saccades
    saccade_indices = find(downsampled_exptData.saccading(trial_indices) == 1);
    for s = 1:length(saccade_indices)
        % Get the absolute index in downsampled_exptData corresponding to this trial
        saccade_index = trial_indices(saccade_indices(s));  % This gives the correct index
        
        % Get time of the saccade
        saccade_time = downsampled_exptData.t(saccade_index);
        
        % Create a rectangle for the saccade, normalized within trial plot
        rectangle('Position', [saccade_time - trial_start_idx - 0.1, trial - 0.4, 0.2, 0.8], ...
                  'FaceColor', 'b', 'EdgeColor', 'none');  % Blue rectangle for saccade
    end
end

% Final adjustments to the plot
xlabel('Time (s)');
ylabel('Trial Number');
title('Raster Plot with Opto Stim and Saccades');
grid on;
set(gca, 'YDir', 'reverse');  % Flip y-axis so trial 1 is at the top
hold off;

save_plot_with_title_as_filename('raster', 'saccade', savepath)

end