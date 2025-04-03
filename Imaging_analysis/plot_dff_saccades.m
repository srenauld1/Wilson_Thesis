function plot_dff_saccades(ts, saccade_both);
    % merge left and right
    % Predefined parameters
    pre_saccade_time = 2;  % Time window (seconds) before the saccade
    post_saccade_time = 2;  % Time window (seconds) after the saccade
    
    figure;
    hold on;
    
    % Loop over each bout to extract dF/F traces
    for i = 1:length(saccade_both.boutStartInd)
        % Current start and end times of the saccade
        start_time = saccade_both.boutStartTimes(i);
        end_time = saccade_both.boutEndTimes(i);
        
        % Get the corresponding indices in the time vector for the full window
        start_idx = find(ts.t >= (start_time - pre_saccade_time), 1, 'first');  % Start index
        end_idx = find(ts.t >= (end_time + post_saccade_time), 1, 'first');      % End index
    
        % Ensure that the indices are valid and within bounds
        if ~isempty(start_idx) && ~isempty(end_idx) && end_idx > start_idx
            % Extract the dF/F trace and corresponding time window
            dff_trace = ts.resp.i2{1}(start_idx:end_idx);
    
            % Normalize to range [0, 1]
            dff_min = min(dff_trace);  % Calculate the minimum value
            dff_max = max(dff_trace);  % Calculate the maximum value
            
            % Prevent division by zero if all values are the same
            if dff_max - dff_min == 0
                normalized_dff_trace = zeros(size(dff_trace));  % If all values are the same, return zeros
            else
                normalized_dff_trace = (dff_trace - dff_min) / (dff_max - dff_min);  % Apply normalization
            end
            
    
            all_traces{i} = normalized_dff_trace; % Extract dF/F trace
            time_window = ts.t(start_idx:end_idx) - start_time;    % Time window relative to saccade start
    
        
            % Plot the aligned dF/F trace
            plot(time_window, normalized_dff_trace, 'Color', [0, 0.5, 0.5], 'LineWidth', 1);  % Using a teal color for the traces
        end
    end
    
    % Calculate the average of all traces, cropping to the shortest length
    % Find the minimum length of traces in the cell array, considering each cell's content
    min_length = min(cellfun(@length, all_traces));  % Get the min length of traces
    
    % Preallocate the array for the average trace
    average_trace = zeros(1, min_length);  % Preallocate with the minimum length
    
    % Initialize a variable to count the number of valid traces
    num_valid_traces = 0;
    
    % Sum the valid traces for averaging
    for i = 1:length(all_traces)
        if ~isempty(all_traces{i}) && length(all_traces{i}) >= min_length
            average_trace = average_trace + all_traces{i}(1:min_length);  % Sum only the first min_length elements
            num_valid_traces = num_valid_traces + 1;  % Count valid traces
        end
    end
    
    % Calculate the average while avoiding division by zero
    if num_valid_traces > 0
        average_trace = average_trace / num_valid_traces;  % Calculate average across valid traces
    else
        error('No valid traces available for averaging.');
    end
    
    % Create a corresponding time vector for the average trace
    average_time = linspace(-pre_saccade_time, post_saccade_time, min_length);  % Create a time vector for the average plot
    
    % Plot the average line
    plot(average_time, average_trace, 'LineWidth', 2);  % Plot average as a dashed line
    
    % Final plot adjustments to show the full range around saccades
    xlabel('Time relative to Saccade Start (s)');
    ylabel('dF/F');
    title('Aligned dF/F Traces Around Saccades');
    grid on;
    xlim([-pre_saccade_time, 3]);  % x-limits to cover the pre and post saccade times
    hold off;

end