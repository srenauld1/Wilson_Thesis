function a2p_data = detect_local_peaks_2p(a2p_data, varargin)
    % Parse inputs
    p = inputParser;
    addRequired(p, 'ts');
    addParameter(p, 'InitTime', 0);
    addParameter(p, 'Prominence', 0.1);
    addParameter(p, 'MinTimeGap', 60);
    parse(p, a2p_data, varargin{:});
    
    init_t = p.Results.InitTime;
    prominence = p.Results.Prominence;
    min_time_gap = p.Results.MinTimeGap;
    
    % Filter data after init_t
    time_values = a2p_data.dq(1).t;
    valid_times = time_values > init_t;
    filtered_time = time_values(valid_times);
    filtered_signal = a2p_data.dq(1).absolute_circular_diff(valid_times);
    
    % Compute time-based distance in number of samples
    avg_sampling_interval = mean(diff(filtered_time));
    min_samples_gap = round(min_time_gap / avg_sampling_interval);
    
    % Find peaks
    [~, peak_locs] = findpeaks(filtered_signal, ...
                              'MinPeakProminence', prominence, ...
                              'MinPeakDistance', min_samples_gap);
    
    % Convert peak locations back to original time indices
    orig_indices = find(valid_times);
    peak_indices = orig_indices(peak_locs);
    
    % Create jump_detected field
    a2p_data.jump_detected = zeros(size(time_values));
    a2p_data.jump_detected(peak_indices) = 1;
end