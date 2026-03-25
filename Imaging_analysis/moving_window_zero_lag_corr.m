function [results, a2p_data] = moving_window_zero_lag_corr(a2p_data, win1_secs, win2_secs, measurement, activity_field)
% For each point, computes windowed Pearson correlation (lag 0) between measure & activity
% Inputs:
%   - a2p_data: data struct
%   - win_sec1/ 2: window size in seconds (scalar) for behavior (1) &
%   xcorr(2)
%   - measurement: 'fwd', 'rho', or 'sinuosity'
%   - activity_field: e.g. 'dff'
% Output:
%   - xcorr0_trace: one value per timepoint (sliding window lag 0 correlation)
    % Get activity and timebases

    % Activity and timebases
    activity_t = a2p_data.dq(1).t(:);           % neural/roi timebase
    if strcmp(activity_field, 'dff')
        activity = a2p_data.roi.ts;
    end


    behavior_t = a2p_data.dq(2).t(:);

    dt_behav = median(diff(behavior_t));
    dt_activity = median(diff(activity_t));
    n = numel(activity_t);

   results = struct;

    for w1 = 1:length(win1_secs)
        win1_sec = win1_secs(w1);

        % -- Measurement computation with win1 smoothing --
        switch lower(measurement)
            case {'forward_velocity','fwd','fwd_vel'}
                raw_signal = a2p_data.dq(2).bvf(:);
                win1_pts = max(3, round(win1_sec / dt_behav));
                smoothed_signal = movmean(raw_signal, win1_pts, 'omitnan');
            case 'rho'
                % Compute with chosen window
                [~, ~, raw_signal, ~, ~, ~] = SegmentMenovsNotMeno_2p(a2p_data, '', win1_sec, 1, 0.88, 0.88, 0, 0);
                smoothed_signal = raw_signal(:);
            case 'sinuosity'
                win1_pts = max(3, round(win1_sec / dt_behav));
                a2p_data = add_sinuosity_sliding(a2p_data, win1_pts);
                smoothed_signal = a2p_data.sinuosity(:);
        end

        smoothed_signal_interp = interp1(behavior_t, smoothed_signal, activity_t, 'linear', 'extrap');
        
        for w2 = 1:length(win2_secs)
            win2_sec = win2_secs(w2);
            halfwin2_pts = round(win2_sec/(2*dt_activity));
            xcorr0_trace = nan(n,1);

            for i = 1:n
                w_start = max(1, i-halfwin2_pts);
                w_end   = min(n, i+halfwin2_pts);
                x  = smoothed_signal_interp(w_start:w_end);
                y  = activity(w_start:w_end)';
                if nnz(x)>3
                    % Always supply column vectors to corr!
                    xcorr0_trace(i) = corr(x, y);
                end
            end

            % Store result for this window combination
            fieldname = strrep(sprintf('corr_win1_%gs_win2_%gs', win1_sec, win2_sec), '.', '_');
            results.(fieldname) = xcorr0_trace;
            
            % Print mean value for this combo:
            fprintf('win1 = %.2fs, win2 = %.2fs: mean corr = %.4f\n', win1_sec, win2_sec, nanmean(xcorr0_trace));
        end
    end
    % Also store timebase
    results.t = activity_t;
    a2p_data.(['xcorr_' measurement]) = results;
end