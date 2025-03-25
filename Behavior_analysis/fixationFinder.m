% FIXATIONFINDER - This function identifies periods of fixation and running behavior in flies
% based on panel position (yaw) and forward velocity data. Fixation is determined by the moving
% mean of unwrapped panel position, with thresholds applied to ensure that fixations meet duration
% and gap criteria. The function also identifies periods where the fly is both fixating and running.

% INPUTS:
%   panelps    - 3D array of panel position data, where each slice represents a condition, 
%                and each column represents a trial.
%   forward    - 3D array of forward velocity data corresponding to each trial.
%   expttime   - Time vector corresponding to the temporal resolution of the data.
%   optPlot    - Optional flag (1/0) to plot the results for a sanity check (1 for yes, 0 for no).

% OUTPUT:
%   fixation   - Structure containing fixation indices and filtered panel/forward data:
%                - fixation.panelps_all : Panel position data where the fly was fixating.
%                - fixation.forward_all : Forward velocity data where the fly was fixating.
%                - fixation.idx_all     : Logical indices indicating fixation periods.
%                - fixation.panelps_run : Panel position data where the fly was both fixating and running.
%                - fixation.forward_run : Forward velocity data where the fly was both fixating and running.
%                - fixation.idx_run     : Logical indices indicating fixation and running periods.

% Created: N/A by MC
% Updated: N/A

% The function:
% - Unwraps the panel position to avoid discontinuities and computes a moving mean.
% - Applies thresholds for fixation based on the range of the moving mean.
% - Identifies fixation periods that meet minimum duration criteria and merges adjacent
%   fixations if the gap between them is below the maximum allowed gap.
% - Detects running behavior based on forward velocity and combines it with fixation indices
%   to find periods where the fly is both fixating and running.
% - Optionally plots the panel position and fixation periods for visualization.
%
function fixation = fixationFinder(panelps, forward, expttime, optPlot)
    %% initialize
    durTrial = size(panelps, 1);
    nTrial = size(panelps, 2);
    nCond = size(panelps, 3);

    % fetch settings
    settings = processSettings();

    % set boxcar settings
    boxcar_window = 2; %s, size of window

    % set fixation requirements
    mean_max = 35; % max range of boxcar means to be considered fixation
    min_fix_dur = 5; % Minimum fixation duration in seconds
    max_gap_dur = 2; % maximum non-fixation gap duration to merge adjacent fixations
    max_run_gap_dur = 2; % maximum non-run gap duration to merge adjacent run bouts
    %mrl_threshold = 0.5; % Minimum MRL for valid fixation
        
    %% threshold for when the fly was running
    runIdx = zeros(durTrial, nTrial, nCond);
    for c = 1:nCond
        runIdx(:, :, c) = schmittTrigger(forward(:, :, c), settings.runThreshB, 0.1);
    end

    %% threshold for when the fly was fixating
    % Unwrap panelps to prevent jumps between -180 and 180 degrees
    panelps_unwrapped = unwrap(deg2rad(panelps)); % Convert to radians and unwrap

    % Calculate the moving mean on unwrapped data
    bcIdx_win = fetchTimeIdx(expttime, boxcar_window);
    meanpanelps_unwrapped = movmean(panelps_unwrapped, bcIdx_win, 'Endpoints', 'fill');

    % Rewrap the result back into [-pi, pi] range and then convert to degrees
    meanpanelps = rad2deg(wrapToPi(meanpanelps_unwrapped)); % Rewrap to [-180, 180]

    % Calculate MRL over the same moving window
    % mean_cos = movmean(cos(panelps_unwrapped), bcIdx_win, 'Endpoints', 'fill');
    % mean_sin = movmean(sin(panelps_unwrapped), bcIdx_win, 'Endpoints', 'fill');
    % mrl = sqrt(mean_cos.^2 + mean_sin.^2); % MRL calculation

    % Set fixation indices based on the circular moving mean
    fixIdx = (abs(meanpanelps) < mean_max);
    
    % Process fixation indices to apply duration criteria
    fixIdx = processFixationPeriods(fixIdx, expttime, max_gap_dur, min_fix_dur);

    %% store for output
    % store ALL data where fly was fixating
    panelps_fix = panelps;
    panelps_fix(~fixIdx) = nan;
    forward_fix = forward;
    forward_fix(~fixIdx) = nan;

    fixation.panelps_all = panelps_fix;
    fixation.forward_all = forward_fix;
    fixation.idx_all = logical(fixIdx);

    % store data where fly was fixating AND running
    fixIdxR = fixIdx & runIdx;
    fixIdxR = processFixationPeriods(fixIdxR, expttime, max_run_gap_dur, min_fix_dur);

    panelps_fix_r = panelps;
    panelps_fix_r(~fixIdxR) = nan;
    forward_fix_r = forward;
    forward_fix_r(~fixIdxR) = nan;

    fixation.panelps_run = panelps_fix_r;
    fixation.forward_run = forward_fix_r;
    fixation.idx_run = logical(fixIdxR);

    % estimate time spent fixating
    timeFixating = reshape(sum(expttime(sum(~isnan(panelps_fix_r)) + 1)),1,nCond);

    %% Plot results for sanity check
    if optPlot
        % initialize
        figure; set(gcf, 'Position', [50 50 4000 900])
        tiledlayout(5, 5, 'TileSpacing', 'tight')
        x = 1;
        if timeFixating > settings.minFixationTime
            fixcolor = "#77AC30";
        else
            fixcolor = "#A2142F";
        end
        % Fixation Indices Plot
        for c = 1:5
            % estimate HD bias
            biasHD = mean(panelps_fix_r(:,:,c), 'all', 'omitnan');

            for t = 1:5
                nexttile; hold on

                % Subtract biasHD and wrap to -180 to 180
                wrapped_panelps = mod(panelps(:, t, c) - biasHD + 180, 360) - 180;
                wrapped_panelps_fix = mod(panelps_fix(:, t, c) - biasHD + 180, 360) - 180;
                wrapped_panelps_fix_r = mod(panelps_fix_r(:, t, c) - biasHD + 180, 360) - 180;

                % Plot each adjusted trace
                plot(expttime, wrapped_panelps, 'Color', [0.5 0.5 0.5]);
                plot(expttime, wrapped_panelps_fix, 'k', 'LineWidth', 1);
                plot(expttime, wrapped_panelps_fix_r, 'Color', fixcolor, 'LineWidth', 1);

                axis tight;
                ylim([-180 180]);
                yline(0);
                set(gca, 'XTick', []);
                set(gca, 'YTick', []);
                x = x + 1;
                if t == 1
                    ylabel(num2str(settings.pursuitGain(c)))
                end
            end
            legend(['Fixation Time: ' num2str(round(timeFixating(c)))])
        end
    end

end
