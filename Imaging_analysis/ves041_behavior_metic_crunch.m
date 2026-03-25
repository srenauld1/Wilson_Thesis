%%% FOR TESTING DIFFERENT WAYS OF CHARACTERIZING VES041 BEHAVIOR AND
%%% FUNCTION
%%% SHOULD HAVE LOADED 2 FLIES, ONE MENOTAXIS AND ONE WIGGLY

% For this calculation, I will be doing short time intervals of:
    % [0.25, 0.5, 1, 2, 4, 8]
% And long (xcorr) intervals of:
    % [15, 30, 60, 120]

%% CALCULATE FWD VELOCITY X DFF XCORR
% meno
[results_meno_fwd, menotaxis_ves041fly] = moving_window_zero_lag_corr(menotaxis_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'fwd', 'dff');
% wiggle
[results_wiggle_fwd, wiggly_ves041fly] = moving_window_zero_lag_corr(wiggly_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'fwd', 'dff');

%% CALCULATE RHO X DFF XCORR
% meno
[results_meno_rho, menotaxis_ves041fly] = moving_window_zero_lag_corr(menotaxis_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'rho', 'dff');
% wiggle
[results_wiggle_rho, wiggly_ves041fly] = moving_window_zero_lag_corr(wiggly_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'rho', 'dff');

%% CALCULATE SINUOSIY X DFF XCORR
% meno
[results_meno_sinuosity, menotaxis_ves041fly] = moving_window_zero_lag_corr(menotaxis_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'sinuosity', 'dff');
% wiggle
[results_wiggle_sinuosity, wiggly_ves041fly] = moving_window_zero_lag_corr(wiggly_ves041fly, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'sinuosity', 'dff');

%% PLOT THE XCORR
% % fwd
% plot_xcorr_summary(menotaxis_ves041fly.xcorr_fwd, 'fwd')
% plot_xcorr_summary(wiggly_ves041fly.xcorr_fwd, 'fwd')
% % rho
% plot_xcorr_summary(menotaxis_ves041fly.xcorr_rho, 'rho')
% plot_xcorr_summary(wiggly_ves041fly.xcorr_rho, 'rho')
% %sinuosity
% plot_xcorr_summary(menotaxis_ves041fly.xcorr_sinuosity, 'sinuosity')
% plot_xcorr_summary(wiggly_ves041fly.xcorr_sinuosity, 'sinuosity')

plot_xcorr_summary_compare(menotaxis_ves041fly.xcorr_fwd, wiggly_ves041fly.xcorr_fwd, 'fwd', 'menotaxis', 'wiggly')

plot_xcorr_summary_compare(menotaxis_ves041fly.xcorr_rho, wiggly_ves041fly.xcorr_rho, 'rho', 'menotaxis', 'wiggly')

plot_xcorr_summary_compare(menotaxis_ves041fly.xcorr_sinuosity, wiggly_ves041fly.xcorr_sinuosity, 'sinuosity', 'menotaxis', 'wiggly')

%% best options
% xcorr_fwd, win1 = 0.5, win2 = any, we'll choose 15
% xcorr_sinuosity, win1 = 1s, win2 = any but best is 15s

%% plot
savepath_meno = '/Users/sophiarenauld/stacks/ves041_behavior_quantification/menotaxis';
% FWD
menotaxis_ves041fly = plot_flat_path_colored(menotaxis_ves041fly, 'crosscorr', 'whole_range', savepath_meno);
% SINUOSITY
menotaxis_ves041fly = plot_flat_path_colored(menotaxis_ves041fly, 'crosscorr', 'whole_range', savepath_meno);
savepath_wiggle = '/Users/sophiarenauld/stacks/ves041_behavior_quantification/wiggly';
% FWD
wiggly_ves041fly = plot_flat_path_colored(wiggly_ves041fly, 'crosscorr', 'whole_range', savepath_wiggle);
% SINUOSITY
wiggly_ves041fly = plot_flat_path_colored(wiggly_ves041fly, 'crosscorr', 'whole_range', savepath_wiggle);
% DFF
menotaxis_ves041fly = plot_flat_path_colored(menotaxis_ves041fly, 'dff', 'data_only', savepath_meno);
wiggly_ves041fly = plot_flat_path_colored(wiggly_ves041fly, 'dff', 'data_only', savepath_wiggle);



%% FUNCTIONS
function plot_xcorr_summary(xcorr_struct, meas_name)
% xcorr_struct: e.g. menotaxis_ves041fly.xcorr_fwd
% meas_name    : string for plot title/legend, e.g. 'fwd', 'rho', 'sinuosity'

    fns = fieldnames(xcorr_struct);
    fns = fns(~strcmp(fns,'t'));        % remove time field

    nF = numel(fns);
    win1 = nan(nF,1);
    win2 = nan(nF,1);
    mval = nan(nF,1);

    % Extract win1, win2, and mean corr for each field
    for i = 1:nF
        fn = fns{i};  % e.g. 'corr_win1_0_25s_win2_15s'

        % parse using regexp: capture the two numbers after win1_ and win2_
        tokens = regexp(fn, 'corr_win1_([\d_]+)s_win2_([\d_]+)s', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end

        % Convert '0_25' -> 0.25, '15' -> 15, etc.
        w1_str = strrep(tokens{1}, '_', '.');
        w2_str = strrep(tokens{2}, '_', '.');
        win1(i) = str2double(w1_str);
        win2(i) = str2double(w2_str);

        vec = xcorr_struct.(fn);
        mval(i) = mean(vec, 'omitnan');
    end

    % Keep only valid rows
    good = ~isnan(win1) & ~isnan(win2) & ~isnan(mval);
    win1 = win1(good);
    win2 = win2(good);
    mval = mval(good);

    % Unique win2 values (these will define line colors)
    u_w2 = unique(win2);
    nW2 = numel(u_w2);
    colors = lines(nW2);  % one color per win2

    figure; hold on
    for j = 1:nW2
        w2_val = u_w2(j);
        idx = (win2 == w2_val);
        [w1_sorted, order] = sort(win1(idx));
        mv_sorted = mval(idx);
        mv_sorted = mv_sorted(order);

        plot(w1_sorted, mv_sorted, '-o', ...
             'Color', colors(j,:), ...
             'MarkerFaceColor', colors(j,:), ...
             'DisplayName', sprintf('win2 = %.0f s', w2_val));
    end

    ylim([-1 1]);
    xlabel('Window 1 (s)');
    ylabel('Mean lag-0 correlation');
    title(sprintf('Mean xcorr vs win1 (%s)', meas_name));
    legend('Location','best');
    grid on;
end

function plot_xcorr_summary_compare(xcorr_struct1, xcorr_struct2, meas_name, label1, label2)
% Compare mean xcorr vs win1 for two datasets on same plot.
%
% xcorr_struct1  : e.g. fly1.xcorr_fwd
% xcorr_struct2  : e.g. fly2.xcorr_fwd
% meas_name      : e.g. 'fwd', 'rho', 'sinuosity'
% label1, label2 : legend labels for the two datasets

    % --- helper to extract win1, win2, mean for one struct ---
    function [win1, win2, mval] = extract_means(S)
        fns = fieldnames(S);
        fns = fns(~strcmp(fns,'t'));  % remove 't'
        nF = numel(fns);
        win1 = nan(nF,1);
        win2 = nan(nF,1);
        mval = nan(nF,1);
        for i = 1:nF
            fn = fns{i};
            tokens = regexp(fn,'corr_win1_([\d_]+)s_win2_([\d_]+)s','tokens','once');
            if isempty(tokens), continue; end
            w1_str = strrep(tokens{1},'_','.');
            w2_str = strrep(tokens{2},'_','.');
            win1(i) = str2double(w1_str);
            win2(i) = str2double(w2_str);
            vec = S.(fn);
            mval(i) = mean(vec,'omitnan');
        end
        good = ~isnan(win1) & ~isnan(win2) & ~isnan(mval);
        win1 = win1(good);
        win2 = win2(good);
        mval = mval(good);
    end

    % --- Extract from both structs ---
    [w1_1, w2_1, m1] = extract_means(xcorr_struct1);
    [w1_2, w2_2, m2] = extract_means(xcorr_struct2);

    % For safety, restrict to common (win1, win2) combos
    combos1 = [w1_1, w2_1];
    combos2 = [w1_2, w2_2];
    [common, ia, ib] = intersect(combos1, combos2, 'rows');

    w1_1 = w1_1(ia);  w2_1 = w2_1(ia);  m1 = m1(ia);
    w1_2 = w1_2(ib);  w2_2 = w2_2(ib);  m2 = m2(ib);

    u_w2 = unique(common(:,2));
    colors = lines(numel(u_w2));

    figure; hold on
    for j = 1:numel(u_w2)
        this_w2 = u_w2(j);

        % dataset 1
        idx1 = (w2_1 == this_w2);
        [w1s_1, ord1] = sort(w1_1(idx1));
        mvals_1 = m1(idx1); mvals_1 = mvals_1(ord1);

        % dataset 2
        idx2 = (w2_2 == this_w2);
        [w1s_2, ord2] = sort(w1_2(idx2));
        mvals_2 = m2(idx2); mvals_2 = mvals_2(ord2);

        % Plot: same color per win2, but different style/marker per dataset
        plot(w1s_1, mvals_1, '-o', 'Color', colors(j,:), ...
             'MarkerFaceColor', colors(j,:), ...
             'DisplayName', sprintf('%s, win2 = %.0fs', label1, this_w2));
        plot(w1s_2, mvals_2, '--s', 'Color', colors(j,:), ...
             'MarkerFaceColor', 'w', ...
             'DisplayName', sprintf('%s, win2 = %.0fs', label2, this_w2));
    end

    ylim([-1 1]);
    xlabel('Window 1 (s)');
    ylabel('Mean lag-0 correlation');
    title(sprintf('Mean xcorr vs win1 (%s): %s vs %s', meas_name, label1, label2));
    legend('Location','bestoutside');
    grid on;
end

