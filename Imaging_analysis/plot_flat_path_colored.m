function a2p_data = plot_flat_path_colored(a2p_data, color_option, scale_option, savepath, odor)
% Flat path colored by chosen variable
% scale_option - 'whole_range' or 'data_only'

% ---------------- Validate scale_option ----------------
if isstring(scale_option)
    scale_option = char(scale_option);
end
valid_scales = {'whole_range','data_only'};
if ~ismember(scale_option, valid_scales)
    error('scale_option must be either ''whole_range'' or ''data_only''.');
end

% ---------------- Colormap: blue->red->yellow ----------------
n = 100;
half = round(n/2);
cmap1 = [linspace(0,1,half)' zeros(half,1) linspace(1,0,half)'];
cmap2 = [ones(n-half,1) linspace(0,1,n-half)' zeros(n-half,1)];
cmap = [cmap1; cmap2];

% ---------------- Select variable to color by ----------------
if strcmp(color_option, 'dff')
    color_scale = a2p_data.roi.ts(:);
    if odor
        ball_px = a2p_data.dq(1).px(:);
        ball_py = a2p_data.dq(1).py(:);
    else
        ball_px = a2p_data.dq(1).pxb(:);
        ball_py = a2p_data.dq(1).pyb(:);
    end

elseif strcmp(color_option, 'rho')
    color_scale = a2p_data.menotaxis.rho(:);
    if odor
        ball_px = a2p_data.dq(2).px(:);
        ball_py = a2p_data.dq(2).py(:);
    else
        ball_px = a2p_data.dq(2).pxb(:);
        ball_py = a2p_data.dq(2).pyb(:);
    end

elseif strcmp(color_option, 'sinuosity')
    color_scale = a2p_data.sinuosity(:);
    if odor
        ball_px = a2p_data.dq(2).px(:);
        ball_py = a2p_data.dq(2).py(:);
    else
        ball_px = a2p_data.dq(2).pxb(:);
        ball_py = a2p_data.dq(2).pyb(:);
    end

elseif strcmp(color_option, 'fwd')
    color_scale = a2p_data.dq(2).bvf(:);
    if odor
        ball_px = a2p_data.dq(2).px(:);
        ball_py = a2p_data.dq(2).py(:);
    else
        ball_px = a2p_data.dq(2).pxb(:);
        ball_py = a2p_data.dq(2).pyb(:);
    end

elseif strcmp(color_option, 'crosscorr')
    % ---------- 1. Choose which xcorr set ----------
    xcorr_choices = {'fwd', 'sinuosity', 'rho'};
    [indx, tf] = listdlg('PromptString', 'Select xcorr type to plot:', ...
                         'SelectionMode', 'single', ...
                         'ListString', xcorr_choices);
    if ~tf
        error('No cross-correlation measurement selected.');
    end
    chosen_type = xcorr_choices{indx};

    fieldname = ['xcorr_' chosen_type]; % e.g. 'xcorr_fwd'
    if ~isfield(a2p_data, fieldname)
        error('No results found for chosen type: %s', chosen_type);
    end
    results_struct = a2p_data.(fieldname);

    % ---------- 2. Parse available win1/win2 combos ----------
    fnames = fieldnames(results_struct);
    fnames = fnames(~strcmp(fnames,'t'));   % remove time field, if present

    win1_all = [];
    win2_all = [];
    for i = 1:numel(fnames)
        tok = regexp(fnames{i}, 'corr_win1_([\d_]+)s_win2_([\d_]+)s', 'tokens', 'once');
        if isempty(tok), continue; end
        w1 = str2double(strrep(tok{1}, '_', '.'));
        w2 = str2double(strrep(tok{2}, '_', '.'));
        win1_all(end+1,1) = w1; %#ok<AGROW>
        win2_all(end+1,1) = w2; %#ok<AGROW>
    end
    if isempty(win1_all)
        error('No corr_win1_* fields found for %s.', fieldname);
    end

    u_win1 = unique(win1_all);
    u_win2 = unique(win2_all);

    % ---------- 3. Ask user for window1 & window2 ----------
    win1_str = arrayfun(@(x) sprintf('%.3g s', x), u_win1, 'UniformOutput', false);
    [i1, tf1] = listdlg('PromptString', 'Select Window 1 (sec):', ...
                        'SelectionMode','single', ...
                        'ListString', win1_str);
    if ~tf1, error('No Window1 selected.'); end
    chosen_win1 = u_win1(i1);

    win2_str = arrayfun(@(x) sprintf('%.0f s', x), u_win2, 'UniformOutput', false);
    [i2, tf2] = listdlg('PromptString', 'Select Window 2 (sec):', ...
                        'SelectionMode','single', ...
                        'ListString', win2_str);
    if ~tf2, error('No Window2 selected.'); end
    chosen_win2 = u_win2(i2);

    % ---------- 4. Find the matching field for chosen win1/win2 ----------
    target_pattern = sprintf('corr_win1_%gs_win2_%gs', chosen_win1, chosen_win2);
    target_pattern = strrep(target_pattern, '.', '_'); % field name uses '_' for decimals
    match_idx = contains(fnames, target_pattern, 'IgnoreCase', true);
    if ~any(match_idx)
        error('No field found for win1=%.3g s, win2=%.3g s', chosen_win1, chosen_win2);
    end
    data_fn = fnames{find(match_idx,1,'first')};

    % ---------- 5. Extract the selected xcorr trace ----------
    color_scale = results_struct.(data_fn)(:);
    a2p_data.xcorr_answer = color_scale;  % store for later use

    if odor
        ball_px = a2p_data.dq(1).px(:);
        ball_py = a2p_data.dq(1).py(:);
    else
        ball_px = a2p_data.dq(1).pxb(:);
        ball_py = a2p_data.dq(1).pyb(:);
    end


else
    error('Unknown color_option: %s', color_option);
end




% ---------------- Clip extremes at N·sigma (keep all points) ----------------
mu = mean(color_scale, 'omitnan');
sigma = std(color_scale, 'omitnan');
Nsig = 2;  % cutoff in sigmas

low_thr  = mu - Nsig*sigma;
high_thr = mu + Nsig*sigma;

color_scale_clipped = color_scale;
color_scale_clipped(color_scale_clipped < low_thr)  = low_thr;
color_scale_clipped(color_scale_clipped > high_thr) = high_thr;

color_scale_f = color_scale_clipped;   % used for mapping
pxb_f = ball_px;
pyb_f = ball_py;

% ---------------- Compute mapping to colormap ----------------
use_whole_range = strcmp(scale_option, 'whole_range') && strcmp(color_option, 'crosscorr');

if use_whole_range
    % For crosscorr “whole_range”: use raw values as CData, fix [-1 1]
    cdata = color_scale_f;
else
    % Linear scaling of (possibly clipped) data -> [0,1] -> colormap indices
    min_f = min(color_scale_f,[],'omitnan');
    max_f = max(color_scale_f,[],'omitnan');
    if max_f == min_f
        linear_scaled = zeros(size(color_scale_f));
    else
        linear_scaled = (color_scale_f - min_f) / (max_f - min_f);
    end
    color_inds = round(linear_scaled * (size(cmap,1)-1)) + 1;
    color_inds = max(min(color_inds, size(cmap,1)), 1);
end

% ---------------- Plot ----------------
figure;
hold on;

if use_whole_range
    scatter(pxb_f, pyb_f, 10, cdata, 'filled', 'MarkerEdgeColor', 'none');
    colormap(cmap);
    c = colorbar;
    caxis([-0.3 0.3]);
    ylabel(c, [color_option ' (r)']);
else
    scatter(pxb_f, pyb_f, 10, cmap(color_inds,:), 'filled', 'MarkerEdgeColor', 'none');
    colormap(cmap);
    c = colorbar;
    caxis([1 n]);  % index space

    yticks_vals = get(c,'YTick');
    actual_scaled = (yticks_vals - 1) / (n - 1);  % 0..1
    min_f = min(color_scale_f,[],'omitnan');
    max_f = max(color_scale_f,[],'omitnan');
    actual_values = actual_scaled * (max_f - min_f) + min_f;

    ylabel(c, color_option);
    set(c,'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false));
end

min_value = min([ball_px(:); ball_py(:)]);
max_value = max([ball_px(:); ball_py(:)]);
xlim([min_value max_value]);
ylim([min_value max_value]);

% Start point
plot(ball_px(1), ball_py(1), 'r.', 'MarkerSize', 20);

% ---------------- Title and filenames ----------------
if isstring(color_option)
    color_str = char(color_option);
else
    color_str = color_option;
end

if strcmp(color_option, 'crosscorr')
    ttl = sprintf('Flat path, colored by %s %s (win1=%.3gs, win2=%.0fs)', ...
                  color_option, chosen_type, chosen_win1, chosen_win2);
    title(ttl);

    w1_str_save = strrep(sprintf('%.3g', chosen_win1), '.', 'p');
    w2_str_save = strrep(sprintf('%.0f', chosen_win2), '.', 'p');

    xfilename = sprintf('x_crosscorr_%s_w1_%ss_w2_%ss', ...
                        chosen_type, w1_str_save, w2_str_save);
    yfilename = sprintf('y_crosscorr_%s', chosen_type);
else
    title("Flat path, colored by " + color_option + " (linear scale)");
    xfilename = ['x_color_linear_' color_str];
    yfilename = ['y_color_linear_' color_str];
end

save_plot_with_title_as_filename(xfilename, yfilename, savepath);
end