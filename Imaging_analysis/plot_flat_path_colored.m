function a2p_data = plot_flat_path_colored(a2p_data, color_option, savepath)
% Flat path colored by chosen variable (linear scaling of color variable)
n = 100;
half = round(n/2);

% Blue to red (first half)
cmap1 = [linspace(0,1,half)' zeros(half,1) linspace(1,0,half)'];

% Red to yellow (second half)
cmap2 = [ones(n-half,1) linspace(0,1,n-half)' zeros(n-half,1)];

% Combine colormaps
cmap = [cmap1; cmap2];

if strcmp(color_option, 'dff')
    color_scale = a2p_data.roi.ts;
    ball_px = a2p_data.dq(1).pxb;
    ball_py = a2p_data.dq(1).pyb;
elseif strcmp(color_option, 'rho')
    color_scale = a2p_data.menotaxis.rho;
    ball_px = a2p_data.dq(2).pxb;
    ball_py = a2p_data.dq(2).pyb;
elseif strcmp(color_option, 'sinuosity')
    color_scale = a2p_data.sinuosity;
    ball_px = a2p_data.dq(2).pxb;
    ball_py = a2p_data.dq(2).pyb;
end

% == Linear scaling (no sqrt, no log) ==
linear_scaled = (color_scale - min(color_scale)) / (max(color_scale) - min(color_scale));

% Convert scaled values to colormap indices
color_inds = round(linear_scaled * (n - 1)) + 1;

figure
hold on
scatter(ball_px, ball_py, 10, cmap(color_inds, :), 'filled', 'MarkerEdgeColor', 'none');
min_value = min([ball_px(:); ball_py(:)]);
max_value = max([ball_px(:); ball_py(:)]);
xlim([min_value max_value]);
ylim([min_value max_value]);
colormap(cmap); % <-- THIS SETS THE FIGURE COLORBAR TO MATCH YOUR DATA

% Create colorbar with actual input variable values
c = colorbar;
caxis([1 n])  % Set colorbar limits to match normalized range
yticks_vals = get(c, 'YTick');  % Get current tick positions

% Convert colorbar ticks back to actual variable values (linear)
actual_scaled = (yticks_vals - 1) / (n - 1);  % Map ticks back to [0,1]
actual_values = actual_scaled * (max(color_scale) - min(color_scale)) + min(color_scale);

ylabel(c, color_option)
set(c, 'YTickLabel', arrayfun(@(x) sprintf('%.3f', x), actual_values, 'UniformOutput', false))

% Start point
plot(ball_px(1), ball_py(1), 'r.', 'MarkerSize', 20)

title("Flat path, colored by " + color_option + " (linear scale)")
% Make sure color_option is a character vector or string for concatenation
if isstring(color_option)
    color_str = char(color_option);
else
    color_str = color_option;
end

% Use it in file name
save_plot_with_title_as_filename(['x_color_linear_' color_str], ...
                                 ['y_color_linear_' color_str], ...
                                 savepath);
end