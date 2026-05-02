function plot_flatpath_opto(exptData, savepath)
% Plot flat path: opto ON = RED, opto OFF = BLACK, start = GREEN
% Expects exptData as output from process_fictrac_panels:
%   x, y, optoStim are already aligned at the same length (800K)

% Validate required fields
required_fields = {'x', 'y'};
for f = required_fields
    if ~isfield(exptData, f{1})
        error('plot_flatpath_opto: exptData missing required field: %s', f{1});
    end
end

x = exptData.x;
y = exptData.y;

if isfield(exptData, 'optoStim_ds')
    opto_timing = exptData.optoStim_ds;
else
    warning('plot_flatpath_opto: No optoStim field found. Plotting all as opto OFF.');
    opto_timing = zeros(size(x));
end

% Downsample just x, y, opto_timing for plotting performance (800K → 80K)
ds_factor   = 10;
x           = x(1:ds_factor:end);
y           = y(1:ds_factor:end);
opto_timing = opto_timing(1:ds_factor:end);

% NaN-mask approach: much faster than per-segment loop
on_mask  = opto_timing ~= 0;
off_mask = ~on_mask;

x_on  = x;  x_on(off_mask) = NaN;
y_on  = y;  y_on(off_mask) = NaN;
x_off = x;  x_off(on_mask) = NaN;
y_off = y;  y_off(on_mask) = NaN;

figure; hold on;

h_noopto = plot(x_off, y_off, 'k-', 'LineWidth', 1.3);
h_opto   = plot(x_on,  y_on,  '-',  'Color', [1 0 0], 'LineWidth', 1.3);
h_start  = plot(x(1),  y(1),  'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');

xlabel('X Position');
ylabel('Y Position');
title('Flat Path: Opto ON (Red), Opto OFF (Black), Start (Green)');

% Build legend only from handles that actually have data
legend_handles = [];
legend_labels  = {};
if any(~isnan(x_on))
    legend_handles(end+1) = h_opto;
    legend_labels{end+1}  = 'Opto ON';
end
if any(~isnan(x_off))
    legend_handles(end+1) = h_noopto;
    legend_labels{end+1}  = 'Opto OFF';
end
legend_handles(end+1) = h_start;
legend_labels{end+1}  = 'Start';

legend(legend_handles, legend_labels, 'Location', 'best');

axis equal; axis tight; grid on; hold off;

save_plot_with_title_as_filename('flatpath', 'opto', savepath)
end