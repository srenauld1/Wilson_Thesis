function save_plot_with_title_as_filename(ylabel_left, ylabel_right, savepath)
    filename = sprintf('%s_vs_%s.fig', strrep(ylabel_left, ' ', '_'), strrep(ylabel_right, ' ', '_'));
    filename_png = sprintf('%s_vs_%s.png', strrep(ylabel_left, ' ', '_'), strrep(ylabel_right, ' ', '_'));
    fullfilename = fullfile(savepath, filename);
    fullfile_png = fullfile(savepath, filename_png);

    % --- Extract last non-leaf directory from savepath ---
    parentdir = fileparts(savepath); % removes 'smr_analysis'
    [~, id_label] = fileparts(parentdir); % gets '20251001-21_ves041_fly3_clbbar'

    % -- Place id_label as annotation in lower-right of figure window --
    ax = gca;
    % Normalized position: right edge, below axes (x=0.98, y=0.02: tweak as needed)
    annotation(gcf, 'textbox', [0.80 0.01 0.18 0.04], ...
        'String', id_label, ...
        'FontSize', 11, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', ...           % No border
        'BackgroundColor', 'none', ...     % No background fill
        'Color', 'k', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none', ...
        'FitBoxToText', 'on');             % Resize box to fit text

    savefig(gcf, fullfilename);
    saveas(gcf, fullfile_png);
    fprintf('Plot saved as: %s\n', fullfilename);
end