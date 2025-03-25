% Function to save the plot with y-axis labels as title
function save_plot_with_title_as_filename(ylabel_left, ylabel_right, savepath)
    % Get current figure title
    current_title = get(get(gca, 'Title'), 'String');
    
    % Generate a valid filename using y-axis labels
    filename = sprintf('%s_vs_%s.fig', strrep(ylabel_left, ' ', '_'), strrep(ylabel_right, ' ', '_'));
    filename_png = sprintf('%s_vs_%s.png', strrep(ylabel_left, ' ', '_'), strrep(ylabel_right, ' ', '_'));
    fullfilename = fullfile(savepath, filename);
    % Save the figure
    savefig(gcf, fullfilename);
    %also save as PNG
    fullfile_png = fullfile(savepath, filename_png);
    saveas(gcf, fullfile_png);
    
    % Print a message indicating the file was saved
    fprintf('Plot saved as: %s\n', fullfilename);
end