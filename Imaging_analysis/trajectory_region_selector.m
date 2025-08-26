function trajectory_region_selector(daq, dff, jump, varargin)
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'moving_only', false, @islogical);
    addParameter(p, 'title_prefix', 'Trajectory', @ischar);
    addParameter(p, 'save_path', pwd, @ischar);
    parse(p, varargin{:});
    
    % Filter data if needed
    if p.Results.moving_only && isfield(daq, 'motion') && isfield(daq.motion, 'moving_or_not')
        moving_mask = logical(daq.motion.moving_or_not);
        px = daq.px(moving_mask);
        py = daq.py(moving_mask);
        dff_plot = dff(moving_mask);
        time_dff_full = daq.t(moving_mask);
        original_indices = find(moving_mask);
    else
        px = daq.px;
        py = daq.py;
        dff_plot = dff;
        time_dff_full = daq.t;
        original_indices = (1:length(px))';
    end
    
    % Create figure
    fig = figure('Position', [100, 100, 1400, 700], 'Name', 'Trajectory Region Selector');
    
    % Main trajectory plot (left side)
    ax_main = subplot('Position', [0.05, 0.15, 0.4, 0.75]);
    scatter(px, py, 36, 'b', 'filled', 'MarkerEdgeColor', 'none');
    hold on
    if jump
    plot(daq.px(logical(daq.jump_detected)), ...
         daq.py(logical(daq.jump_detected)), ...
         'm.', 'MarkerSize', 15)  % Adjust color and size as needed
    end
    title([p.Results.title_prefix ' - Draw polygon regions'])
    xlabel('X Position')
    ylabel('Y Position')
    axis equal
    grid on
    
    % Dual-axis plot area (right side)
    ax_signals = subplot('Position', [0.55, 0.4, 0.4, 0.5]);
    plot(0, 0); % Placeholder
    title('DFF and BYV\_deg vs Time in Selected Regions')
    xlabel('Time (s)')
    grid on
    
    % Summary text area (bottom right)
    ax_summary = subplot('Position', [0.55, 0.15, 0.4, 0.2]);
    axis off
    
    % Initialize storage
    regions = {};
    region_handles = {};
    selected_points = {};
    region_colors = ['r', 'g', 'm', 'c', 'y', 'k']; % Colors for different regions
    
    % Create UI controls (two rows)
    % First row
    uicontrol('Style', 'pushbutton', 'String', 'Add Polygon Region', ...
        'Position', [20, 60, 150, 40], 'Callback', @add_polygon_region, ...
        'FontSize', 12);
    
    uicontrol('Style', 'pushbutton', 'String', 'Clear All Regions', ...
        'Position', [180, 60, 150, 40], 'Callback', @clear_regions, ...
        'FontSize', 12);
    
    % Second row
    uicontrol('Style', 'pushbutton', 'String', 'Save Figure', ...
        'Position', [20, 20, 150, 40], 'Callback', @save_figure, ...
        'FontSize', 12, 'BackgroundColor', [0.8, 0.9, 0.8]);
    
    uicontrol('Style', 'pushbutton', 'String', 'Export Results', ...
        'Position', [180, 20, 150, 40], 'Callback', @export_results, ...
        'FontSize', 12);
    
    % Text for summary
    summary_text = text(0.05, 0.8, 'Draw polygon regions to see DFF and BYV_deg analysis', ...
        'Parent', ax_summary, 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'FontName', 'FixedWidth');
    
    function add_polygon_region(~, ~)
        axes(ax_main);
        
        % Draw interactive polygon
        region_num = length(regions) + 1;
        color = region_colors(mod(region_num-1, length(region_colors)) + 1);
        
        h = drawpolygon('Color', color, 'LineWidth', 2);
        
        if isvalid(h) && size(h.Position, 1) >= 3
            % Get polygon vertices
            vertices = h.Position;
            
            % Find points inside polygon (spatial selection)
            inside_mask = inpolygon(px, py, vertices(:,1), vertices(:,2));
            
            % Get time range of selected spatial points
            selected_times_dff = time_dff_full(inside_mask);
            time_min = min(selected_times_dff);
            time_max = max(selected_times_dff);
            
            % Find corresponding time windows in both datasets
            dff_time_mask = time_dff_full >= time_min & time_dff_full <= time_max;
            supp_time_mask = daq.t_supp >= time_min & daq.t_supp <= time_max;
            
            % Extract data based on time alignment
            region_dff = dff_plot(dff_time_mask);
            region_time_dff = time_dff_full(dff_time_mask);
            region_time_supp = daq.t_supp(supp_time_mask);
            
            % Handle different field names for yaw velocity
            if isfield(daq, 'byv_deg_supp')
                region_byv_deg = daq.byv_deg_supp(supp_time_mask);
            elseif isfield(daq, 'byv_supp')
                region_byv_deg = daq.byv_supp(supp_time_mask);
            else
                error('Neither daq.byv_deg_supp nor daq.byv_supp found');
            end
            
            % Store region info
            regions{region_num} = struct('type', 'polygon', 'vertices', vertices, 'color', color, ...
                'time_range', [time_min, time_max]);
            region_handles{region_num} = h;
            selected_points{region_num} = struct(...
                'spatial_indices', original_indices(inside_mask), ...
                'spatial_x', px(inside_mask), ...
                'spatial_y', py(inside_mask), ...
                'spatial_count', sum(inside_mask), ...
                'dff', region_dff, ...
                'byv_deg', region_byv_deg, ...
                'time_dff', region_time_dff, ...
                'time_supp', region_time_supp, ...
                'time_range', [time_min, time_max], ...
                'dff_count', length(region_dff), ...
                'supp_count', length(region_byv_deg));
            
            % Plot selected spatial points with same color as region
            plot(px(inside_mask), py(inside_mask), 'o', 'Color', color, ...
                'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor', color);
            
            % Update displays
            update_signals_plot();
            update_summary_display();
            
            % Print info to console
            fprintf('Region %d: Spatial points = %d, DFF points = %d, BYV points = %d\n', ...
                region_num, sum(inside_mask), length(region_dff), length(region_byv_deg));
            fprintf('  Time range: %.2f to %.2f s\n', time_min, time_max);
        end
    end
    
    function clear_regions(~, ~)
        % Clear all region handles
        for i = 1:length(region_handles)
            if isvalid(region_handles{i})
                delete(region_handles{i});
            end
        end
        
        % Clear all plotted points (keep original scatter)
        axes(ax_main);
        children = get(gca, 'Children');
        % Keep only the first child (original scatter plot)
        if length(children) > 1
            delete(children(1:end-1));
        end
        
        % Reset storage
        regions = {};
        region_handles = {};
        selected_points = {};
        
        % Properly clear the dual-axis signals plot
        axes(ax_signals);
        cla reset;  % This clears all axes properties and content
        
        % Reset to single axis and add placeholder
        yyaxis left;
        cla;
        yyaxis right;
        cla;
        
        % Reset to left axis and set up clean state
        yyaxis left;
        plot(0, 0, 'Color', 'none'); % Invisible placeholder to establish axes
        title('DFF and BYV\_deg vs Time in Selected Regions')
        xlabel('Time (s)')
        ylabel('')  % Clear any axis labels
        
        % Clear right axis
        yyaxis right;
        ylabel('')  % Clear any axis labels
        
        % Reset to left axis as default
        yyaxis left;
        
        grid on
        legend off  % Remove any existing legend
        
        % Clear summary text
        set(summary_text, 'String', 'Draw polygon regions to see DFF and BYV_deg analysis');
        end
    
    function save_figure(~, ~)
        % Generate filename with timestamp and region info
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        num_regions = length(selected_points);
        
        if num_regions == 0
            filename = sprintf('trajectory_analysis_%s_no_regions', timestamp);
        else
            total_spatial_points = sum(cellfun(@(x) x.spatial_count, selected_points));
            filename = sprintf('trajectory_analysis_%s_%d_regions_%d_spatial_points', ...
                timestamp, num_regions, total_spatial_points);
        end
        
        % Ask user for save location and format
        [file, path, filterindex] = uiputfile({...
            '*.png', 'PNG files (*.png)';...
            '*.fig', 'MATLAB Figure files (*.fig)';...
            '*.pdf', 'PDF files (*.pdf)';...
            '*.eps', 'EPS files (*.eps)';...
            '*.tiff', 'TIFF files (*.tiff)';...
            '*.svg', 'SVG files (*.svg)'}, ...
            'Save Figure As', fullfile(p.Results.save_path, filename));
        
        if file ~= 0  % User didn't cancel
            full_filename = fullfile(path, file);
            
            % Get file extension to determine format
            [~, ~, ext] = fileparts(file);
            
            try
                switch lower(ext)
                    case '.fig'
                        savefig(fig, full_filename);
                    case '.png'
                        saveas(fig, full_filename, 'png');
                    case '.pdf'
                        saveas(fig, full_filename, 'pdf');
                    case '.eps'
                        saveas(fig, full_filename, 'epsc');
                    case '.tiff'
                        saveas(fig, full_filename, 'tiff');
                    case '.svg'
                        saveas(fig, full_filename, 'svg');
                    otherwise
                        saveas(fig, full_filename, 'png');  % Default to PNG
                end
                
                % Show success message
                total_spatial = sum(cellfun(@(x) x.spatial_count, selected_points));
                msgbox(sprintf('Figure saved successfully!\n\nFile: %s\nRegions: %d\nSpatial points: %d', ...
                    full_filename, num_regions, total_spatial), 'Save Complete');
                
                fprintf('Figure saved: %s\n', full_filename);
                
            catch ME
                % Show error message
                msgbox(sprintf('Error saving figure:\n%s', ME.message), 'Save Error', 'error');
            end
        end
    end
    
    function update_signals_plot()
        axes(ax_signals);
        cla;
        
        if isempty(selected_points)
            plot(0, 0); % Placeholder
            title('DFF and BYV\_deg vs Time in Selected Regions')
            xlabel('Time (s)')
        else
            % Set up dual y-axes
            yyaxis left
            hold on;
            
            legend_entries = {};
            dff_handles = [];
            
            % Plot DFF on left axis
            for i = 1:length(selected_points)
                if selected_points{i}.dff_count > 0
                    base_color = regions{i}.color;
                    
                    % Sort by time for proper line plotting (DFF)
                    [sorted_time_dff, sort_idx_dff] = sort(selected_points{i}.time_dff);
                    sorted_dff = selected_points{i}.dff(sort_idx_dff);
                    
                    h_dff = plot(sorted_time_dff, sorted_dff, ...
                        'o-', 'Color', base_color, 'LineWidth', 2, ...
                        'MarkerSize', 5, 'MarkerFaceColor', base_color);
                    
                    dff_handles = [dff_handles, h_dff];
                    legend_entries{end+1} = sprintf('Region %d DFF (n=%d)', i, selected_points{i}.dff_count);
                end
            end
            
            ylabel('DFF', 'Color', 'b')
            ax_left = gca;
            ax_left.YColor = 'b';
            
            % Plot BYV_deg on right axis
            yyaxis right
            hold on;
            
            byv_handles = [];
            
            for i = 1:length(selected_points)
                if selected_points{i}.supp_count > 0
                    base_color = regions{i}.color;
                    
                    % Get darker shade for BYV
                    if base_color == 'r'
                        byv_color = [0.6, 0, 0]; % Dark red
                    elseif base_color == 'g'
                        byv_color = [0, 0.6, 0]; % Dark green
                    elseif base_color == 'm'
                        byv_color = [0.6, 0, 0.6]; % Dark magenta
                    elseif base_color == 'c'
                        byv_color = [0, 0.6, 0.6]; % Dark cyan
                    elseif base_color == 'y'
                        byv_color = [0.8, 0.6, 0]; % Dark yellow/orange
                    else % 'k'
                        byv_color = [0.3, 0.3, 0.3]; % Dark gray
                    end
                    
                    % Sort by time for proper line plotting (BYV)
                    [sorted_time_supp, sort_idx_supp] = sort(selected_points{i}.time_supp);
                    sorted_byv_deg = selected_points{i}.byv_deg(sort_idx_supp);
                    
                    h_byv = plot(sorted_time_supp, sorted_byv_deg, ...
                        's--', 'Color', byv_color, 'LineWidth', 1.5, ...
                        'MarkerSize', 4, 'MarkerFaceColor', byv_color);
                    
                    byv_handles = [byv_handles, h_byv];
                    legend_entries{end+1} = sprintf('Region %d BYV_deg (n=%d)', i, selected_points{i}.supp_count);
                end
            end
            
            ylabel('BYV\_deg (°)', 'Color', 'r')
            ax_right = gca;
            ax_right.YColor = 'r';
            
            title('DFF (left, circles) and BYV\_deg (right, squares) vs Time')
            xlabel('Time (s)')
            
            % Create legend with both types of handles
            if ~isempty(legend_entries)
                all_handles = [dff_handles, byv_handles];
                legend(all_handles, legend_entries, 'Location', 'best', 'FontSize', 9);
            end
        end
        
        grid on;
    end
    
    function update_summary_display()
        if isempty(selected_points)
            summary_str = 'Draw polygon regions to see DFF and BYV_deg analysis';
        else
            summary_str = sprintf('REGION SUMMARY:\n\n');
            
            for i = 1:length(selected_points)
                summary_str = [summary_str, sprintf('Region %d (%.2f-%.2fs):\n', i, ...
                    selected_points{i}.time_range(1), selected_points{i}.time_range(2))];
                summary_str = [summary_str, sprintf('  Spatial points: %d\n', selected_points{i}.spatial_count)];
                summary_str = [summary_str, sprintf('  DFF points: %d\n', selected_points{i}.dff_count)];
                summary_str = [summary_str, sprintf('  BYV points: %d\n', selected_points{i}.supp_count)];
                
                if selected_points{i}.dff_count > 0
                    summary_str = [summary_str, sprintf('  DFF: %.3f ± %.3f\n', ...
                        mean(selected_points{i}.dff), std(selected_points{i}.dff))];
                end
                
                if selected_points{i}.supp_count > 0
                    summary_str = [summary_str, sprintf('  BYV_deg: %.1f ± %.1f°\n', ...
                        mean(selected_points{i}.byv_deg), std(selected_points{i}.byv_deg))];
                end
                
                summary_str = [summary_str, sprintf('\n')];
            end
        end
        
        set(summary_text, 'String', summary_str);
    end
    
    function export_results(~, ~)
        if isempty(selected_points)
            msgbox('No regions selected to export!', 'Warning');
            return;
        end
        
        % Create results structure
        results = struct();
        results.regions = regions;
        results.selected_points = selected_points;
        results.total_regions = length(selected_points);
        results.summary_stats = struct();
        
        % Calculate overall statistics
        all_dff = [];
        all_byv_deg = [];
        for i = 1:length(selected_points)
            if selected_points{i}.dff_count > 0
                all_dff = [all_dff; selected_points{i}.dff];
            end
            if selected_points{i}.supp_count > 0
                all_byv_deg = [all_byv_deg; selected_points{i}.byv_deg];
            end
            
            % Individual region stats
            results.summary_stats.(sprintf('region_%d', i)) = struct(...
                'spatial_count', selected_points{i}.spatial_count, ...
                'dff_count', selected_points{i}.dff_count, ...
                'supp_count', selected_points{i}.supp_count, ...
                'time_range', selected_points{i}.time_range);
            
            if selected_points{i}.dff_count > 0
                results.summary_stats.(sprintf('region_%d', i)).mean_dff = mean(selected_points{i}.dff);
                results.summary_stats.(sprintf('region_%d', i)).std_dff = std(selected_points{i}.dff);
            end
            
            if selected_points{i}.supp_count > 0
                results.summary_stats.(sprintf('region_%d', i)).mean_byv_deg = mean(selected_points{i}.byv_deg);
                results.summary_stats.(sprintf('region_%d', i)).std_byv_deg = std(selected_points{i}.byv_deg);
            end
        end
        
        if ~isempty(all_dff) && ~isempty(all_byv_deg)
            results.summary_stats.overall = struct(...
                'mean_dff', mean(all_dff), ...
                'std_dff', std(all_dff), ...
                'mean_byv_deg', mean(all_byv_deg), ...
                'std_byv_deg', std(all_byv_deg), ...
                'total_dff_points', length(all_dff), ...
                'total_byv_points', length(all_byv_deg));
        end
        
        % Export to base workspace
        assignin('base', 'region_selection_results', results);
        
        % Also create individual variables for each region
        for i = 1:length(selected_points)
            var_name = sprintf('region_%d_data', i);
            assignin('base', var_name, selected_points{i});
        end
        
        % Show summary
        total_spatial = sum(cellfun(@(x) x.spatial_count, selected_points));
        msgbox(sprintf('Results exported!\n\nRegions: %d\nSpatial points: %d\nDFF points: %d\nBYV points: %d', ...
            length(selected_points), total_spatial, length(all_dff), length(all_byv_deg)), ...
            'Export Complete');
        
        % Print to command window
        fprintf('\n=== REGION SELECTION RESULTS (TIME-ALIGNED) ===\n');
        for i = 1:length(selected_points)
            fprintf('Region %d (%.2f-%.2fs):\n', i, ...
                selected_points{i}.time_range(1), selected_points{i}.time_range(2));
            fprintf('  Spatial: %d points\n', selected_points{i}.spatial_count);
            fprintf('  DFF: %d points', selected_points{i}.dff_count);
            if selected_points{i}.dff_count > 0
                fprintf(', mean = %.3f ± %.3f', mean(selected_points{i}.dff), std(selected_points{i}.dff));
            end
            fprintf('\n  BYV: %d points', selected_points{i}.supp_count);
            if selected_points{i}.supp_count > 0
                fprintf(', mean = %.1f ± %.1f°', mean(selected_points{i}.byv_deg), std(selected_points{i}.byv_deg));
            end
            fprintf('\n');
        end
    end
end