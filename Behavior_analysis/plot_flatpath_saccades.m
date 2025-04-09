%% a function for plotting the flat path with detected saccades

function plot_flatpath_saccades(daq, jump, savepath, upsampled_x_data, upsampled_y_data, turn)

if isfield(daq, 'optoStim')
    % downsample
    downsampled_exptData = struct();
    fields = fieldnames(daq);
    for i = 1:length(fields)
        field = fields{i};
    
        % Check if the field is numeric
        if isnumeric(daq.(field))
            % Perform downsampling by taking every 100th element
            downsampled_exptData.(field) = daq.(field)(1:50:end);
        end
    end
    opto=1;
    x = downsampled_exptData.x;  % x-coordinates
    y = downsampled_exptData.y ;  % y-coordinates
    opto_timing = downsampled_exptData.optoStim;
    saccading = downsampled_exptData.saccading;
else
    opto=0;
    x = upsampled_x_data; %daq.px  ;  % x-coordinates
    y = upsampled_y_data; %daq.py;  % y-coordinates
    if turn
        saccading=daq.turning_supp;
    else
        saccading = daq.saccading_supp;
    end
end

% Create a figure for plotting
figure; 
hold on;

% Initialize handles for legend
hSaccading = [];  % For storing saccading graphic handles
hNormal = [];     % For storing normal graphic handles
hOptoonly = [];
hSaccading_opto =[];
hJump = [];       % For storing jump graphic handles

total=0;
% Plot all segments
for i = 1:length(x)-1
    if opto  % If the current index is part of a saccade
        if opto_timing(i) ~= 0  % Check if optostim is on
            if saccading(i)
                % Highlight saccade with a pink line when opto is on
                hSaccading_opto = plot(x(i:i+1), y(i:i+1), 'm', 'LineWidth', 2);  % Pink for saccading with opto
            else
                % Highlight saccade with a red line when opto is off
                hOptoonly = plot(x(i:i+1), y(i:i+1), 'g', 'LineWidth', 2);  % green for opto
            end
        else
            if saccading(i)
                % Highlight saccade with a pink line when opto is on
                hSaccading = plot(x(i:i+1), y(i:i+1), 'r', 'LineWidth', 2);  % red for saccading
            else
                % Highlight saccade with a red line when opto is off
                hNormal = plot(x(i:i+1), y(i:i+1), 'b', 'LineWidth', 2);  % blue for nothing
            end
        end
    else
        if saccading(i)
            hSaccading = plot(x(i:i+1), y(i:i+1), 'r', 'LineWidth', 2);  % Red for saccading
            total = total+1;
        else
            % Plot non-saccading segments in blue
            hNormal = plot(x(i:i+1), y(i:i+1), 'b', 'LineWidth', 1);  % Blue for non-saccading
        end
    end
end

% Plot the start and end points
hStart = plot(x(1), y(1), 'Color', [1, 0.5, 0], 'Marker', '.', 'MarkerSize', 25);  % Start point (green dot)
hEnd = plot(x(end), y(end), 'y.', 'MarkerSize', 25);  % End point (yellow dot)

% If jump is detected, plot the corresponding markers
if jump
    if opto
        jump_indices = find(daq.jump_detected == 1);  % Find indices where jump detected is 1
        hJump = plot(daq.x(jump_indices), daq.y(jump_indices), 'kx', 'MarkerSize', 10);  % Black 'X' markers for jumps
    else
        jump_indices = find(daq.vis.jump_detected == 1);  % Find indices where jump detected is 1
        hJump = plot(x(jump_indices), y(jump_indices), 'kx', 'MarkerSize', 10);  % Black 'X' markers for jumps
    end
end

% Add labels and title
xlabel('X Position');
ylabel('Y Position');
title('Fly Trajectory with Saccades Highlighted');

% Create the legend using the handles declared above
if jump
    if opto
        legend([hSaccading, hSaccading_opto, hOptoonly, hNormal, hStart, hEnd, hJump], {'Saccading', 'Saccading with Opto', 'Opto On', 'Normal', 'Start', 'End', 'Jump'}, 'Location', 'Best');
    else
        legend([hSaccading, hNormal, hStart, hEnd, hJump], {'Saccading', 'Normal', 'Start', 'End', 'Jump'}, 'Location', 'Best');
    end
else
    if opto
        legend([hSaccading, hOptoonly, hNormal, hStart, hEnd], {'Saccading with Opto', 'Opto On', 'Normal', 'Start', 'End'}, 'Location', 'Best');
    else
        legend([hSaccading, hNormal, hStart, hEnd], {'Saccading', 'Normal', 'Start', 'End'}, 'Location', 'Best');
    end
end
grid on;
hold off;


save_plot_with_title_as_filename('flatpath', 'colors', savepath)
end