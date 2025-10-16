function [daq, dat, dff_motion] = openloop_gratings_process(daq, dat, savepath, box, splitgrating)
%% calculate normalized dff if exists 2 rois
    a=size(dat.ts);
    if a(1)>1
        dff = dat.ts(1,:)-dat.ts(2,:) ;
    else
        dff = dat.ts(1,:);
    end

%% --- LABEL SWITCH ---
if splitgrating
    label_cw  = 'Back-to-Front';
    label_ccw = 'Front-to-Back';
else
    label_cw  = 'Clockwise';
    label_ccw = 'Counterclockwise';
end

%% Plot the dff and pattern
    % Plot the dff and fwd
    figure;
    yyaxis left;
    plot(daq.t_supp, daq.vh_supp, '-k');
    ylabel('pattern');  % Label for the left y-axis
    
    
    % Plot the second time series on the right y-axis
    yyaxis right;
    plot(daq.t, dff, '-b');  % Plot second time series in red
    ylabel('dff');  % Label for the right y-axis
    
    % Add title and x-axis label
    xlabel('Time');
    title('Pattern and DFF');
    
    % Optional: Add legend
    legend('Pattern', 'dff');
    save_plot_with_title_as_filename('Pattern', 'dff', savepath);
    %% variables
    % frequency
    freq = daq.t(end)/length(daq.t);
    % Define the number of timepoints to include before and after motion
    len_sec = 2; % length in seconds of time to buffer with
    stimlen = 2; % how long the pattern moves for
    pre_motion_points = round(len_sec/freq);
    post_motion_points = round(len_sec/freq);
    remove_if_shorter_than = round(1.5/freq);

   
    %value for total speed (moving or no)
    still_cutoff = 1;
    
    %% using velocity to determine moving clockwise or counter
    % for 3 Hz

    
    positive_slope = zeros(size(daq.vh));
    positive_slope(daq.vh >= -3.1 & daq.vh <= 0) = 1;
    % take care of 6hz case
    if max(daq.vvy)>2.5
         positive_slope(daq.vh >= -3.11 & daq.vvy >= 0) = 1;
    end
    diff_array = diff([0, positive_slope, 0]);
    run_starts = find(diff_array == 1);
    run_ends = find(diff_array == -1) - 1;
    
    for i = 1:length(run_starts)
        run_length = run_ends(i) - run_starts(i) + 1;
        if run_length < remove_if_shorter_than
            positive_slope(run_starts(i):run_ends(i)) = 0;
        end
    end
    % ccw
    negative_slope = zeros(size(daq.vh));
    negative_slope(daq.vh >= 0 & daq.vh <= 4) = 1;
    if max(daq.vvy)>2.5
         negative_slope(daq.vh >= -3.11 & daq.vvy <= 0) = 1;
    end
    diff_array = diff([0, negative_slope, 0]);
    run_starts = find(diff_array == 1);
    run_ends = find(diff_array == -1) - 1;
    
    for i = 1:length(run_starts)
        run_length = run_ends(i) - run_starts(i) + 1;
        if run_length < remove_if_shorter_than
            negative_slope(run_starts(i):run_ends(i)) = 0;
        end
    end

    % now ensure there are no tiny gaps
    vec = positive_slope; % Your original 1D row vector (1 x N)
    vec_padded = [1, vec, 1]; % Pad with ones at both ends to allow flanked checking at edges
    out = vec;
    
    z_idx = find(vec_padded == 0); % Indices of zeros (in the padded array)
    
    for i = 1:length(z_idx)
        idx = z_idx(i);
        % Check: if previous and next are 1s (in the padded array)
        if vec_padded(idx-1)==1 && vec_padded(idx+1)==1
            % It's a singleton, copy to original array
            out(idx-1)=1;
        elseif vec_padded(idx-1)==1 && vec_padded(idx+2)==1 && vec_padded(idx+1)==0
            % It's doublet (two zeros), check both are zero and both side are 1s
            out(idx-1:idx)=1;
        end
    end
    positive_slope = out;
    vec = negative_slope; % Your original 1D row vector (1 x N)
    vec_padded = [1, vec, 1]; % Pad with ones at both ends to allow flanked checking at edges
    out = vec;
    
    z_idx = find(vec_padded == 0); % Indices of zeros (in the padded array)
    
    for i = 1:length(z_idx)
        idx = z_idx(i);
        % Check: if previous and next are 1s (in the padded array)
        if vec_padded(idx-1)==1 && vec_padded(idx+1)==1
            % It's a singleton, copy to original array
            out(idx-1)=1;
        elseif vec_padded(idx-1)==1 && vec_padded(idx+2)==1 && vec_padded(idx+1)==0
            % It's doublet (two zeros), check both are zero and both side are 1s
            out(idx-1:idx)=1;
        end
    end
    negative_slope = out;


    %% now chop up with preface & plot - CLOCKWISE

    % start with clockwise motion
    % Initialize cell arrays to store segments of velocity
    dff_cw_motion = {};  
    pattern_cw_motion = {}; 
    total_speed_cw_motion ={};
    forward_speed_cw_motion ={};
    rotational_speed_cw_motion ={};
    times_cw_motion ={};

    motion_on_markers =[];
    motion_off_markers =[];

    % find changes of onset of cw motion. this is opposite between bars and
    % boxes
    if box
        change_points = [1, find(diff(positive_slope) ~= 0) + 1, length(positive_slope)];
        clockwise_motion = positive_slope;
    else
        change_points = [1, find(diff(negative_slope) ~= 0) + 1, length(negative_slope)];
        clockwise_motion = negative_slope;
    end

    % out is the modified binary vector with isolated/sandwiched zeros filled with 1s
    % Loop through each segment and separate based on motion value
    % skip first and last presentation due to issues with timing
    % start at 3 bc 1 is just 1 (start of trial) and first presentation
    % doesnt have enough pre-data

    %determine how long to cut at the end
    % pi/2 = 3
    % pi = 2
    for i = 3:length(change_points) - 3
        % Get the indices for the current segment
        segment_indices = change_points(i):change_points(i+1) - 1;
    
        % Check the motion value in this segment
        if clockwise_motion(segment_indices(1)) > 0
            % Calculate the start index for pre-motion points (ensuring it doesn’t go below 1)
            pre_motion_point = max(1, segment_indices(1) - pre_motion_points);
            % Calculate the end index for post-motion points (ensuring it doesn’t go below 1)
            post_motion_point = min(length(clockwise_motion), segment_indices(end) + post_motion_points);
            %select entire peirod
            extended_segment_indices = pre_motion_point:post_motion_point;

            % only motion on fwd and rot
            derivative_of_position = diff(daq.vh);
            dff_cw_motion{end+1} = dff(extended_segment_indices);
            pattern_cw_motion{end+1} = daq.vh(extended_segment_indices);

            % changed to totalspeed_supp
            upsample_factor = length(daq.t_supp)/length(daq.t);
            extended_segment_indices_supp = round(pre_motion_point*upsample_factor):round(post_motion_point*upsample_factor);
            total_speed_cw_motion{end+1} = daq.totalspeed_supp(extended_segment_indices_supp);
            forward_speed_cw_motion{end+1} = daq.bvf(extended_segment_indices);
            rotational_speed_cw_motion{end+1} = abs(daq.bvy_supp(extended_segment_indices_supp));
            times_cw_motion{end+1} = daq.t(extended_segment_indices);
            
            % Record the exact point where motion switches on for plotting
            motion_on_markers(end+1) = segment_indices(1) - pre_motion_point + 1;
            motion_off_markers(end+1) = segment_indices(end)-segment_indices(1)+pre_motion_point;% Relative position within segment
        end
    end

    % Ensure all "motionon" segments are the same length by trimming to the shortest segment
    min_length = min(cellfun(@length, dff_cw_motion));
    dff_cw_motion = cellfun(@(x) x(1:min_length), dff_cw_motion, 'UniformOutput', false);
    pattern_cw_motion = cellfun(@(x) x(1:min_length), pattern_cw_motion, 'UniformOutput', false);
    min_length_speed = min(cellfun(@length, total_speed_cw_motion));
    total_speed_cw_motion = cellfun(@(x) x(1:min_length_speed), total_speed_cw_motion, 'UniformOutput', false);
    forward_speed_cw_motion = cellfun(@(x) x(1:min_length), forward_speed_cw_motion, 'UniformOutput', false);
    rotational_speed_cw_motion = cellfun(@(x) x(1:min_length_speed), rotational_speed_cw_motion, 'UniformOutput', false);
    times_cw_motion= cellfun(@(x) x(1:min_length), times_cw_motion, 'UniformOutput', false);
    
    % Convert to a matrix for calculating the average dff
    dff_ave = cell2mat(reshape(dff_cw_motion, length(dff_cw_motion), []));
    average_dff = mean(dff_ave, 1);

    % Normalize each dff cell so values are between 0 and 1
    dff_cw_motion_norm = cellfun(@(x) (x - min(x)) / (max(x) - min(x)), dff_cw_motion, 'UniformOutput', false);
    cwdff_ave_norm = cell2mat(reshape(dff_cw_motion_norm, length(dff_cw_motion_norm), []));
    average_cw_dff_norm = mean(cwdff_ave_norm, 1);

    % Plot the "motionon" segments with pre-motion points and vertical lines
    figure;
    
    % Subplot 1: Velocity
    subplot(3, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(dff_cw_motion)
        plot(daq.t(1:length(dff_cw_motion{i})), dff_cw_motion{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t(motion_on_markers(1))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_cw_motion{i})), average_dff, 'g', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    % Customize the first subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title([label_cw ' DFF']);
    legend(label_cw);
    hold off;

    % Subplot 2: norm dff
    subplot(3, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot each raw "motionon" speed segment in lighter blue
    for i = 1:length(dff_cw_motion_norm)
        plot(daq.t(1:length(dff_cw_motion_norm{i})), dff_cw_motion_norm{i}, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_cw_motion{i})), average_cw_dff_norm, 'black', 'LineWidth', 3, 'DisplayName', 'Average DFF');

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('DFF (norm)');
    title('DFF Normalized');
    legend('Average motion-On Angular Speed');
    hold off;
    
    % Subplot 3: pattern
    subplot(3, 1, 3);  
    hold on;
    % Plot each raw "motionon" speed segment in lighter blue
    for i = 1:length(pattern_cw_motion)
        plot(daq.t(1:length(pattern_cw_motion{i})), pattern_cw_motion{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('pattern position');
    title('pattern');
    legend('Average motion-On Angular Speed');
    hold off;
    save_plot_with_title_as_filename('pattern', label_cw, savepath)

    %% ccw
    % start with clockwise motion
    % Initialize cell arrays to store segments of velocity
    dff_ccw_motion = {};  
    pattern_ccw_motion = {}; 
    total_speed_ccw_motion ={};
    forward_speed_ccw_motion ={};
    rotational_speed_ccw_motion ={};
    times_ccw_motion ={};
    
    ccw_motion_on_markers =[];
    ccw_motion_off_markers =[];


    % find changes of onset of ccw motion - box and everything else are
    % opposite
    if box
        change_points = [1, find(diff(negative_slope) ~= 0) + 1, length(negative_slope)];
        counterclockwise_motion = negative_slope;
    else
        change_points = [1, find(diff(positive_slope) ~= 0) + 1, length(positive_slope)];
        counterclockwise_motion = positive_slope;
    end


    % Loop through each segment and separate based on motion value
    % skip first and last presentation due to issues with timing
    for i = 3:length(change_points) - 3
        % Get the indices for the current segment
        segment_indices = change_points(i):change_points(i+1) - 1;
    
        % Check the motion value in this segment
        if counterclockwise_motion(segment_indices(1)) > 0
            % Calculate the start index for pre-motion points (ensuring it doesn’t go below 1)
            pre_motion_point = max(1, segment_indices(1) - pre_motion_points);
            % Calculate the end index for post-motion points (ensuring it doesn’t go below 1)
            post_motion_point = min(length(counterclockwise_motion), segment_indices(end) + post_motion_points);

            %select entire peirod
            extended_segment_indices = pre_motion_point:post_motion_point;

            % only motion on fwd and rot
            derivative_of_position = diff(daq.vh);
            dff_ccw_motion{end+1} = dff(extended_segment_indices);
            pattern_ccw_motion{end+1} = daq.vh(extended_segment_indices);
            
            % changed to totalspeed_supp
            upsample_factor = length(daq.t_supp)/length(daq.t);
            extended_segment_indices_supp = round(pre_motion_point*upsample_factor):round(post_motion_point*upsample_factor);
            total_speed_ccw_motion{end+1} = daq.totalspeed_supp(extended_segment_indices_supp);
            forward_speed_ccw_motion{end+1} = daq.bvf(extended_segment_indices);
            rotational_speed_ccw_motion{end+1} = abs(daq.bvy_supp(extended_segment_indices_supp));
            times_ccw_motion{end+1} = daq.t(extended_segment_indices);
            
            % Record the exact point where motion switches on for plotting
            ccw_motion_on_markers(end+1) = segment_indices(1) - pre_motion_point + 1;
            ccw_motion_off_markers(end+1) = segment_indices(end)-segment_indices(1)+pre_motion_point;% Relative position within segment
        end
    end

    % Ensure all "motionon" segments are the same length by trimming to the shortest segment
    min_length = min(cellfun(@length, dff_ccw_motion));
    dff_ccw_motion = cellfun(@(x) x(1:min_length), dff_ccw_motion, 'UniformOutput', false);
    pattern_ccw_motion = cellfun(@(x) x(1:min_length), pattern_ccw_motion, 'UniformOutput', false);
    min_length_speed = min(cellfun(@length, total_speed_ccw_motion));
    total_speed_ccw_motion = cellfun(@(x) x(1:min_length_speed), total_speed_ccw_motion, 'UniformOutput', false);
    forward_speed_ccw_motion = cellfun(@(x) x(1:min_length), forward_speed_ccw_motion, 'UniformOutput', false);
    rotational_speed_ccw_motion = cellfun(@(x) x(1:min_length_speed), rotational_speed_ccw_motion, 'UniformOutput', false);
    times_ccw_motion= cellfun(@(x) x(1:min_length), times_ccw_motion, 'UniformOutput', false);
    
    % Convert to a matrix for calculating the average dff
    ccw_dff_ave = cell2mat(reshape(dff_ccw_motion, length(dff_ccw_motion), []));
    CCW_average_dff = mean(ccw_dff_ave, 1);

    % Normalize each dff cell so values are between 0 and 1
    dff_ccw_motion_norm = cellfun(@(x) (x - min(x)) / (max(x) - min(x)), dff_ccw_motion, 'UniformOutput', false);
    ccwdff_ave_norm = cell2mat(reshape(dff_ccw_motion_norm, length(dff_ccw_motion_norm), []));
    average_ccw_dff_norm = mean(ccwdff_ave_norm, 1);

    % Plot the "motionon" segments with pre-motion points and vertical lines
    figure;
    
    % Subplot 1: dff ccw 
    subplot(3, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(dff_ccw_motion)
        plot(daq.t(1:length(dff_ccw_motion{i})), dff_ccw_motion{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_ccw_motion{i})), CCW_average_dff, 'g', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    % Customize the first subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_ccw);
    legend(label_ccw);
    hold off;

    % Subplot 2: norm dff ccw
    subplot(3, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot each raw "motionon" speed segment in lighter blue
    for i = 1:length(dff_ccw_motion_norm)
        plot(daq.t(1:length(dff_ccw_motion_norm{i})), dff_ccw_motion_norm{i}, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_ccw_motion{i})), average_ccw_dff_norm, 'black', 'LineWidth', 3, 'DisplayName', 'Average DFF');

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('DFF (norm)');
    title('DFF Normalized');
    legend('Average motion-On Angular Speed');
    hold off;
    
    % Subplot 3: pattern
    subplot(3, 1, 3);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot each raw "motionon" speed segment in lighter blue
    for i = 1:length(pattern_ccw_motion)
        plot(daq.t(1:length(pattern_ccw_motion{i})), pattern_ccw_motion{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('pattern position');
    title('pattern');
    legend(label_ccw);
    hold off;
    save_plot_with_title_as_filename('pattern', label_ccw, savepath)

    %% plot cw and ccw on one chart
    % Plot the "motionon" segments with pre-motion points and vertical lines
    figure;
    
    % Subplot 1: CCW motion dff
    subplot(3, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(dff_ccw_motion)
        plot(daq.t(1:length(dff_ccw_motion{i})), dff_ccw_motion{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_ccw_motion{i})), CCW_average_dff, 'g', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    % Customize the first subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_ccw);
    legend([label_ccw ' DFF']);
    hold off;
    
    % Subplot 2: CW motion dff
    subplot(3, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot cw rotation
    for i = 1:length(dff_cw_motion)
        plot(daq.t(1:length(dff_cw_motion{i})), dff_cw_motion{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Plot the average velocity for cw in blue
    plot(daq.t(1:length(dff_cw_motion{i})), average_dff, 'b', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_cw);
    legend([label_cw ' DFF']);
    hold off;
    
    % Subplot 3: fly speed
    subplot(3, 1, 3);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(total_speed_ccw_motion)
        plot(daq.t_supp(1:length(total_speed_ccw_motion{i})), total_speed_ccw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(total_speed_cw_motion{i})), total_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    for i = 1:length(total_speed_cw_motion)
        plot(daq.t_supp(1:length(total_speed_cw_motion{i})), total_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(total_speed_cw_motion{i})), total_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    xline(len_sec, 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
    xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);

    total_speed = [total_speed_ccw_motion, total_speed_cw_motion];
    min_length = min(cellfun(@length, total_speed));
    total_speed = cellfun(@(x) x(1:min_length), total_speed, 'UniformOutput', false);
    ave_speed = cell2mat(reshape(total_speed, length(total_speed), []));
    ave_spd = mean(ave_speed, 1);
    plot(daq.t_supp(1:size(ave_speed, 2)), ave_spd, 'm', 'LineWidth', 3, 'DisplayName', 'Average Speed');

    % Customize the plot
    xlabel('Time (s)');
    ylabel('Total speed (mm/s)');
    title('Total Speed');
    save_plot_with_title_as_filename(label_cw, [label_ccw 'andtotspeed'], savepath)
    
    %% plot cw and ccw dff with fwd and rot speed
    % Plot the "motionon" segments with pre-motion points and vertical lines
    figure;
    
    % Subplot 1: CCW motion dff
    subplot(4, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(dff_ccw_motion)
        plot(daq.t(1:length(dff_ccw_motion{i})), dff_ccw_motion{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    % Plot the average velocity for "motionon" in bold green
    plot(daq.t(1:length(dff_ccw_motion{i})), CCW_average_dff, 'g', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    % Customize the first subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_ccw);
    legend([label_ccw ' DFF']);
    hold off;
    
    % Subplot 2: CW motion dff
    subplot(4, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot cw rotation
    for i = 1:length(dff_cw_motion)
        plot(daq.t(1:length(dff_cw_motion{i})), dff_cw_motion{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
        xline(daq.t((motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end

    % Plot the average velocity for cw in blue
    plot(daq.t(1:length(dff_cw_motion{i})), average_dff, 'b', 'LineWidth', 3, 'DisplayName', 'Average DFF');
    

    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_cw);
    legend([label_cw ' DFF']);
    hold off;
    
    % Subplot 3: fly speed
    subplot(4, 1, 3);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(forward_speed_ccw_motion)
        plot(daq.t(1:length(forward_speed_ccw_motion{i})), forward_speed_ccw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(forward_speed_cw_motion{i})), forward_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    for i = 1:length(forward_speed_cw_motion)
        plot(daq.t(1:length(forward_speed_cw_motion{i})), forward_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(forward_speed_cw_motion{i})), forward_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    xline(len_sec, 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
    xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);

    forward_speed = [forward_speed_ccw_motion, forward_speed_cw_motion];
    min_length = min(cellfun(@length, forward_speed));
    forward_speed = cellfun(@(x) x(1:min_length), forward_speed, 'UniformOutput', false);
    ave_speed = cell2mat(reshape(forward_speed, length(forward_speed), []));
    ave_spd = mean(ave_speed, 1);
    plot(daq.t_supp(1:size(ave_speed, 2)), ave_spd, 'm', 'LineWidth', 3, 'DisplayName', 'Average Speed');

    % Customize the plot
    xlabel('Time (s)');
    ylabel('Forward speed (mm/s)');
    title('forward Speed');

    % Subplot 4: fly rotational speed
    subplot(4, 1, 4);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "motionon" velocity segment in lighter green
    for i = 1:length(rotational_speed_ccw_motion)
        plot(daq.t_supp(1:length(rotational_speed_ccw_motion{i})), rotational_speed_ccw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(rotational_speed_cw_motion{i})), rotational_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    for i = 1:length(rotational_speed_cw_motion)
        plot(daq.t_supp(1:length(rotational_speed_cw_motion{i})), rotational_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green
        %plot(daq.t(1:length(rotational_speed_cw_motion{i})), rotational_speed_cw_motion{i}, 'Color', [0.8 0 1], 'LineWidth', 1);  % Lighter green

    end
    xline(len_sec, 'k--', 'LineWidth', 1.5);  % Dashed black line at motion-on marker
    xline(len_sec + stimlen, 'k--', 'LineWidth', 1.5);

    rotational_speed = [rotational_speed_ccw_motion, rotational_speed_cw_motion];
    min_length = min(cellfun(@length, rotational_speed));
    rotational_speed = cellfun(@(x) x(1:min_length), rotational_speed, 'UniformOutput', false);
    ave_speed = cell2mat(reshape(rotational_speed, length(rotational_speed), []));
    ave_spd = mean(ave_speed, 1);
    plot(daq.t_supp(1:size(ave_speed, 2)), ave_spd, 'm', 'LineWidth', 3, 'DisplayName', 'Average Speed');

    % Customize the plot
    xlabel('Time (s)');
    ylabel('rotational speed (mm/s)');
    title('rotational Speed');
    save_plot_with_title_as_filename([label_cw label_ccw], 'fwd_rot', savepath)



    %% just cw and ccw dff on one axis
    figure
    hold on;
    
    % Plot each raw "motionon" velocity segment in lighter green
    h1 = [];  % To store handle for legend
    for i = 1:length(dff_ccw_motion)
        h = plot(daq.t(1:length(dff_ccw_motion{i})), dff_ccw_motion{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);
        if i == 1  % Only keep first handle for legend
            h1 = h;
        end
    end
    
    % Plot each raw "motionon" velocity segment in lighter blue
    h2 = [];  % To store handle for legend
    for i = 1:length(dff_cw_motion)
        h = plot(daq.t(1:length(dff_cw_motion{i})), dff_cw_motion{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);
        if i == 1  % Only keep first handle for legend
            h2 = h;
        end
    end

    length_cut = min(length(average_dff), length(CCW_average_dff));
    
    bidirectional_dff = (average_dff(1:length_cut)+CCW_average_dff(1:length_cut))/2;
    h3 = plot(daq.t(1:length_cut), bidirectional_dff, 'm', 'LineWidth', 3);
    
    % Create legend with specific handles and labels
    legend([h1, h2, h3], {label_ccw, label_cw, 'Average'}, 'Location', 'best');
    
    xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
    xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    
    % Customize the plot
    xlabel('Time (s)');
    ylabel('DFF');
    title('Bi-Directional Visual Rotation');
    save_plot_with_title_as_filename(label_cw, [label_cw '_samechart'], savepath)


    %% now separate walking vs not - cw

    % first merge cw and ccw
    dff_total= dff_cw_motion;
    speed= forward_speed_cw_motion;
    still = {}; 
    moving={};
    ave_still = {}; 
    ave_moving ={};
    
    for i = 1:length(dff_total)
        pre_speed = mean(speed{i}(pre_motion_points:length(speed{i})-post_motion_points));
        pre_speed
        if pre_speed <still_cutoff
            disp("Still")
            still{end+1} = dff_total{i};
            ave_still{end+1} = dff_total{i};  
        else
            moving{end+1} = dff_total{i};
            ave_moving{end+1} = dff_total{i};
            disp("moving")
        end
    end

    % take average fwd
    % Convert to a matrix for calculating the average
    min_length_moving = min(cellfun(@length, ave_moving));
    min_length_still= min(cellfun(@length, ave_still));
    min_length = min(min_length_moving, min_length_still);
    ave_moving = cellfun(@(x) x(1:min_length), ave_moving, 'UniformOutput', false);
    ave_still = cellfun(@(x) x(1:min_length), ave_still, 'UniformOutput', false);
    moving_ave = cell2mat(reshape(ave_moving, length(ave_moving), []));
    moving_ave = mean(moving_ave, 1);
    % Convert to a matrix for calculating the average
    still_ave = cell2mat(reshape(ave_still, length(ave_still), []));
    still_ave = mean(still_ave, 1);

    %% plot cw moving not moving
    figure;
    
    % ---- Subplot 1: moving ----
    ax(1) = subplot(2, 1, 1);  
    hold on;
    for i = 1:length(moving)
        plot(daq.t(1:length(moving{i})), moving{i}, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
    end
    xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
    xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    plot(daq.t(1:length(moving_ave)), moving_ave, 'black', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
    xlabel('Time (s)');
    ylabel('Moving DFF during rotation');
    title('Moving during visual panel motion');
    hold off;
    
    % ---- Subplot 2: still ----
    ax(2) = subplot(2, 1, 2);
    hold on;
    for i = 1:length(still)
        plot(daq.t(1:length(still{i})), still{i}, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
    end
    plot(daq.t(1:length(still_ave)), still_ave, 'black', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
    xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
    xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(sprintf('DFF still (still = forward <%.2f)', still_cutoff));
    title('Still during visual panel motion');
    hold off;
    
    % ---- Link axes ----
    linkaxes(ax, 'xy'); % Now both axes have the same y and x limits
    
    % ---- Add the patch to both axes using the FINAL y-limits ----
    final_ylims = ylim(ax(1)); % Or ax(2), they're the same now
    
    opto_on = daq.t(ccw_motion_on_markers(1));
    opto_off = daq.t((ccw_motion_on_markers(1)))+stimlen;
    X_patch = [opto_on, opto_off, opto_off, opto_on];
    Y_patch = [final_ylims(1), final_ylims(1), final_ylims(2), final_ylims(2)];
    
    for i = 1:2
        axes(ax(i)); % make current
        hold on;
        patch(X_patch, Y_patch, [0.5, 0, 0.5], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.2);
        hold off;
    end
    save_plot_with_title_as_filename([label_cw '_dff_moving'], [label_cw '_dff_still'], savepath)

     %% now separate walking vs not - ccw
    dff_total= dff_ccw_motion; %[dff_cw_motion, dff_ccw_motion];
    speed= forward_speed_ccw_motion; %[total_speed_cw_motion, total_speed_ccw_motion];
    still = {}; 
    moving={};
    ave_still = {}; 
    ave_moving ={};
    
    for i = 1:length(dff_total)
        pre_speed = mean(speed{i}(pre_motion_points:length(speed{i})-post_motion_points));
        if pre_speed <still_cutoff
            disp("Still")
            still{end+1} = dff_total{i};
            ave_still{end+1} = dff_total{i};  
        else
            moving{end+1} = dff_total{i};
            ave_moving{end+1} = dff_total{i};
            disp("moving")
        end
    end

    % take average fwd
    % Convert to a matrix for calculating the average
    min_length_moving = min(cellfun(@length, ave_moving));
    min_length_still= min(cellfun(@length, ave_still));
    min_length = min(min_length_moving, min_length_still);
    ave_moving = cellfun(@(x) x(1:min_length), ave_moving, 'UniformOutput', false);
    ave_still = cellfun(@(x) x(1:min_length), ave_still, 'UniformOutput', false);
    moving_ave = cell2mat(reshape(ave_moving, length(ave_moving), []));
    moving_ave = mean(moving_ave, 1);
    % Convert to a matrix for calculating the average
    still_ave = cell2mat(reshape(ave_still, length(ave_still), []));
    still_ave = mean(still_ave, 1);

    %% plot ccw moving not moving
    figure;
    
    % ---- Subplot 1: moving ----
    ax(1) = subplot(2, 1, 1);  
    hold on;
    for i = 1:length(moving)
        plot(daq.t(1:length(moving{i})), moving{i}, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
    end
    xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
    xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    plot(daq.t(1:length(moving_ave)), moving_ave, 'black', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
    xlabel('Time (s)');
    ylabel('Moving DFF during rotation');
    title('Moving during visual panel motion');
    hold off;
    
    % ---- Subplot 2: still ----
    ax(2) = subplot(2, 1, 2);
    hold on;
    for i = 1:length(still)
        plot(daq.t(1:length(still{i})), still{i}, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
    end
    plot(daq.t(1:length(still_ave)), still_ave, 'black', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
    xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
    xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(sprintf('DFF still (still = forward <%.2f)', still_cutoff));
    title('Still during visual panel motion');
    hold off;
    
    % ---- Link axes ----
    linkaxes(ax, 'xy'); % Now both axes have the same y and x limits
    
    % ---- Add the patch to both axes using the FINAL y-limits ----
    final_ylims = ylim(ax(1)); % Or ax(2), they're the same now
    
    opto_on = daq.t(ccw_motion_on_markers(1));
    opto_off = daq.t((ccw_motion_on_markers(1)))+stimlen;
    X_patch = [opto_on, opto_off, opto_off, opto_on];
    Y_patch = [final_ylims(1), final_ylims(1), final_ylims(2), final_ylims(2)];
    
    for i = 1:2
        axes(ax(i)); % make current
        hold on;
        patch(X_patch, Y_patch, [1, 0.5, 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.2);
        hold off;
    end
    save_plot_with_title_as_filename([label_ccw '_dff_moving'], [label_ccw '_dff_still'], savepath)
    %% save dff motions in struct
    % Create the combined structure
    dff_motion = struct();
    dff_motion.ccw = dff_ccw_motion;
    dff_motion.cw = dff_cw_motion;
    dff_motion.cw_ave = average_dff;
    dff_motion.ccw_ave = CCW_average_dff;
    
    % Add metadata
    dff_motion.num_ccw_trials = length(dff_ccw_motion);
    dff_motion.num_cw_trials = length(dff_cw_motion);
    dff_motion.total_trials = length(dff_ccw_motion) + length(dff_cw_motion);
    
%% now plotting the dff and kinematics with shaded visual motion
    % Create the figure and subplots
    time = daq.t;
    time_supp = daq.t_supp;
    fwd = daq.bvf_supp;
    [fwd,winsize] = smoothdata(fwd,"gaussian", 20);
    winsize
    rot = daq.bvy_deg_supp;
    figure;
    % Subplot 1: dFF
    subplot(3, 1, 1);  % First subplot
    hold on
    % Calculate your scaling
    total_padding = max(dff)-min(dff);
    clockwise_motion_fordff = (clockwise_motion*total_padding)+min(dff);
    counterclockwise_motion_fordff = (counterclockwise_motion*total_padding)+min(dff);
    
    % Find non-zero periods
    cw_mask = clockwise_motion ~= 0;
    ccw_mask = counterclockwise_motion ~= 0;
    
    % Create rectangles with scaled heights
    for i = 1:length(time)
        if cw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(dff), min(dff), clockwise_motion_fordff(i), clockwise_motion_fordff(i)], [0.5, 0, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        if ccw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(dff), min(dff), counterclockwise_motion_fordff(i), counterclockwise_motion_fordff(i)], [1, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end
    h_dff = plot(time, dff, 'b', 'LineWidth', 1.5);  % Plot dFF in red
    
    % Create invisible patches for legend
    h_cw = patch(NaN, NaN, [0.5, 0, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h_ccw = patch(NaN, NaN, [1, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    legend([h_dff, h_cw, h_ccw], {'dF/F', label_cw, label_ccw}, 'Location', 'best');
    
    xlabel('Time (s)');
    ylabel('dF/F');
    title('dF/F');
    grid on;
    
    % Subplot 2: Forward Velocity
    subplot(3, 1, 2);  % Second subplot
    hold on
    % Calculate your scaling
    total_padding = max(fwd)-min(fwd);
    clockwise_motion_fordff = (clockwise_motion*total_padding)+min(fwd);
    counterclockwise_motion_fordff = (counterclockwise_motion*total_padding)+min(fwd);
    
    % Find non-zero periods
    cw_mask = clockwise_motion ~= 0;
    ccw_mask = counterclockwise_motion ~= 0;
    
    % Create rectangles with scaled heights
    for i = 1:length(time)
        if cw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(fwd), min(fwd), clockwise_motion_fordff(i), clockwise_motion_fordff(i)], [0.5, 0, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        if ccw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(fwd), min(fwd), counterclockwise_motion_fordff(i), counterclockwise_motion_fordff(i)], [1, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end
    plot(time_supp, fwd, 'k', 'LineWidth', 1.5);  % Plot forward velocity in green
    xlabel('Time (s)');
    ylabel('Forward Velocity (mm/s)');
    title('Forward Velocity');
    grid on;
    upperbound = min(20, max(fwd));
    lowerbound = max(-10, min(fwd));
    ylim([lowerbound, upperbound]);
    
    % Subplot 3: Rotational Velocity
    subplot(3, 1, 3);  % Third subplot
    hold on
    % Calculate your scaling
    total_padding = max(rot)-min(rot);
    clockwise_motion_fordff = (clockwise_motion*total_padding)+min(rot);
    counterclockwise_motion_fordff = (counterclockwise_motion*total_padding)+min(rot);
    
    % Find non-zero periods
    cw_mask = clockwise_motion ~= 0;
    ccw_mask = counterclockwise_motion ~= 0;
    
    % Create rectangles with scaled heights
    for i = 1:length(time)
        if cw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(rot), min(rot), clockwise_motion_fordff(i), clockwise_motion_fordff(i)], [0.5, 0, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        if ccw_mask(i)
            patch([time(i)-0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)+0.5*mean(diff(time)), time(i)-0.5*mean(diff(time))], ...
                  [min(rot), min(rot), counterclockwise_motion_fordff(i), counterclockwise_motion_fordff(i)], [1, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end
    plot(time_supp, rot, 'r', 'LineWidth', 2);  % Plot rotational velocity in blue
    xlabel('Time (s)');
    ylabel('Rotational Velocity (rad/s)');
    title('Rotational Velocity');
    grid on;
    upperbound = min(300, max(rot));
    lowerbound = max(-300, min(rot));
    ylim([lowerbound, upperbound]);
    
    % Link the x-axes of all subplots
    linkaxes(findall(gcf, 'Type', 'axes'), 'x');
    
    % Add overall labels
    sgtitle('dF/F, Forward, and Rotational Velocity');  % Super title
    
    save_plot_with_title_as_filename('fwd_rot_opto', 'dff_separate', savepath);


%% now plot ccw and cw motion on same axis scale
    % Calculate overall y-limits from all data first
    all_ccw_data = cell2mat(dff_ccw_motion');
    all_cw_data = cell2mat(dff_cw_motion');
    all_data = [all_ccw_data(:); all_cw_data(:); CCW_average_dff(:); average_dff(:)];
    overall_ylim = [min(all_data), max(all_data)];
    
    % Add some padding
    padding = 0.05 * (overall_ylim(2) - overall_ylim(1));
    overall_ylim = [overall_ylim(1) - padding, overall_ylim(2) + padding];
    overall_ylim = [-3, 10];
    % Plot the "motionon" segments with pre-motion points and vertical lines
    figure;
    
    % Subplot 1: CCW motion dff (ORANGE)
    subplot(2, 1, 1);
    hold on;
    % Plot each raw "motionon" velocity segment in lighter orange
    for i = 1:length(dff_ccw_motion)
        plot(daq.t(1:length(dff_ccw_motion{i})), dff_ccw_motion{i}, 'Color', [1 0.7 0.3], 'LineWidth', 1);
        xline(daq.t(ccw_motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
        xline(daq.t((ccw_motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    plot(daq.t(1:length(dff_ccw_motion{i})), CCW_average_dff, 'Color', [1 0.5 0], 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    ylim(overall_ylim);  % Set consistent limits
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_ccw);
    legend(label_ccw);
    hold off;
    
    % Subplot 2: CW motion dff (PURPLE)
    subplot(2, 1, 2);
    hold on;
    for i = 1:length(dff_cw_motion)
        plot(daq.t(1:length(dff_cw_motion{i})), dff_cw_motion{i}, 'Color', [0.7 0.3 1], 'LineWidth', 1);
        xline(daq.t(motion_on_markers(1)), 'k--', 'LineWidth', 1.5);
        xline(daq.t((motion_on_markers(1)))+stimlen, 'k--', 'LineWidth', 1.5);
    end
    
    plot(daq.t(1:length(dff_cw_motion{i})), average_dff, 'Color', [0.5 0 1], 'LineWidth', 3, 'DisplayName', 'Average DFF');
    
    ylim(overall_ylim);  % Set consistent limits
    xlabel('Time (s)');
    ylabel('DFF');
    title(label_cw);
    legend(label_cw);
    hold off;
    
    save_plot_with_title_as_filename(label_cw, label_ccw, savepath)
end