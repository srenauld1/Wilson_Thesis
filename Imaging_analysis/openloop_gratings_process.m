function [daq, ts] = openloop_gratings_process(daq, ts, savepath)
    % Find the time points when visual pattern is moving
    movingGratingClockwise = zeros(size(daq.vyv));
    movingGratingCCW = zeros(size(daq.vyv));
    % Check conditions: greater than 0.1 and less than 5
    movingGratingClockwise(daq.vyv > 0.1 & daq.vyv < 3) = 1;
    movingGratingCCW(daq.vyv < -0.1 & daq.vyv >- 3) = 1;
    
    % Initialize cell arrays to store segments of velocity
    velocity_segments_optoon = {};  % Velocity when opto is 1
    velocity_segments_optooff = {};  % Velocity when opto is 0
    % now take care off immediate on off for 10 sec stim windows
    if sum(opto)==36000
        % Assuming your array is called change_points
        keep_idx = false(size(change_points,2));
        for i = 1:4:length(change_points)
            if i+1 <= length(change_points)
                keep_idx(i) = true;    % Keep first number
                keep_idx(i+1) = true;  % Keep second number
            end
        end
        change_points_filtered = change_points(keep_idx);
        change_points= change_points_filtered;
        on_off = 1;
        disp("you have stop start opto here 4.5 then 1 off then 4.5")

    end
    % Loop through each segment and separate based on opto value
    for i = 1:length(change_points) - 1
        % Get the indices for the current segment
        segment_indices = change_points(i):change_points(i+1) - 1;
        % Check the opto value in this segment
        if opto(segment_indices(1)) >0
            % Store velocity data where opto is 1
            velocity_segments_optoon{end+1} = exptData.forwardVelocity(segment_indices);
        else
            % Store velocity data where opto is 0
            velocity_segments_optooff{end+1} = exptData.forwardVelocity(segment_indices);
        end
    end

    %% now chop up with preface & plot
    % Define the number of timepoints to include before opto-on
    pre_opto_points = 1000;
    post_opto_points = 1000;
    
    % Initialize cell arrays to store segments of velocity
    fwdvelocity_segments_optoon_extra = {};  % Velocity when opto is 1 (including pre-opto points
    rotspeed_segments_optoon_extra = {}; 
    opto_on_markers = [];  % To store the position of opto-on points for vertical lines
    opto_off_markers = [];
    opto_averages_fwd = [];
    opto_averages_rot = [];
    no_averages_fwd = [];
    no_averages_rot = [];
    % Loop through each segment and separate based on opto value
    for i = 1:length(change_points) - 1
        % Get the indices for the current segment
        segment_indices = change_points(i):change_points(i+1) - 1;
    
        % Check the opto value in this segment
        if opto(segment_indices(1)) > 0

            angular_speed = min(exptData.angularSpeed, 600);
            % Calculate the start index for pre-opto points (ensuring it doesn’t go below 1)
            pre_opto_start = max(1, segment_indices(1) - pre_opto_points);
            % Calculate the end index for post-opto points (ensuring it doesn’t go below 1)
            post_opto_end = min(length(opto), segment_indices(end) + post_opto_points);


            % only opto on fwd and rot
            opto_on_fwd = exptData.forwardVelocity(segment_indices);
            opto_on_rot = angular_speed(segment_indices);

            % opto off fwd and rot
            no_opto = [(pre_opto_start:segment_indices(1)), (segment_indices(end):post_opto_points)];
            off_fwd = exptData.forwardVelocity(no_opto);
            off_rot = angular_speed(no_opto);

            % take average for trial
            average_fwd = mean(opto_on_fwd);
            average_rot = mean(opto_on_rot);
            off_fwd_ave = mean(off_fwd);
            off_rot_ave = mean(off_rot);
        
            
            % Include the 200 points before the "optoon" segment
            extended_segment_indices = pre_opto_start:post_opto_end;
            
            % Store the extended velocity data for opto-on
            fwdvelocity_segments_optoon_extra{end+1} = exptData.forwardVelocity(extended_segment_indices);
            rotspeed_segments_optoon_extra{end+1} = angular_speed(extended_segment_indices);

            % for calculating averages
            opto_averages_fwd(end+1) = average_fwd;
            opto_averages_rot(end+1) = average_rot;
            no_averages_fwd(end+1) = off_fwd_ave;
            no_averages_rot(end+1) = off_rot_ave;
            
            % Record the exact point where opto switches on for plotting
            opto_on_markers(end+1) = segment_indices(1) - pre_opto_start + 1;
            opto_off_markers(end+1) = segment_indices(end)-segment_indices(1)+pre_opto_points;% Relative position within segment
        end
    end

    average_fwd = mean(opto_averages_fwd);
    %disp(["OPTO ON FWD AVE:", num2str(average_fwd)])
    average_rot = mean(opto_averages_rot);
    %disp(["OPTO ON ROT AVE:", num2str(average_rot)])
    average_fwd_no = mean(no_averages_fwd);
    %disp(["OPTO OFF FWD AVE:", num2str(average_fwd_no)])
    average_rot_no = mean(no_averages_rot);
    %disp(["OPTO OFF ROT AVE:", num2str(average_rot_no)])

    % Ensure all "optoon" segments are the same length by trimming to the shortest segment
    min_length = min(cellfun(@length, fwdvelocity_segments_optoon_extra));
    if length(change_points) >100
        fwdvelocity_segments_optoon_extra_endremoved = fwdvelocity_segments_optoon_extra(1:end-1);
        rotspeed_segments_optoon_extra_endremoved = rotspeed_segments_optoon_extra(1:end-1);
        min_length = min(cellfun(@length, fwdvelocity_segments_optoon_extra_endremoved));
        rotspeed_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), rotspeed_segments_optoon_extra_endremoved, 'UniformOutput', false);
        fwdvelocity_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), fwdvelocity_segments_optoon_extra_endremoved, 'UniformOutput', false);
    else
        fwdvelocity_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), fwdvelocity_segments_optoon_extra, 'UniformOutput', false);
        % the shortest segmennt
        rotspeed_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), rotspeed_segments_optoon_extra, 'UniformOutput', false);
        
    end
    
    % Convert to a matrix for calculating the average
    velocity_matrix_optoon = cell2mat(reshape(fwdvelocity_segments_optoon_trimmed, length(fwdvelocity_segments_optoon_trimmed), []));
    average_velocity_optoon = mean(velocity_matrix_optoon, 1);
    
    % Ensure all "optoon" segments are the same length by trimming to
    
    % Convert to a matrix for calculating the average
    speed_matrix_optoon = cell2mat(reshape(rotspeed_segments_optoon_trimmed, length(rotspeed_segments_optoon_trimmed), []));
    average_speed_optoon = mean(speed_matrix_optoon, 1);

    % Plot the "optoon" segments with pre-opto points and vertical lines
    figure;
    
    % Subplot 1: Velocity
    subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
    hold on;
    % Plot each raw "optoon" velocity segment in lighter green
    for i = 1:length(fwdvelocity_segments_optoon_extra)
        plot(time(1:length(fwdvelocity_segments_optoon_extra{i})), fwdvelocity_segments_optoon_extra{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
    end
    if on_off
        xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
        xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
    end
    
    % Plot the average velocity for "optoon" in bold green
    plot(time(1:length(fwdvelocity_segments_optoon_extra{1})), average_velocity_optoon, 'g', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
    
    % Customize the first subplot
    xlabel('Time (s)');
    ylabel('Forward Velocit (mm/s)');
    title('Opto-On Forward Velocity');
    legend('Average Opto-On Forward Velocity');
    hold off;
    
    % Subplot 2: Speed
    subplot(2, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
    hold on;
    % Plot each raw "optoon" speed segment in lighter blue
    for i = 1:length(rotspeed_segments_optoon_extra)
        plot(time(1:length(rotspeed_segments_optoon_extra{i})), rotspeed_segments_optoon_extra{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
    end

    if on_off
        xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
        xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
    end
    
    % Calculate and plot the average speed for "optoon" in bold blue
    plot(time(1:length(rotspeed_segments_optoon_trimmed{1})), average_speed_optoon, 'b', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Speed');
    
    % Customize the second subplot
    xlabel('Time (s)');
    ylabel('Speed');
    title('Opto-On Angular Speed');
    legend('Average Opto-On Angular Speed');
    hold off;
    save_plot_with_title_as_filename('forward', 'rotation', savepath)

end