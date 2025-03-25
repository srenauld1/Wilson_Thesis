% plotExpt_patternsliced.m
%
% Function that generates a series of summary plots containing all relevant
% data acquired during a given trial. Specifically, plots kinematic
% variables broken up by opto stim, visual stim, or both
%
% INPUTS:
%   exptData - processed data, can contain ephys, fictrac, and/or output
%   exptMeta - processed ephys meta
%
%
% OUTPUTS:

% not strictly an output but it saves 3 files
% rv_by_visual (rotational velocities split by visual stim)
% rv_by_optoandvis (rotational velocities split by visual and opto)
% rs_ffw_by_opto (forward velocity and angular speed split by opto on/off)
%
% Original: 10/31/2024 - SMR

function [] = plotExptpatternsliced(exptData,exptMeta)
    
    % reset counters
    n = 0;
    
    % set number of sublots based on number of expt variables
    checkEphys = contains(exptMeta.exptCond,'ephys','IgnoreCase',true);
    checkIInj = contains(exptMeta.exptCond,'inj','IgnoreCase',true);
    checkG4 = contains(exptMeta.exptCond,'g4','IgnoreCase',true);
    checkFicTrac = contains(exptMeta.exptCond,'fictrac','IgnoreCase',true);
    checkOpto = contains(exptMeta.exptCond,'stim','IgnoreCase',true); %not separate
    checkPython = contains(exptMeta.exptCond,'jump','IgnoreCase',true); %not separate
    checkOpenLoop = isfield(exptMeta, 'func');
    
    %% adjust panel display
    if checkG4
        if std(exptData.g4displayXPos) %and only if not stationary
            n = n+1; % update counter
            
            % pull xpos data for pre-plot processing
            g4Pos_mod = exptData.g4displayXPos;
            
            % if function used center and remove noise
            if isfield(exptMeta,'func')
                % set midpoint (in front of fly) based on object size
                midPos = (88 - (exptMeta.objSize/2 - 1))/192 *360; %center position, in degrees
                % zero data across midpoing, with right + and left -
                g4Pos_mod = g4Pos_mod - midPos;

                % some experiments use a "hidden" position for when the target
                % is behind the fly in the empty column, and therefore not visible
                if contains(exptMeta.func,'pulse')
                    hiddenPos = (184 - (exptMeta.objSize/2))/192 * 360; %hidden position, in degrees
                    g4Pos_mod(exptData.g4displayXPos>hiddenPos) = nan;
                end

                % hide noise caused by data acquisition or motion across sides
                g4Pos_mod(abs(diff(g4Pos_mod))>2) = nan;
                g4Pos_mod(isoutlier(g4Pos_mod)) = nan;

                % else no function used, center but do not remove noise
            else
                g4Pos_mod = g4Pos_mod - 180;
            end
        end
    end
    %% Behavior plotting
    %downsample time and optoopto and g4
    time = downsample(exptData.t, 30);
    opto = downsample(exptData.optoStim, 30);
    %chop up by opto blocks
    if checkOpto
        % Find the change points in the opto signal
        change_points = [1; find(diff(opto) ~= 0) + 1; length(opto)];
        % Initialize cell arrays to store segments of velocity
        velocity_segments_optoon = {};  % Velocity when opto is 1
        velocity_segments_optooff = {};  % Velocity when opto is 0
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
        disp(["OPTO ON FWD AVE:", num2str(average_fwd)])
        average_rot = mean(opto_averages_rot);
        disp(["OPTO ON ROT AVE:", num2str(average_rot)])
        average_fwd_no = mean(no_averages_fwd);
        disp(["OPTO OFF FWD AVE:", num2str(average_fwd_no)])
        average_rot_no = mean(no_averages_rot);
        disp(["OPTO OFF ROT AVE:", num2str(average_rot_no)])

        % Ensure all "optoon" segments are the same length by trimming to the shortest segment
        min_length = min(cellfun(@length, fwdvelocity_segments_optoon_extra));
        fwdvelocity_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), fwdvelocity_segments_optoon_extra, 'UniformOutput', false);
        
        % Convert to a matrix for calculating the average
        velocity_matrix_optoon = cell2mat(reshape(fwdvelocity_segments_optoon_trimmed, 1, []));
        average_velocity_optoon = mean(velocity_matrix_optoon, 2);
        
                % Ensure all "optoon" segments are the same length by trimming to the shortest segment
        min_length = min(cellfun(@length, rotspeed_segments_optoon_extra));
        rotspeed_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), rotspeed_segments_optoon_extra, 'UniformOutput', false);
        
        % Convert to a matrix for calculating the average
        speed_matrix_optoon = cell2mat(reshape(rotspeed_segments_optoon_trimmed, 1, []));
        average_speed_optoon = mean(speed_matrix_optoon, 2);

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
        
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(fwdvelocity_segments_optoon_trimmed{i})), average_velocity_optoon, 'g', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
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
        
        % Calculate and plot the average speed for "optoon" in bold blue
        plot(time(1:length(rotspeed_segments_optoon_trimmed{i})), average_speed_optoon, 'b', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Speed');
        
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Speed');
        title('Opto-On Angular Speed');
        legend('Average Opto-On Angular Speed');
        hold off;
        %saveas(gcf, ['rs_fw_by_opto.png']);

        %% now separate walking vs not
        still_pre_opto = {}; 
        moving_pre_opto ={};
        ave_still = {}; 
        ave_moving ={};
        for i = 1:length(fwdvelocity_segments_optoon_trimmed)
            pre_speed = mean(fwdvelocity_segments_optoon_trimmed{i}((pre_opto_points-100):pre_opto_points));
            if pre_speed < 1
                disp("Still")
                still_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_still{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                
            else
                moving_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_still{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                disp("moving")
            end
        end

        % take average
                % Convert to a matrix for calculating the average
        moving = cell2mat(reshape(ave_moving, 1, []));
        moving_ave = mean(moving, 2);
                % Convert to a matrix for calculating the average
        still = cell2mat(reshape(ave_still, 1, []));
        still_ave = mean(still, 2);

        %% plot this
         % Plot the "optoon" segments with pre-opto points and vertical lines
        figure;
        
        % Subplot 1: moving
        subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
        hold on;
        % Plot each raw "optoon" velocity segment in lighter green
        for i = 1:length(moving_pre_opto)
            plot(time(1:length(moving_pre_opto{i})), moving_pre_opto{i}, 'Color', [0.75, 0, 0.75], 'LineWidth', 1);  % Lighter green
            xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
            xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
        end

        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(moving_ave)), moving_ave, 'Color', [1 0 1], 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
        
        % Customize the first subplot
        xlabel('Time (s)');
        ylabel('Forward Velocit (mm/s)');
        title('Moving pre-opto');
        hold off;
        
        % Subplot 2: Speed
        subplot(2, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
        hold on;
        % Plot each raw "optoon" speed segment in lighter blue
        for i = 1:length(still_pre_opto)
            plot(time(1:length(still_pre_opto{i})), still_pre_opto{i}, 'Color', [1 0.5 0.5], 'LineWidth', 1);  % Lighter blue
            xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
            xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
        end
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(still_ave)), still_ave, 'r', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
       
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Forward Velocity (mm/s)');
        title('Still pre-opto');
        hold off;
        %saveas(gcf, ['rs_fw_by_opto.png']);


    end
%% chop up by visual stim blocks
    if checkOpenLoop
        g4Pos_mod_ds = downsample(g4Pos_mod, 30);
        
        % Define stationary ranges (around 90 and 270)
        stationary_range_1 = [-94, -93];
        stationary_range_2 = [90, 92];
        
        % Identify stationary states
        is_stationary = (g4Pos_mod_ds >= stationary_range_1(1) & g4Pos_mod_ds <= stationary_range_1(2)) | ...
                        (g4Pos_mod_ds >= stationary_range_2(1) & g4Pos_mod_ds <= stationary_range_2(2));
        
        % Identify moving clockwise (90 to 270) and moving counterclockwise (270 to 90)
        is_moving = (g4Pos_mod_ds > stationary_range_1(2) & g4Pos_mod_ds < stationary_range_2(1));
        %separate into 2
        
        % Initialize two arrays of zeros with the same length as the input is_moving
        is_moving_clockwise = zeros(size(is_moving));
        is_moving_counterclockwise = zeros(size(is_moving));
        
        % Initialize a flag to alternate between sequences
        use_is_moving_clockwise = true;
        
        % Loop through each continuous segment of 1s
        start_idx = 1;
        while start_idx <= length(is_moving)
            % Find the next segment of 1s
            while start_idx <= length(is_moving) && is_moving(start_idx) == 0
                start_idx = start_idx + 1;
            end
            
            % If we are at the end of the array, break out of the loop
            if start_idx > length(is_moving)
                break;
            end
            
            % Find the end of the current segment of 1s
            end_idx = start_idx;
            while end_idx <= length(is_moving) && is_moving(end_idx) == 1
                end_idx = end_idx + 1;
            end
            end_idx = end_idx - 1; % Adjust for the last increment
            
            % Assign the segment to the appropriate sequence
            if use_is_moving_clockwise
                is_moving_clockwise(start_idx:end_idx) = 1;
            else
                is_moving_counterclockwise(start_idx:end_idx) = 1;
            end
            
            % Toggle the flag to alternate sequences
            use_is_moving_clockwise = ~use_is_moving_clockwise;
            
            % Move the start index to the end of the current segment
            start_idx = end_idx + 1;
        end
    
        % Initialize cell arrays to store segmented velocity data
        velocity_segments_stationary = {};
        velocity_segments_moving = {};
        velocity_segments_clockwise = {};
        velocity_segments_counterclockwise = {};
        
        % Detect change points where state switches
        %change_points = [1; find(diff(is_stationary | is_moving_clockwise | is_moving_counterclockwise)) + 1; length(g4Pos_mod_ds) + 1];
        change_points = [1; find(diff(is_stationary)) + 1; length(g4Pos_mod_ds) + 1];
        % Loop through each segment and separate based on the state
        for i = 1:length(change_points) - 1
            % Get indices for the current segment
            segment_indices = change_points(i):change_points(i+1) - 1;
            
            % Check if the current segment is stationary, moving clockwise, or moving counterclockwise
            if is_stationary(segment_indices(1))
                % Store the velocity segment for stationary condition
                velocity_segments_stationary{end+1} = exptData.angularVelocity(segment_indices);
            elseif is_moving(segment_indices(1))
                % Store the velocity segment for clockwise movement
                velocity_segments_moving{end+1} = exptData.angularVelocity(segment_indices);
            if is_moving_clockwise(segment_indices(1))
                % Store the velocity segment for clockwise movement
                velocity_segments_clockwise{end+1} = exptData.angularVelocity(segment_indices);
            end
            if is_moving_counterclockwise(segment_indices(1))
                % Store the velocity segment for counterclockwise movement
                velocity_segments_counterclockwise{end+1} = exptData.angularVelocity(segment_indices);
            end
            end
        end
            % Define the target length
        target_length = 652;
        % crop 
        % Remove the first and last cell from velocity_segments_stationary
        velocity_segments_stationary = velocity_segments_stationary(2:end-1);
        % Filter out elements from velocity_segments_stationary that are 1x1
        velocity_segments_stationary = velocity_segments_stationary( ...
            cellfun(@(x) isnumeric(x) && ~(isequal(size(x), [1, 1])), velocity_segments_stationary));
        velocity_segments_moving = velocity_segments_moving( ...
            cellfun(@(x) isnumeric(x) && ~(isequal(size(x), [1, 1])), velocity_segments_moving));
        velocity_segments_clockwise = velocity_segments_clockwise( ...
            cellfun(@(x) isnumeric(x) && ~(isequal(size(x), [1, 1])), velocity_segments_clockwise));
        velocity_segments_counterclockwise = velocity_segments_counterclockwise( ...
            cellfun(@(x) isnumeric(x) && ~(isequal(size(x), [1, 1])), velocity_segments_counterclockwise));

        % Trim each segment in the stationary, clockwise, and counterclockwise cell arrays
        velocity_segments_stationary_trimmed = cellfun(@(x) x(1:target_length), velocity_segments_stationary, 'UniformOutput', false);
        velocity_segments_moving_trimmed = cellfun(@(x) x(1:target_length), velocity_segments_moving, 'UniformOutput', false);
        velocity_segments_counterclockwise_trimmed = cellfun(@(x) x(1:target_length), velocity_segments_counterclockwise, 'UniformOutput', false);
        velocity_segments_clockwise_trimmed = cellfun(@(x) x(1:target_length), velocity_segments_clockwise, 'UniformOutput', false);
        
        % Convert each trimmed cell array to a matrix for mean calculation
        stationary_matrix = cell2mat(reshape(velocity_segments_stationary_trimmed, 1, []));
        moving_matrix = cell2mat(reshape(velocity_segments_moving_trimmed, 1, []));
        clockwise_matrix = cell2mat(reshape(velocity_segments_clockwise_trimmed, 1, []));
        counterclockwise_matrix = cell2mat(reshape(velocity_segments_counterclockwise_trimmed, 1, []));
        
        % Calculate the mean velocity across segments at each time point
        mean_velocity_stationary = mean(stationary_matrix, 2);
        mean_velocity_moving = mean(moving_matrix, 2);
        mean_velocity_cw = mean(clockwise_matrix, 2);
        mean_velocity_ccw = mean(counterclockwise_matrix, 2);
    
        % Plot the segmented velocity traces
        figure;
        hold on;
        
        % Plot each stationary segment in a lighter blue
        for i = 1:length(velocity_segments_stationary_trimmed)
            plot(velocity_segments_stationary_trimmed{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Light blue
        end
        
        % Plot each clockwise segment in a lighter green
        for i = 1:length(velocity_segments_clockwise)
            plot(velocity_segments_clockwise{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Light green
        end
        
        % Plot each counterclockwise segment in a lighter magenta
        for i = 1:length(velocity_segments_counterclockwise)
            plot(velocity_segments_counterclockwise{i}, 'Color', [1 0.5 1], 'LineWidth', 1);  % Light magenta
        end
        
        % Plot the mean velocity traces with solid colors and thicker lines
        h1 = plot(mean_velocity_stationary, 'b', 'LineWidth', 3, 'DisplayName', 'Stationary');
        h2 = plot(mean_velocity_cw, 'g', 'LineWidth', 3, 'DisplayName', 'Clockwise');
        h3 = plot(mean_velocity_ccw, 'm', 'LineWidth', 3, 'DisplayName', 'Counterclockwise');
        
        % Customize the plot
        xlabel('Time (samples)');
        ylabel('Angular Velocity (deg/s)');
        title('Angular Velocity Segmented by Optomotor Stim Direction');
        legend([h1, h2, h3]);  % Only include the mean traces in the legend
        hold off;
        %saveas(gcf, ['rv_by_visual.png']);

            %% separate flat path based on moving or stationary       
        % % Plot the path with different colors for each opto state
        % figure;
        % hold on;
        % 
        % % Loop through the points to plot segments based on opto values
        % for i = 1:length(x)-1
        %     if is_stationary(i) == 1
        %         % Plot black segment when pattern is stationary
        %         plot(x(i:i+1), y(i:i+1), 'k-', 'LineWidth', 1.5);
        %     else
        %         % Plot red segment when pattern is moving
        %         plot(x(i:i+1), y(i:i+1), 'b-', 'LineWidth', 1.5);
        %     end
        % end
        % 
        % % Customize the plot
        % xlabel('Sideways Position');
        % ylabel('forward Position');
        % title('Object Path Colored by Pattern Movement State');
        % legend({'Pattern Stationary', 'Pattern moving (blue)'});
        % grid on;
        % hold off;

    end


    %% separate by both opto and condition
    if checkOpenLoop
        if checkOpto
            % Identify stationary states
            is_stationary = (g4Pos_mod_ds >= stationary_range_1(1) & g4Pos_mod_ds <= stationary_range_1(2)) | ...
                            (g4Pos_mod_ds >= stationary_range_2(1) & g4Pos_mod_ds <= stationary_range_2(2));
            
            % Identify moving clockwise (90 to 270) and moving counterclockwise (270 to 90)
            is_moving = (g4Pos_mod_ds > stationary_range_1(2) & g4Pos_mod_ds < stationary_range_2(1));
            %separate into 2
            
            % Initialize two arrays of zeros with the same length as the input is_moving
            is_moving_clockwise = zeros(size(is_moving));
            is_moving_counterclockwise = zeros(size(is_moving));
            
            % Initialize a flag to alternate between sequences
            use_is_moving_clockwise = true;
            
            % Loop through each continuous segment of 1s
            start_idx = 1;
            while start_idx <= length(is_moving)
                % Find the next segment of 1s
                while start_idx <= length(is_moving) && is_moving(start_idx) == 0
                    start_idx = start_idx + 1;
                end
                
                % If we are at the end of the array, break out of the loop
                if start_idx > length(is_moving)
                    break;
                end
                
                % Find the end of the current segment of 1s
                end_idx = start_idx;
                while end_idx <= length(is_moving) && is_moving(end_idx) == 1
                    end_idx = end_idx + 1;
                end
                end_idx = end_idx - 1; % Adjust for the last increment
                
                % Assign the segment to the appropriate sequence
                if use_is_moving_clockwise
                    is_moving_clockwise(start_idx:end_idx) = 1;
                else
                    is_moving_counterclockwise(start_idx:end_idx) = 1;
                end
                
                % Toggle the flag to alternate sequences
                use_is_moving_clockwise = ~use_is_moving_clockwise;
                
                % Move the start index to the end of the current segment
                start_idx = end_idx + 1;
            end
    
            % separate by opto on or off (normalize to 0's and 1's for
            % opto)
            norm = 1/max(opto);
            opto_norm = opto*norm;
            cw_opto = is_moving_clockwise.*opto_norm;
            ccw_opto = is_moving_counterclockwise.*opto_norm;
            cw_noopto = is_moving_clockwise.* ~opto_norm;
            ccw_noopto = is_moving_counterclockwise.* ~opto_norm;

            % Remove short segments from opto off

            % Minimum length of consecutive 1's to keep
            min_length = 200;
            
            % For ccw_noopto
            data = ccw_noopto(:);  % Ensure data is a column vector
            % Find sequences of 1's and their lengths
            diff_data = [0; diff(data); 0];  % Add padding to detect start and end in column format
            start_indices = find(diff_data == 1);  % Start of sequences of 1's
            end_indices = find(diff_data == -1) - 1;  % End of sequences of 1's
            
            % Loop over each sequence and set sequences of length < min_length to 0
            for i = 1:length(start_indices)
                % Length of the current sequence
                sequence_length = end_indices(i) - start_indices(i) + 1;
                
                % Replace with 0's if the sequence length is less than min_length
                if sequence_length < min_length
                    data(start_indices(i):end_indices(i)) = 0;
                end
            end
            ccw_noopto = data;  % Ensure output is a column vector
            
            % For cw_noopto
            data = cw_noopto(:);  % Ensure data is a column vector
            % Find sequences of 1's and their lengths
            diff_data = [0; diff(data); 0];  % Add padding to detect start and end in column format
            start_indices = find(diff_data == 1);  % Start of sequences of 1's
            end_indices = find(diff_data == -1) - 1;  % End of sequences of 1's
            
            % Loop over each sequence and set sequences of length < min_length to 0
            for i = 1:length(start_indices)
                % Length of the current sequence
                sequence_length = end_indices(i) - start_indices(i) + 1;
                
                % Replace with 0's if the sequence length is less than min_length
                if sequence_length < min_length
                    data(start_indices(i):end_indices(i)) = 0;
                end
            end
            cw_noopto = data;  % Ensure output is a column vector


        
            % Initialize cell arrays to store segmented velocity data
            velocity_segments_stationary = {};
            velocity_segments_clockwise_opto = {};
            velocity_segments_counterclockwise_opto = {};
            velocity_segments_clockwise_noopto = {};
            velocity_segments_counterclockwise_noopto = {};
            
            % Detect change points where state switches
            %change_points = [1; find(diff(is_stationary | cw_noopto | cw_opto | ccw_opto | ccw_noopto)) + 1; length(g4Pos_mod_ds) + 1];
            change_points = [1; find(diff(is_stationary)) + 1; length(g4Pos_mod_ds) + 1];
            % Loop through each segment and separate based on the state
            for i = 1:length(change_points) - 1
                % Get indices for the current segment
                segment_indices = change_points(i):change_points(i+1) - 1;
                
                % Check if the current segment is stationary, moving clockwise, or moving counterclockwise
                if is_stationary(segment_indices(1))
                    % Store the velocity segment for stationary condition
                    velocity_segments_stationary{end+1} = exptData.angularVelocity(segment_indices);
                end
                if cw_opto(segment_indices(1))
                    % Store the velocity segment for clockwise movement
                    velocity_segments_clockwise_opto{end+1} = exptData.angularVelocity(segment_indices);
                elseif ccw_opto(segment_indices(1))
                    % Store the velocity segment for counterclockwise movement
                    velocity_segments_counterclockwise_opto{end+1} = exptData.angularVelocity(segment_indices);
                elseif ccw_noopto(segment_indices(1))
                    % Store the velocity segment for counterclockwise movement
                    velocity_segments_clockwise_noopto{end+1} = exptData.angularVelocity(segment_indices);
                elseif cw_noopto(segment_indices(1))
                    % Store the velocity segment for counterclockwise movement
                    velocity_segments_counterclockwise_noopto{end+1} = exptData.angularVelocity(segment_indices);
                
                
                end
            end
                % Define the target length
            target_length = 652;
            % crop 
            % Remove the first and last cell from velocity_segments_stationary
            velocity_segments_stationary = velocity_segments_stationary(2:end-1);
            
            % Trim each segment in the stationary, clockwise, and counterclockwise cell arrays
            velocity_segments_stationary_trimmed = cellfun(@(x) x(1:target_length), velocity_segments_stationary, 'UniformOutput', false);
            velocity_segments_clockwise_opto_trim = cellfun(@(x) x(1:target_length), velocity_segments_clockwise_opto, 'UniformOutput', false);
            velocity_segments_counterclockwise_opto_trim = cellfun(@(x) x(1:target_length), velocity_segments_counterclockwise_opto, 'UniformOutput', false);
            velocity_segments_clockwise_noopto_trim = cellfun(@(x) x(1:target_length), velocity_segments_clockwise_noopto, 'UniformOutput', false);
            velocity_segments_counterclockwise_noopto_trim = cellfun(@(x) x(1:target_length), velocity_segments_counterclockwise_noopto, 'UniformOutput', false);

            % Convert each trimmed cell array to a matrix for mean calculation
            stationary_matrix = cell2mat(reshape(velocity_segments_stationary_trimmed, 1, []));
            cw_opto_matrix = cell2mat(reshape(velocity_segments_clockwise_opto_trim, 1, []));
            ccw_opto_matrix = cell2mat(reshape(velocity_segments_counterclockwise_opto_trim, 1, []));
            cw_noopto_matrix = cell2mat(reshape(velocity_segments_clockwise_noopto_trim, 1, []));
            ccw_noopto_matrix = cell2mat(reshape(velocity_segments_counterclockwise_noopto_trim, 1, []));
            
            % Calculate the mean velocity across segments at each time point
            mean_velocity_stationary = mean(stationary_matrix, 2);
            mean_velocity_cw_opto = mean(cw_opto_matrix, 2);
            mean_velocity_ccw_opto= mean(ccw_opto_matrix, 2);
            mean_velocity_cw_noopto = mean(cw_noopto_matrix, 2);
            mean_velocity_ccw_noopto = mean(ccw_noopto_matrix, 2);

            % Initialize the figure
            figure;
            hold on;
            
            % Plot each raw velocity trace with lighter colors
            
            % Stationary segments (lighter blue)
            for i = 1:length(velocity_segments_stationary)
                plot(velocity_segments_stationary{i}, 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Light blue
            end
            
            % CW with opto segments (lighter green)
            for i = 1:length(velocity_segments_clockwise_opto)
                plot(velocity_segments_clockwise_opto{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Light green
            end
            
            % CCW with opto segments (lighter magenta)
            for i = 1:length(velocity_segments_counterclockwise_opto)
                plot(velocity_segments_counterclockwise_opto{i}, 'Color', [1 0.5 1], 'LineWidth', 1);  % Light magenta
            end
            
            % CW without opto segments (lighter red)
            for i = 1:length(velocity_segments_clockwise_noopto)
                plot(velocity_segments_clockwise_noopto{i}, 'Color', [1 0.5 0.5], 'LineWidth', 1);  % Light red
            end
            
            % CCW without opto segments (lighter cyan)
            for i = 1:length(velocity_segments_counterclockwise_noopto)
                plot(velocity_segments_counterclockwise_noopto{i}, 'Color', [0.5 1 1], 'LineWidth', 1);  % Light cyan
            end
            
            % Overlay the average velocity traces with solid colors and thicker lines
            % NEED TO DOUBLE CHECK THE PLOTTING AND ASSIGNMENT OF CW NO
            % OPTO AND CCW NO OPTO
            
            % Plot each mean velocity trace
            h1 = plot(mean_velocity_stationary, 'b', 'LineWidth', 5, 'DisplayName', 'Stationary');
            h2 = plot(mean_velocity_cw_opto, 'g', 'LineWidth', 5, 'DisplayName', 'CW Opto');
            h3 = plot(mean_velocity_ccw_opto, 'm', 'LineWidth', 5, 'DisplayName', 'CCW Opto');
            h4 = plot(mean_velocity_cw_noopto, 'r', 'LineWidth', 5, 'DisplayName', 'CCW No Opto');
            h5 = plot(mean_velocity_ccw_noopto, 'c', 'LineWidth', 5, 'DisplayName', 'CW No Opto');
            
            % Customize the plot
            xlabel('Time (samples)');
            ylabel('Angular Velocity (deg/s)');
            title('Angular Velocity by Opto and Visual Condition');
            legend([h1, h2, h3, h4, h5]);  % Only include the mean traces in the legend
            grid on;
            hold off;

            %saveas(gcf, ['rv_by_optoandvis.png']);

        end
    end