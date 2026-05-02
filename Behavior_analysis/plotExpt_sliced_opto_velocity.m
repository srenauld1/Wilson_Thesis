% plotExpt_sliced_opto_velocity.m
%
% Function that generates a series of summary plots containing all relevant
% data acquired during a given trial. Specifically, plots kinematic
% variables broken up by opto stim, visual stim, or both
%
% INPUTS:
%   exptData - processed data, can contain ephys, fictrac, and/or output
%   exptMeta - processed ephys meta
%   savepath - where to save the files to on the Wilson lab server
%
%
% OUTPUTS:
%   this saves 3 files to savepath
%   returns exptData and exptMeta 
% Original: 10/31/2024 - SMR (plotExptpatternsliced)
%           3/19/2025 - SMR adjusted to only plot closed loop data]=
%           3/25/2026 - SMR adjusted to account for minor errors using LLM

function [exptData, exptMeta, fwdvelocity_segments_optoon_extra, rotvel_segments_optoon_extra] = plotExpt_sliced_opto_velocity(exptData,exptMeta, savepath, minVel)  
    % set number of sublots based on number of expt variables
    checkOpto = contains(exptMeta.exptCond,'stim','IgnoreCase',true); %not separate

    %% define variables
    % define angular velocity
    angular_vel = exptData.angularVelocity_raw;
    fwd_vel = exptData.forwardVelocity_raw;
    n_angular = size(angular_vel, 2);

    % original time vector
    % Original opto and its time base
    opto = exptData.optoStim;    % [1 x N] or [N x 1]
    t_opto = exptData.t;     % [1 x N] or [N x 1]
    
    % Target time points (upsampled), matching angularVelocity_raw
    n_angular = size(exptData.angularVelocity_raw, 2);
    time_up = linspace(t_opto(1), t_opto(end), n_angular);
    
    % Upsample opto using nearest neighbor (step-wise, no interpolation)
    opto_up = interp1(t_opto(:), opto(:), time_up(:), 'nearest', 0);
    % opto_up is now [n_angular x 1] or [1 x n_angular]

    time = time_up;
    opto = opto_up';

    % initialize container for sliced data
    exptData.opto_sliced = struct([]);

    %% chop up by opto blocks
    if checkOpto
        % Find the change points in the opto signal
        change_points = [1; find(diff(opto(:)) ~= 0) + 1; length(opto)];
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
                velocity_segments_optoon{end+1} = fwd_vel(segment_indices);
            else
                % Store velocity data where opto is 0
                velocity_segments_optooff{end+1} = fwd_vel(segment_indices);
            end
        end

        %% now chop up with preface & plot
        % Define the number of timepoints to include before opto-on
        % look at the sampling rate
        rate = length(angular_vel)/ max(time);
        pre_opto_points = rate;
        post_opto_points = rate;


        % Initialize cell arrays to store segments of velocity
        fwdvelocity_segments_optoon_extra = {};  % Velocity when opto is 1 (including pre-opto points
        rotvel_segments_optoon_extra = {}; 
        opto_on_markers = [];  % To store the position of opto-on points for vertical lines
        opto_off_markers = [];
        opto_averages_fwd = [];
        opto_averages_rot = [];
        no_averages_fwd = [];
        no_averages_rot = [];
        opto_averages_fwdaccel = [];
        opto_averages_rotaccel = [];
        no_averages_fwdaccel = [];
        no_averages_rotaccel = [];
        % Loop through each segment and separate based on opto value
        for i = 1:length(change_points) - 1
            % Get the indices for the current segment
            segment_indices = change_points(i):change_points(i+1) - 1;
        
            % Check the opto value in this segment
            if opto(segment_indices(1)) > 0


                % Calculate the start index for pre-opto points (ensuring it doesn’t go below 1)
                pre_opto_start = max(1, segment_indices(1) - pre_opto_points);
                % Calculate the end index for post-opto points (ensuring it doesn’t go below 1)
                post_opto_end = min(length(opto), segment_indices(end) + post_opto_points);


                % only opto on fwd and rot
                opto_on_fwd = fwd_vel(segment_indices);
                opto_on_rot = abs(angular_vel(segment_indices));

                % look at the sampling rate
                rate = round(length(angular_vel) / max(time));  % <-- round here
                pre_opto_points = rate;
                post_opto_points = rate;

                % opto off fwd and rot
                pre_opto_start = max(1, segment_indices(1) - pre_opto_points);
                pre_opto = pre_opto_start : segment_indices(1) - 1;
                
                if isempty(pre_opto)
                    continue  % skip this segment, no pre-opto baseline available
                end
                
                pre_opto_fwd = fwd_vel(pre_opto);
                pre_opto_rot = abs(angular_vel(pre_opto));


                % take average for trial
                average_fwd = mean(opto_on_fwd);
                average_opto_fwd_accel = mean(diff(opto_on_fwd));
                average_rot = mean(opto_on_rot);
                average_opto_rot_accel = mean(diff(opto_on_rot));
                off_fwd_ave = mean(pre_opto_fwd);
                average_preopto_fwd_accel = mean(diff(pre_opto_fwd));
                off_rot_ave = mean(pre_opto_rot);
                average_preopto_rot_accel = mean(diff(pre_opto_rot));
            
                
                % Include the 1000 points before the "optoon" segment
                extended_segment_indices = pre_opto_start:post_opto_end;
                extended_segment_indices = round(extended_segment_indices);
                
                % Store the extended velocity data for opto-on
                fwdvelocity_segments_optoon_extra{end+1} = fwd_vel(extended_segment_indices);
                rotvel_segments_optoon_extra{end+1} = angular_vel(extended_segment_indices);

                % ---------------------------------------------------------------------
                % Slice ALL exptData time-series fields for this opto segment
                % ---------------------------------------------------------------------
                seg_idx = numel(exptData.opto_sliced) + 1;
                
                % store indices + time + opto explicitly
                exptData.opto_sliced(seg_idx).indices = extended_segment_indices(:)';
                exptData.opto_sliced(seg_idx).time    = time(extended_segment_indices);
                exptData.opto_sliced(seg_idx).opto    = opto(extended_segment_indices);
                
                fn = fieldnames(exptData);
                for f = 1:numel(fn)
                    thisField = fn{f};
                    % skip the container we are building
                    if strcmp(thisField, 'opto_sliced')
                        continue;
                    end
                
                    val = exptData.(thisField);
                
                    % Only slice 1D vectors that match the length of opto/time
                    if isnumeric(val) || islogical(val)
                        s = size(val);
                        if isvector(val) && numel(val) == numel(time)
                            exptData.opto_sliced(seg_idx).(thisField) = val(extended_segment_indices);
                        else
                            % leave non-matching numeric as-is (e.g. scalars, 2D arrays)
                            % You can comment this out if you don't want them at all:
                            % exptData.opto_sliced(seg_idx).(thisField) = val;
                        end
                    else
                        % Non-numeric fields are copied unchanged (meta info, etc.)
                        exptData.opto_sliced(seg_idx).(thisField) = val;
                    end
                end
                % ---------------------------------------------------------------------

                % for calculating averages
                opto_averages_fwd(end+1) = average_fwd;
                opto_averages_rot(end+1) = average_rot;
                no_averages_fwd(end+1) = off_fwd_ave;
                no_averages_rot(end+1) = off_rot_ave;

                % average accelerations
                opto_averages_fwdaccel(end+1) = average_opto_fwd_accel;
                opto_averages_rotaccel(end+1) = average_opto_rot_accel;
                no_averages_fwdaccel(end+1) = average_preopto_fwd_accel;
                no_averages_rotaccel(end+1) = average_preopto_rot_accel;
                
                % Record the exact point where opto switches on for plotting
                opto_on_markers(end+1) = segment_indices(1) - pre_opto_start + 1;
                opto_off_markers(end+1) = segment_indices(end)-segment_indices(1)+pre_opto_points +1;% Relative position within segment
            end
        end

        average_fwd = mean(opto_averages_fwd);
        disp(["OPTO ON FWD AVE:", num2str(average_fwd)])
        average_rot = mean(opto_averages_rot);
        disp(["OPTO ON ROT AVE:", num2str(average_rot)])
        average_fwd_no = mean(no_averages_fwd);
        disp(["PRE-OPTO FWD AVE:", num2str(average_fwd_no)])
        average_rot_no = mean(no_averages_rot);
        disp(["PRE-OPTO ROT AVE:", num2str(average_rot_no)])
        % acceleration
        average_fwd_accel = mean(opto_averages_fwdaccel);
        %disp(["OPTO ON FWD ACCELERATION AVE:", num2str(average_fwd_accel)])
        average_rot_accel = mean(opto_averages_rotaccel);
        %disp(["OPTO ON ROT ACCELERATION AVE:", num2str(average_rot_accel)])
        average_fwd_no_accel = mean(no_averages_fwdaccel);
        %disp(["PRE-OPTO FWD ACCELERATION AVE:", num2str(average_fwd_no_accel)])
        average_rot_no_accel = mean(no_averages_rotaccel);
        %disp(["PRE-OPTO ROT ACCELERATION AVE:", num2str(average_rot_no_accel)])

        % Ensure all "optoon" segments are the same length by trimming to the shortest segment
        min_length = min(cellfun(@length, fwdvelocity_segments_optoon_extra));
        if length(change_points) >100
            fwdvelocity_segments_optoon_extra_endremoved = fwdvelocity_segments_optoon_extra(1:end-1);
            rotvel_segments_optoon_extra_endremoved = rotvel_segments_optoon_extra(1:end-1);
            min_length = min(cellfun(@length, fwdvelocity_segments_optoon_extra_endremoved));
            rotvel_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), rotvel_segments_optoon_extra_endremoved, 'UniformOutput', false);
            fwdvelocity_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), fwdvelocity_segments_optoon_extra_endremoved, 'UniformOutput', false);
        else
            fwdvelocity_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), fwdvelocity_segments_optoon_extra, 'UniformOutput', false);
            % the shortest segmennt
            rotvel_segments_optoon_trimmed = cellfun(@(x) x(1:min_length), rotvel_segments_optoon_extra, 'UniformOutput', false);
            
        end
        
        % Convert to a matrix for calculating the average
        velocity_matrix_optoon = cell2mat(reshape(fwdvelocity_segments_optoon_trimmed, length(fwdvelocity_segments_optoon_trimmed), []));
        average_velocity_optoon = mean(velocity_matrix_optoon, 1);
        
        
        % Convert to a matrix for calculating the average
        speed_matrix_optoon = cell2mat(reshape(rotvel_segments_optoon_trimmed, length(rotvel_segments_optoon_trimmed), []));
        speed_matrix_optoon = abs(speed_matrix_optoon);
        average_speed_optoon = mean(speed_matrix_optoon, 1);

        % Plot the "optoon" segments with pre-opto points and vertical lines
        figure;
        
        % Subplot 1: Velocity
        subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
        hold on;
        % Plot each raw "optoon" velocity segment in lighter green
        for i = 1:length(fwdvelocity_segments_optoon_extra)
            plot(time(1:length(fwdvelocity_segments_optoon_extra{i})), fwdvelocity_segments_optoon_extra{i}, 'Color', [0.5 1 0.5], 'LineWidth', 1);  % Lighter green
        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);
        
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(average_velocity_optoon)), average_velocity_optoon, 'g', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
        % Customize the first subplot
        xlabel('Time (s)');
        ylabel('Forward Velocity (mm/s)');
        title('Opto-On Forward Velocity');
        legend('Average Opto-On Forward Velocity');
        hold off;
        
        % Subplot 2: Speed
        subplot(2, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
        hold on;
        % Plot each raw "optoon" speed segment in lighter blue
        for i = 1:length(rotvel_segments_optoon_trimmed)
            plot(time(1:length(rotvel_segments_optoon_trimmed{i})), abs(rotvel_segments_optoon_trimmed{i}), 'Color', [0.5 0.5 1], 'LineWidth', 1);  % Lighter blue
        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);
        
        % Calculate and plot the average speed for "optoon" in bold blue
        plot(time(1:length(rotvel_segments_optoon_trimmed{1})), average_speed_optoon, 'b', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Speed');
        
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Speed');
        title('Opto-On Angular Speed');
        legend('Average Opto-On Angular Speed');
        hold off;
        save_plot_with_title_as_filename('forward', 'rotation', savepath)

        %% now separate walking vs not
        still_pre_opto = {}; 
        moving_pre_opto ={};
        ave_still = {}; 
        ave_moving ={};
        still_pre_opto_rot = {}; 
        moving_pre_opto_rot ={};
        ave_still_rot = {}; 
        ave_moving_rot ={};
        moving_opto_forward_acceleration={};
        moving_opto_rot_acceleration={};
        moving_opto_forward={};
        moving_opto_rotation={};
        moving_preopto_forward_acceleration={};
        moving_preopto_rot_acceleration={};
        moving_preopto_forward={};
        moving_preopto_rotation={};
        for i = 1:length(fwdvelocity_segments_optoon_trimmed)
            % pre-speed should be definded as the 1 second prior to
            % stimulation
            opto_points_prespeed_calc = 1 : round(pre_opto_points);
            pre_speed = mean(fwdvelocity_segments_optoon_trimmed{i}(opto_points_prespeed_calc));
            pre_speed_rot = mean(rotvel_segments_optoon_trimmed{i}(opto_points_prespeed_calc));
            if pre_speed <minVel
                %disp("Still")
                still_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_still{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                still_pre_opto_rot{end+1} = rotvel_segments_optoon_extra{i};
                ave_still_rot{end+1} = abs(rotvel_segments_optoon_trimmed{i});
                
            else
                moving_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_moving{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                moving_pre_opto_rot{end+1} = rotvel_segments_optoon_extra{i};
                ave_moving_rot{end+1} = abs(rotvel_segments_optoon_trimmed{i});
                %disp("moving")'
                % NOW LOOK AT ACCELERATION AND SPEED FOR MOVING ONLY TRIALS
                moving_opto_forward_acceleration{end+1} = opto_averages_fwdaccel(:,i);
                moving_opto_rot_acceleration{end+1} = opto_averages_rotaccel(:,i);
                moving_opto_forward{end+1} = opto_averages_fwd(:,i);
                moving_opto_rotation{end+1} = opto_averages_rot(:,i);
                moving_preopto_forward_acceleration{end+1} = no_averages_fwdaccel(:,i);
                moving_preopto_rot_acceleration{end+1} = no_averages_rotaccel(:,i);
                moving_preopto_forward{end+1} = no_averages_fwd(:,i);
                moving_preopto_rotation{end+1} = no_averages_rot(:,i);
            end
        end

        % speed
        moving_opto_forward_ave = mean(cell2mat(moving_opto_forward));
        disp(["MOVING OPTO ON FWD SPEED AVE:", num2str(moving_opto_forward_ave)])
        moving_opto_rot_ave = mean(cell2mat(moving_opto_rotation));
        disp(["MOVING OPTO ON ROT SPEED AVE:", num2str(moving_opto_rot_ave)])
        moving_preopto_forward_ave = mean(cell2mat(moving_preopto_forward));
        disp(["MOVING PRE-OPTO FWD SPEED AVE:", num2str(moving_preopto_forward_ave)])
        moving_preopto_rot_ave = mean(cell2mat(moving_preopto_rotation));
        disp(["MOVING PRE-OPTO ROT SPEED AVE:", num2str(moving_preopto_rot_ave)])


        % acceleration
        moving_opto_forward_acceleration_ave = mean(cell2mat(moving_opto_forward_acceleration));
        disp(["MOVING OPTO ON FWD ACCELERATION AVE:", num2str(moving_opto_forward_acceleration_ave)])
        moving_opto_rot_acceleration_ave = mean(cell2mat(moving_opto_rot_acceleration));
        disp(["MOVING OPTO ON ROT ACCELERATION AVE:", num2str(moving_opto_rot_acceleration_ave)])
        moving_preopto_forward_acceleration_ave = mean(cell2mat(moving_preopto_forward_acceleration));
        disp(["MOVING PRE-OPTO FWD ACCELERATION AVE:", num2str(moving_preopto_forward_acceleration_ave)])
        moving_preopto_rot_acceleration_ave = mean(cell2mat(moving_preopto_rot_acceleration));
        disp(["MOVING PRE-OPTO ROT ACCELERATION AVE:", num2str(moving_preopto_rot_acceleration_ave)])

        percent_moving = length(moving_pre_opto) / (length(still_pre_opto)+length(moving_pre_opto))*100;
        disp(["MOVING PERCENTAGE:", num2str(percent_moving)])

        % take average fwd
        % Convert to a matrix for calculating the average
        speed_matrix_optoon = cell2mat(reshape(rotvel_segments_optoon_trimmed, length(rotvel_segments_optoon_trimmed), []));  
        moving = cell2mat(reshape(ave_moving, length(ave_moving), []));
        moving_ave = mean(moving, 1);
        % Convert to a matrix for calculating the average
        still = cell2mat(reshape(ave_still, length(ave_still), []));
        still_ave = mean(still, 1);
        % take average rot
                % Convert to a matrix for calculating the average
        moving_rot = cell2mat(reshape(ave_moving_rot, length(ave_moving_rot), []));
        moving_ave_rot = mean(moving_rot, 1);
                % Convert to a matrix for calculating the average
        still_rot = cell2mat(reshape(ave_still_rot, length(ave_still_rot), []));
        still_ave_rot = mean(still_rot, 1);


        %% plot this rotation
         % Plot the "optoon" segments with pre-opto points and vertical lines
        figure;
        
        % Subplot 1: moving
        subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
        hold on;
        % Plot each raw "optoon" velocity segment in lighter green
        for i = 1:length(moving_pre_opto_rot)
            plot(time(1:length(moving_pre_opto_rot{i})), abs(moving_pre_opto_rot{i}), 'Color', [0.8, 0.3, 0], 'LineWidth', 1);  % Lighter green

        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);

        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(moving_ave_rot)), moving_ave_rot, 'Color', [1, 0.5, 0], 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
        
        % Customize the first subplot
        xlabel('Time (s)');
        ylabel('Rotational speed (deg/s)');
        title('Moving pre-opto');
        hold off;
   
        % Subplot 2: Speed
        subplot(2, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
        hold on;
        % Plot each raw "optoon" speed segment in lighter blue
        for i = 1:length(still_pre_opto_rot)
            plot(time(1:length(still_pre_opto_rot{i})), abs(still_pre_opto_rot{i}), 'Color', [1 0.5 0.5], 'LineWidth', 1);  % Lighter blue
        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(still_ave_rot)), still_ave_rot, 'r', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
   
       
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Rotational Speed (deg/s)');
        title('Still pre-opto');
        hold off;
        save_plot_with_title_as_filename('rotation_stop', 'rotation_start', savepath)
        %% plot this forward
         % Plot the "optoon" segments with pre-opto points and vertical lines
        figure;
        
        % Subplot 1: moving
        subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
        hold on;
        % Plot each raw "optoon" velocity segment in lighter green
        for i = 1:length(moving_pre_opto)
            plot(time(1:length(moving_pre_opto{i})), moving_pre_opto{i}, 'Color', [0.75, 0, 0.75], 'LineWidth', 1);  % Lighter green
        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);

        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(moving_ave)), moving_ave, 'Color', [1 0 1], 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');

        
        % Customize the first subplot
        xlabel('Time (s)');
        ylabel('Forward Velocity (mm/s)');
        title('Moving pre-opto');
        hold off;
        
        % Subplot 2: Speed
        subplot(2, 1, 2);  % Create the second subplot (2 rows, 1 column, position 2)
        hold on;
        % Plot each raw "optoon" speed segment in lighter blue
        for i = 1:length(still_pre_opto)
            plot(time(1:length(still_pre_opto{i})), still_pre_opto{i}, 'Color', [1 0.5 0.5], 'LineWidth', 1);  % Lighter blue
        end
        xline(time(round(opto_on_markers(1))), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
        xline(time(round(opto_off_markers(1))), 'k--', 'LineWidth', 1.5);
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(still_ave)), still_ave, 'r', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
       
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Forward Velocity (mm/s)');
        title('Still pre-opto');
        hold off;
        save_plot_with_title_as_filename('fwd_stop', 'fwd_start', savepath)
    end

    end