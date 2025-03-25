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
%           3/19/2025 - SMR adjusted to only plot closed loop data

function [exptData, exptMeta] = plotExpt_sliced_opto_velocity(exptData,exptMeta, savepath)
    
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
    on_off = 0;

    %% define variables
    opto=exptData.optoStim;
    time = exptData.t;

    %% chop up by opto blocks
    if checkOpto
        % Find the change points in the opto signal
        change_points = [1; find(diff(exptData.optoStim) ~= 0) + 1; length(opto)];
        % Initialize cell arrays to store segments of velocity
        velocity_segments_optoon = {};  % Velocity when opto is 1
        velocity_segments_optooff = {};  % Velocity when opto is 0
        % now take care off immediate on off for 10 sec stim windows
        if sum(opto)==36000
            % Assuming your array is called change_points
            keep_idx = false(size(change_points));
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
        disp(["OPTO ON FWD AVE:", num2str(average_fwd)])
        average_rot = mean(opto_averages_rot);
        disp(["OPTO ON ROT AVE:", num2str(average_rot)])
        average_fwd_no = mean(no_averages_fwd);
        disp(["OPTO OFF FWD AVE:", num2str(average_fwd_no)])
        average_rot_no = mean(no_averages_rot);
        disp(["OPTO OFF ROT AVE:", num2str(average_rot_no)])

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
        velocity_matrix_optoon = cell2mat(reshape(fwdvelocity_segments_optoon_trimmed, 1, []));
        average_velocity_optoon = mean(velocity_matrix_optoon, 2);
        
        % Ensure all "optoon" segments are the same length by trimming to
        
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

        %% now separate walking vs not
        still_pre_opto = {}; 
        moving_pre_opto ={};
        ave_still = {}; 
        ave_moving ={};
        still_pre_opto_rot = {}; 
        moving_pre_opto_rot ={};
        ave_still_rot = {}; 
        ave_moving_rot ={};
        for i = 1:length(fwdvelocity_segments_optoon_trimmed)
            pre_speed = mean(fwdvelocity_segments_optoon_trimmed{i}((pre_opto_points-100):pre_opto_points));
            if on_off
                pre_speed = mean(fwdvelocity_segments_optoon_trimmed{i}((pre_opto_points+2900):pre_opto_points+3000));
            end
            pre_speed_rot = mean(rotspeed_segments_optoon_trimmed{i}((pre_opto_points-100):pre_opto_points));
            if pre_speed <3
                disp("Still")
                still_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_still{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                still_pre_opto_rot{end+1} = rotspeed_segments_optoon_extra{i};
                ave_still_rot{end+1} = rotspeed_segments_optoon_trimmed{i};
                
            else
                moving_pre_opto{end+1} = fwdvelocity_segments_optoon_extra{i};
                ave_moving{end+1} = fwdvelocity_segments_optoon_trimmed{i};
                moving_pre_opto_rot{end+1} = rotspeed_segments_optoon_extra{i};
                ave_moving_rot{end+1} = rotspeed_segments_optoon_trimmed{i};
                disp("moving")
            end
        end

        % take average fwd
                % Convert to a matrix for calculating the average
        moving = cell2mat(reshape(ave_moving, 1, []));
        moving_ave = mean(moving, 2);
                % Convert to a matrix for calculating the average
        still = cell2mat(reshape(ave_still, 1, []));
        still_ave = mean(still, 2);
        % take average rot
                % Convert to a matrix for calculating the average
        moving_rot = cell2mat(reshape(ave_moving_rot, 1, []));
        moving_ave_rot = mean(moving_rot, 2);
                % Convert to a matrix for calculating the average
        still_rot = cell2mat(reshape(ave_still_rot, 1, []));
        still_ave_rot = mean(still_rot, 2);


        %% plot this rotation
         % Plot the "optoon" segments with pre-opto points and vertical lines
        figure;
        
        % Subplot 1: moving
        subplot(2, 1, 1);  % Create the first subplot (2 rows, 1 column, position 1)
        hold on;
        % Plot each raw "optoon" velocity segment in lighter green
        for i = 1:length(moving_pre_opto_rot)
            plot(time(1:length(moving_pre_opto_rot{i})), moving_pre_opto_rot{i}, 'Color', [0.75, 0, 0.75], 'LineWidth', 1);  % Lighter green
            xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
            xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
        end

        if on_off
            xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
            xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
        end

        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(moving_ave_rot)), moving_ave_rot, 'Color', [1 0 1], 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
        
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
            plot(time(1:length(still_pre_opto_rot{i})), still_pre_opto_rot{i}, 'Color', [1 0.5 0.5], 'LineWidth', 1);  % Lighter blue
            xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
            xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
        end
        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(still_ave_rot)), still_ave_rot, 'r', 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        
        if on_off
            xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
            xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
        end
       
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Rotational Speed (deg/s)');
        title('Still pre-opto');
        hold off;
        %saveas(gcf, ['rs_fw_by_opto.png']);
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
            xline(time(opto_on_markers(1)), 'k--', 'LineWidth', 1.5);  % Dashed black line at opto-on marker
            xline(time(opto_off_markers(1)), 'k--', 'LineWidth', 1.5);
        end

        % Plot the average velocity for "optoon" in bold green
        plot(time(1:length(moving_ave)), moving_ave, 'Color', [1 0 1], 'LineWidth', 3, 'DisplayName', 'Average Opto-On Velocity');
        if on_off
            xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
            xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
        end
        
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
        if on_off
            xline(time((3000+opto_on_markers(1))), 'k--', 'LineWidth', 1.5)
            xline(time((opto_off_markers(1)-3000)), 'k--', 'LineWidth', 1.5)
        end
       
        % Customize the second subplot
        xlabel('Time (s)');
        ylabel('Forward Velocity (mm/s)');
        title('Still pre-opto');
        hold off;
        save_plot_with_title_as_filename('fwd_stop', 'fwd_start', savepath)
    end

    end