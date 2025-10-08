% process_fictrac_panels_2p.m
%
% Function that adjusts the exptData from behavior trials on the behavior
% patio to make sure all things are the same dimentions, re-center the G4
% panels, detect jumps, and overly smooths fictrac data for finding peaks
%
% INPUTS:
%   exptData - processed data, can contain ephys, fictrac, and/or output
%   exptMeta - processed ephys meta
%
%
% OUTPUTS:
%   returns exptData and exptMeta with adjusted variables
%   Variables:
%       exptData.g4displayXPos adjusted to have 0 be in front of fly
%       exptData.t, optostim, g4panes all downsampled by 30
%       exptData.jump_detected, the exact frame where a jump occurs
%       exptData.omit_jump_blocks, a matrix of 0s and 1s showing when to
%           heading data +/- 2 seconds after a jump
%       exptData.motion.moving_not, a matrix of 0s and 1s where 1 is moving, 0 is
%           not moving
%       exptData.motion.ftNotMoveInd, inidices where the fly isnt moving
%       exptData.smoothedangularVelocity, excessively smoothed angular
%           velocity to pick out turning bouts
%
% Original: 3/19/2025 - SMR


function [daq, yaw_information_right,yaw_information_left,yaw_information_left_supp, yaw_information_right_supp ] = process_fictrac_panels_2p(daq, vel_on,vel_off, jump)

    if jump
        %% detect jumps
        daq = compute_absolute_circular_diff_2p(daq);
        daq = detect_local_peaks_2p(daq);
    
        %% add in jump buffer for later calculations
        % Define the padding size on each side
        one_indices = find(daq.jump_detected == 1);
        jump_idx = daq.jump_detected;
    
        % calculate trial duration
        duration = length(daq.t)/daq.t(end);
        
        % Define the padding size on each side
        buffer_before_and_after_jump = 2;
        padding = round(buffer_before_and_after_jump*duration);
        
        % Expand each 1 into a block of ones
        for i = 1:length(one_indices)
            idx = one_indices(i);
            start_idx = max(1, idx - padding);
            end_idx = min(length(daq.jump_detected), idx + padding);
            jump_idx(start_idx:end_idx) = 1;
        end
    
        daq.omit_jump_blocks = jump_idx;
    end

    %% remove anything that is empty
    % Get the field names of the structure
    fieldNames = fieldnames(daq);
    
    % Loop through each field to check for emptiness
    for i = 1:length(fieldNames)
        field = fieldNames{i};  % Current field name
        
        % Check if the current field is empty
        if isempty(daq.(field))
            daq = rmfield(daq, field);  % Remove the field if it's empty
        end
    end

    %% Define Schmitt trigger thresholds (example: you choose these values!)
    %speed = abs(daq.bvf) + (abs(daq.bvy)*4.5) + (abs(daq.bvs)*4.5);
    %daq.totalspeed = speed;
    speed = daq.bvf;
    
    motion = zeros(size(speed));
    
    % Start state (assume not moving)
    if speed(1) >= vel_on
        motion(1) = 1;
    else
        motion(1) = 0;
    end
    
    % Apply Schmitt trigger logic
    for i = 2:numel(speed)
        if motion(i-1) == 0 && speed(i) >= vel_on
            motion(i) = 1; % start moving
        elseif motion(i-1) == 1 && speed(i) > vel_off
            motion(i) = 1; % keep moving
        else
            motion(i) = 0; % stop moving
        end
    end
    
    daq.motion.moving_or_not = motion;
    
    % Optionally: Indices where fly is NOT moving
    daq.motion.ftNotMoveInd = find(motion==0);

    %% turn yaw into degrees
    daq.bvy_deg = rad2deg(daq.bvy);
    daq.bvy_deg_supp = rad2deg(daq.bvy_supp);
    daq.bvs_deg = rad2deg(daq.bvs);
    daq.bvs_deg_supp = rad2deg(daq.bvs_supp);
    
    
    %% overly smooth rotational velocity
    daq.smoothedangularVelocity = smoothdata(daq.bvy_deg, 'gaussian',30);
    daq.smoothedfwdVelocity = smoothdata(daq.bvf, 'gaussian',30);

    %% create moving not moving binary for 60Hz data
    % Compute the total speed as before
    total_speed_supp = abs(daq.bvf_supp) + (4.5*abs(daq.bvy_supp)) + (4.5*abs(daq.bvs_supp));
    daq.totalspeed_supp = total_speed_supp;
    
    % Initialize movement state
    motion_supp = zeros(size(total_speed_supp));
    
    % Start state (can also force to zero if you prefer)
    if total_speed_supp(1) >= vel_on
        motion_supp(1) = 1;
    else
        motion_supp(1) = 0;
    end
    
    % Schmitt trigger logic
    for i = 2:numel(total_speed_supp)
        if motion_supp(i-1) == 0 && total_speed_supp(i) >= vel_on
            motion_supp(i) = 1; % Enter moving state
        elseif motion_supp(i-1) == 1 && total_speed_supp(i) > vel_off
            motion_supp(i) = 1; % Remain moving
        else
            motion_supp(i) = 0; % Enter not-moving state
        end
    end
    
    % Store results
    daq.motion_supp.moving_not = motion_supp;
    daq.motion_supp.ftNotMoveInd = find(motion_supp==0);

    %% turn yaw into degrees
    daq.bvy_deg_supp = rad2deg(daq.bvy_supp);
    
    
    %% overly smooth rotational velocity
    daq.smoothedangularVelocity_supp = smoothdata(daq.bvy_deg_supp, 'gaussian',100);
    daq.smoothedfwdVelocity_supp = smoothdata(daq.bvf_supp, 'gaussian',100);
    
    %% find turns
    yaw_information_right = findYawVelPeaksFT(daq, 20, [0.2,1], daq.motion, 1,0);
    yaw_information_left = findYawVelPeaksFT(daq, 20, [0.2,1], daq.motion, 0, 0);
    yaw_information_right_supp = findYawVelPeaksFT(daq, 20, [0.2,3], daq.motion_supp, 1, 1);
    yaw_information_left_supp = findYawVelPeaksFT(daq, 20, [0.2,3], daq.motion_supp, 0, 1);
    
    %% merge turn data for assessing algorythm
    yaw_both = struct();
    yaw_both.yawVelPeakTimes = [yaw_information_left.yawVelPeakTimes, yaw_information_right.yawVelPeakTimes];
    yaw_both.boutStartTimes = [yaw_information_left.boutStartTimes, yaw_information_right.boutStartTimes];
    yaw_both.boutEndTimes = [yaw_information_left.boutEndTimes, yaw_information_right.boutEndTimes];
    yaw_both.yawVelPeakInd = [yaw_information_left.yawVelPeakInd, yaw_information_right.yawVelPeakInd];
    yaw_both.boutStartInd = [yaw_information_left.boutStartInd, yaw_information_right.boutStartInd];
    yaw_both.boutEndInd = [yaw_information_left.boutEndInd, yaw_information_right.boutEndInd];
    
    %% add to expt data
    daq.turning = zeros(1, length(daq.smoothedangularVelocity));
    % Loop through each entry for saccades
    for i = 1:length(yaw_both.boutStartInd)
        bout_start = yaw_both.boutStartInd(i);
        bout_end = yaw_both.boutEndInd(i);
    
        % Set saccading to 1 for any index in between start and end
        if bout_start <= length(daq.smoothedangularVelocity) && bout_end <= length(daq.smoothedangularVelocity)
            daq.turning(1, bout_start:bout_end) = 1;  % Mark saccading for left
        end
    end
    
    %% merge turn data for assessing algorythm
    yaw_both_supp = struct();
    yaw_both_supp.yawVelPeakTimes = [yaw_information_left_supp.yawVelPeakTimes, yaw_information_right_supp.yawVelPeakTimes];
    yaw_both_supp.boutStartTimes = [yaw_information_left_supp.boutStartTimes, yaw_information_right_supp.boutStartTimes];
    yaw_both_supp.boutEndTimes = [yaw_information_left_supp.boutEndTimes, yaw_information_right_supp.boutEndTimes];
    yaw_both_supp.yawVelPeakInd = [yaw_information_left_supp.yawVelPeakInd, yaw_information_right_supp.yawVelPeakInd];
    yaw_both_supp.boutStartInd = [yaw_information_left_supp.boutStartInd, yaw_information_right_supp.boutStartInd];
    yaw_both_supp.boutEndInd = [yaw_information_left_supp.boutEndInd, yaw_information_right_supp.boutEndInd];
    
    %% add to expt data
    daq.turning_supp = zeros(1, length(daq.smoothedangularVelocity_supp));
    % Loop through each entry for saccades
    for i = 1:length(yaw_both_supp.boutStartInd)
        bout_start = yaw_both_supp.boutStartInd(i);
        bout_end = yaw_both_supp.boutEndInd(i);
    
        % Set saccading to 1 for any index in between start and end
        if bout_start <= length(daq.smoothedangularVelocity_supp) && bout_end <= length(daq.smoothedangularVelocity_supp)
            daq.turning_supp(1, bout_start:bout_end) = 1;  % Mark saccading for left  % Mark saccading for left
        end
    end

end