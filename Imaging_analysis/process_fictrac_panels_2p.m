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


function [daq] = process_fictrac_panels_2p(daq, minVel, jump)
    
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

    %% create moving not moving binary for downsampled data
    % calculate total speed
    total_speed = abs(daq.bfv) +abs(daq.byv) + abs(daq.bsv);
    daq.totalspeed = total_speed;
    ids_fly_moving = find(daq.totalspeed >= minVel);

    % Create a column for moving/not moving
    daq.motion.moving_not = zeros(size(daq.totalspeed));  % Initialize with zeros
    
    daq.motion.moving_not(ids_fly_moving) = 1;  % Set to 1 where the conditions are met
    
    % Find indices where the fly is not moving (where movement_state is 0)
    ftNotMoveInd = find(daq.motion.moving_not == 0);
    
    daq.motion.ftNotMoveInd = ftNotMoveInd;

    %% turn yaw into degrees
    daq.byv_deg = rad2deg(daq.byv);
    
    
    %% overly smooth rotational velocity
    daq.smoothedangularVelocity = smoothdata(daq.byv_deg, 'gaussian',30);
    daq.smoothedfwdVelocity = smoothdata(daq.bfv, 'gaussian',30);

    %% create moving not moving binary for 60Hz data
    % calculate total speed
    total_speed = abs(daq.bfv_supp) +abs(daq.byv_supp) + abs(daq.bsv_supp);
    daq.totalspeed_supp = total_speed;
    ids_fly_moving = find(daq.totalspeed_supp >= minVel);

    % Create a column for moving/not moving
    daq.motion_supp.moving_not = zeros(size(daq.totalspeed_supp));  % Initialize with zeros
    
    daq.motion_supp.moving_not(ids_fly_moving) = 1;  % Set to 1 where the conditions are met
    
    % Find indices where the fly is not moving (where movement_state is 0)
    ftNotMoveInd = find(daq.motion_supp.moving_not == 0);
    
    daq.motion_supp.ftNotMoveInd = ftNotMoveInd;

    %% turn yaw into degrees
    daq.byv_deg_supp = rad2deg(daq.byv_supp);
    
    
    %% overly smooth rotational velocity
    daq.smoothedangularVelocity_supp = smoothdata(daq.byv_deg_supp, 'gaussian',100);
    daq.smoothedfwdVelocity_supp = smoothdata(daq.bfv_supp, 'gaussian',100);


end