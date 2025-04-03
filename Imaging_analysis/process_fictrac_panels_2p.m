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


function [ts] = process_fictrac_panels_2p(ts, minVel, jump)
    
    if jump
        %% detect jumps
        ts = compute_absolute_circular_diff_2p(ts);
        ts = detect_local_peaks_2p(ts);
    
        %% add in jump buffer for later calculations
        % Define the padding size on each side
        one_indices = find(ts.vis.jump_detected == 1);
        jump_idx = ts.vis.jump_detected;
    
        % calculate trial duration
        duration = length(ts.t)/ts.t(end);
        
        % Define the padding size on each side
        buffer_before_and_after_jump = 2;
        padding = round(buffer_before_and_after_jump*duration);
        
        % Expand each 1 into a block of ones
        for i = 1:length(one_indices)
            idx = one_indices(i);
            start_idx = max(1, idx - padding);
            end_idx = min(length(ts.vis.jump_detected), idx + padding);
            jump_idx(start_idx:end_idx) = 1;
        end
    
        ts.vis.omit_jump_blocks = jump_idx;
    end

    %% create moving not moving binary
    % calculate total speed
    total_speed = abs(ts.ball.forvel) +abs(ts.ball.yaw) + abs(ts.ball.sidevel);
    ts.ball.totalspeed = total_speed;
    ids_fly_moving = find(ts.ball.totalspeed >= minVel);

    % Create a column for moving/not moving
    ts.motion.moving_not = zeros(size(ts.ball.totalspeed));  % Initialize with zeros
    
    ts.motion.moving_not(ids_fly_moving) = 1;  % Set to 1 where the conditions are met
    
    % Find indices where the fly is not moving (where movement_state is 0)
    ftNotMoveInd = find(ts.motion.moving_not == 0);
    
    ts.motion.ftNotMoveInd = ftNotMoveInd;

    %% turn yaw into degrees
    ts.ball.yawdeg = rad2deg(ts.ball.yawvel);
    
    
    %% overly smooth rotational velocity
    ts.smoothedangularVelocity = smoothdata(ts.ball.yawdeg, 'gaussian',10);
    ts.smoothedfwdVelocity = smoothdata(ts.ball.forvel, 'gaussian',10);


end