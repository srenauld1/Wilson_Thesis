% process_fictrac_panels.m
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


function [exptData, exptMeta] = process_fictrac_panels(exptData,exptMeta, minVel)
    
    % reset counters
    n = 0;
    
    % set number of sublots based on number of expt variables
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
            exptData.g4displayXPos = g4Pos_mod;
        end
    end
    %% Behavior variable adjustments
    %downsample time and optoopto and g4
    exptData.t = downsample(exptData.t, 30);
    exptData.optoStim = downsample(exptData.optoStim, 30);
    exptData.g4displayXPos = downsample(exptData.g4displayXPos, 30);

    %% detect jumps
    exptData = compute_absolute_circular_diff(exptData);
    exptData = detect_local_peaks(exptData);

    %% add in jump buffer for later calculations
    % Define the padding size on each side
    one_indices = find(exptData.jump_detected == 1);
    jump_idx = exptData.jump_detected;
    
    % Define the padding size on each side
    buffer_before_and_after_jump = 2;
    padding = buffer_before_and_after_jump*60;
    
    % Expand each 1 into a block of ones
    for i = 1:length(one_indices)
        idx = one_indices(i);
        start_idx = max(1, idx - padding);
        end_idx = min(length(exptData.jump_detected), idx + padding);
        jump_idx(start_idx:end_idx) = 1;
    end

    exptData.omit_jump_blocks = jump_idx;

    %% create moving not moving binary
    ids_fly_moving = find(exptData.totSpeed >= minVel);

    % Create a column for moving/not moving
    exptData.motion.moving_not = zeros(size(exptData.totSpeed));  % Initialize with zeros
    
    exptData.motion.moving_not(ids_fly_moving) = 1;  % Set to 1 where the conditions are met
    
    % Find indices where the fly is not moving (where movement_state is 0)
    ftNotMoveInd = find(exptData.motion.moving_not == 0);
    
    exptData.motion.ftNotMoveInd = ftNotMoveInd;
    
    
    %% overly smooth rotational velocity
    exptData.smoothedangularVelocity = smoothdata(exptData.angularVelocity, 'gaussian',500);

    %% flip everything to be rows
    % Initialize a new structure for the flipped data
    flippedStruct = struct();
    
    % Get the field names of the original structure
    fieldNames = fieldnames(exptData);
    
    % Loop through each field and transpose the data
    for i = 1:length(fieldNames)
        field = fieldNames{i};              % Get the current field name
        flippedStruct.(field) = exptData.(field)';  % Transpose the data to a row vector
    end
    exptData = flippedStruct;


end