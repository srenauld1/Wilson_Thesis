function [daq, saccade_both] = find_saccades_2p(daq, yaw_information_right,yaw_information_left,forvel_cutoff, upsample)
% 
%   This function takes in information on turning (from findYawPeaksFT.m)
%   and returns only those turns that qualify as saccades. this means that
%   the forward velocity at the peak of the turn is less than the forward
%   velocity at the beginning of the turn (less by forvel_cutoff)
%
%
% INPUTS:
%   daq - struct with kinematic data (both imaging frame rate and 60hz)
%   yaw_information_right - struct from findYawVelPeaksFT.m for extracting
%       all right turns in a trial
%   yaw_information_left - struct from findYawVelPeaksFT.m for extracting
%       all left turns in a trial
%   forvel_cutoff - velocity value set in mater_imaging_analysis to
%       determine the amount slow down required to be a saccade (1 is a good
%       start)
%   upsample - boolean for whether imaging rate (0) or 60Hz fictrac (1)
%
% OUTPUTS:
%   daq - returns kinematic data with structs added (?)
%   saccade_both - a struct containing information from yaw_information
%   left and right for only turns defined as saccades
%
%
% LAST UPDATED
%   3/25 - SMR CREATED
%


% initialize kinematic data
if (upsample)
    forward = daq.smoothedfwdVelocity_supp;
    yaw = daq.bvy_deg_supp;
    time = daq.t_supp;
    % Initialize the saccading array with zeros
    numEntries = size(time,2);  % The total number of observations you want to evaluate
    daq.saccading_supp = zeros(1, numEntries);  % Initialize as non-saccading
else
    forward = daq.smoothedfwdVelocity;
    yaw = daq.bvy_deg;
    time = daq.t;
    % Initialize the saccading array with zeros
    numEntries = size(time,2);  % The total number of observations you want to evaluate
    daq.saccading = zeros(1, numEntries);  % Initialize as non-saccading
end


%% left saccades

% Initialize an array to hold the differences
numPairs = min(length(yaw_information_left.boutStartInd), length(yaw_information_left.yawVelPeakInd));  % To avoid index out of bounds
velocityDifferences_left = zeros(numPairs, 1);  % Preallocate for efficiency

% Calculate the velocity difference between the start and peak indices
for i = 1:numPairs
    startIdx = yaw_information_left.boutStartInd(i);
    peakIdx = yaw_information_left.yawVelPeakInd(i);
    
    % Ensure the indices are within bounds of the forwardVelocity array
    if startIdx <= length(forward) && peakIdx <= length(forward)
        velocityDifferences_left(i) = forward(peakIdx) - forward(startIdx);
    else
        warning('Index out of bounds for forwardVelocity at pair %d', i);
        velocityDifferences_left(i) = NaN;  % Use NaN for out-of-bounds
    end
end

% Display the differences
% disp('Differences in forward velocity at start and peak indices:');
% disp(velocityDifferences_left);

% threshold differences for a decrease in forward velocity
thresholded_velocities_left = zeros(size(velocityDifferences_left));  % Initialize an array of zeros

% Apply the thresholding
thresholded_velocities_left(velocityDifferences_left >= -forvel_cutoff) = 0;  % Set to 0 if >= -forvel_cutoff
thresholded_velocities_left(velocityDifferences_left < -forvel_cutoff) = 1;   % Set to 1 if < -forvel_cutoff

% Initialize a new structure to hold the filtered yaw information
saccade_left = struct();

% Filter each field based on the threshold condition
saccade_left.yawVelPeakTimes = yaw_information_left.yawVelPeakTimes(thresholded_velocities_left == 1);
saccade_left.boutStartTimes = yaw_information_left.boutStartTimes(thresholded_velocities_left == 1);
saccade_left.boutEndTimes = yaw_information_left.boutEndTimes(thresholded_velocities_left == 1);
saccade_left.yawVelPeakInd = yaw_information_left.yawVelPeakInd(thresholded_velocities_left == 1);
saccade_left.boutStartInd = yaw_information_left.boutStartInd(thresholded_velocities_left == 1);
saccade_left.boutEndInd = yaw_information_left.boutEndInd(thresholded_velocities_left == 1);


%% right saccades

% Initialize an array to hold the differences
numPairs = min(length(yaw_information_right.boutStartInd), length(yaw_information_right.yawVelPeakInd));  % To avoid index out of bounds
velocityDifferences_right = zeros(numPairs, 1);  % Preallocate for efficiency

% Calculate the velocity difference between the start and peak indices
for i = 1:numPairs
    startIdx = yaw_information_right.boutStartInd(i);
    peakIdx = yaw_information_right.yawVelPeakInd(i);
    
    % Ensure the indices are within bounds of the forwardVelocity array
    if startIdx <= length(forward) && peakIdx <= length(forward)
        velocityDifferences_right(i) = forward(peakIdx) - forward(startIdx);
    else
        warning('Index out of bounds for forwardVelocity at pair %d', i);
        velocityDifferences_right(i) = NaN;  % Use NaN for out-of-bounds
    end
end

% Display the differences
%disp('Differences in forward velocity at start and peak indices:');
%disp(velocityDifferences_right);

% threshold differences for a decrease in forward velocity
thresholded_velocities_right = zeros(size(velocityDifferences_right));  % Initialize an array of zeros

% Apply the thresholding
thresholded_velocities_right(velocityDifferences_right >= -forvel_cutoff) = 0;  % Set to 0 if >= -forvel_cutoff
thresholded_velocities_right(velocityDifferences_right < -forvel_cutoff) = 1;   % Set to 1 if < -forvel_cutoff

% Initialize a new structure to hold the filtered yaw information
saccade_right = struct();

% Filter each field based on the threshold condition
saccade_right.yawVelPeakTimes = yaw_information_right.yawVelPeakTimes(thresholded_velocities_right == 1);
saccade_right.boutStartTimes = yaw_information_right.boutStartTimes(thresholded_velocities_right == 1);
saccade_right.boutEndTimes = yaw_information_right.boutEndTimes(thresholded_velocities_right == 1);
saccade_right.yawVelPeakInd = yaw_information_right.yawVelPeakInd(thresholded_velocities_right == 1);
saccade_right.boutStartInd = yaw_information_right.boutStartInd(thresholded_velocities_right == 1);
saccade_right.boutEndInd = yaw_information_right.boutEndInd(thresholded_velocities_right == 1);

%% merge left and right
saccade_both = struct();
saccade_both.yawVelPeakTimes = [saccade_left.yawVelPeakTimes, saccade_right.yawVelPeakTimes];
saccade_both.boutStartTimes = [saccade_left.boutStartTimes, saccade_right.boutStartTimes];
saccade_both.boutEndTimes = [saccade_left.boutEndTimes, saccade_right.boutEndTimes];
saccade_both.yawVelPeakInd = [saccade_left.yawVelPeakInd, saccade_right.yawVelPeakInd];
saccade_both.boutStartInd = [saccade_left.boutStartInd, saccade_right.boutStartInd];
saccade_both.boutEndInd = [saccade_left.boutEndInd, saccade_right.boutEndInd];

%% add to expt data

% Loop through each entry for saccades
for i = 1:length(saccade_both.boutStartInd)
    bout_start = saccade_both.boutStartInd(i);
    bout_end = saccade_both.boutEndInd(i);

    % Set saccading to 1 for any index in between start and end
    if bout_start <= numEntries && bout_end <= numEntries
        if (upsample)
            daq.saccading_supp(1, bout_start:bout_end) = 1;  % Mark saccading for left
        else
            daq.saccading(1, bout_start:bout_end) = 1;  % Mark saccading for left
        end
    end
end


end