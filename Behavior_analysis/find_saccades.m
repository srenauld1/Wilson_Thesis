function [exptData, saccade_both] = find_saccades(exptData, yaw_information_right,yaw_information_left, forvel_cutoff)
% 
%
% INPUTS:
%
% OUTPUTS:
%

% initialize kinematic data
forward = exptData.forwardVelocity;
yaw = exptData.angularVelocity;


%% left saccades

% Initialize an array to hold the differences
numPairs = min(length(yaw_information_left.boutStartInd), length(yaw_information_left.yawVelPeakInd));  % To avoid index out of bounds
velocityDifferences_left = zeros(1, numPairs);  % Preallocate for efficiency

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
thresholded_velocities_left(velocityDifferences_left >= -forvel_cutoff) = 0;  % Set to 0 if >= -1
thresholded_velocities_left(velocityDifferences_left < -forvel_cutoff) = 1;   % Set to 1 if < -1

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
thresholded_velocities_right(velocityDifferences_right >= -forvel_cutoff) = 0;  % Set to 0 if >= -cutoff
thresholded_velocities_right(velocityDifferences_right < -forvel_cutoff) = 1;   % Set to 1 if < -cutoff

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

% Initialize the saccading array with zeros
numEntries = max(size(exptData.t, 2), size(exptData.t, 1));  % The total number of observations you want to evaluate
exptData.saccading = zeros(1, numEntries);  % Initialize as non-saccading

% Loop through each entry for saccades
for i = 1:length(saccade_both.boutStartInd)
    bout_start = saccade_both.boutStartInd(i);
    bout_end = saccade_both.boutEndInd(i);

    % Set saccading to 1 for any index in between start and end
    if bout_start <= numEntries && bout_end <= numEntries
        exptData.saccading(bout_start:bout_end) = 1;  % Mark saccading for left
    end
end

end