%% This script is an amalgamation of all of the calcium data processing scripts:

% This first function re-orients the G4 panels to the true front (rather
% than where the panels decide is 0) 0 is in front of the fly
% This function also downsamples time, opto, and the panels such that they
% have the same length as the fictrac variables
% This function also uses the functions compute_absolute_circular_diff and 
% detect_local_peaks to determine where jumps are and add in a 2 second
% bloc on either side to remove following a bar jump
% [exptData, exptMeta] = process_fictrac_panels(exptData,exptMeta)

% This function plots the forward and rotational velocities of the fly
% chopped up by the specific optostim going on
% [exptData, exptMeta] = plotExpt_sliced_opto_velocity(expData, exptMeta, savepath);

close all
clear all

%% initialize variables
neck = 0;
ves041 = 0;
an4 = 0;
jump = 0;
lplc4=0;
visual = 0;

%% load in daq and ts
% Define the base folder path
baseFolder = '/Users/sophiarenauld/stacks/'; % Change this to the desired directory

% Prompt the user to select a file from the base folder
[fileName, pathName] = uigetfile('*daq_.mat', 'Select Processed Behavior Data', baseFolder);
[fileName2, pathName2] = uigetfile('*roi_.mat', 'Select Processed ROI Data', baseFolder);
% Check if the user clicked "Cancel"
if isequal(fileName, 0)
    disp('User canceled the file selection.');
else
    % Load the selected file into the workspace
    fullFilePath = fullfile(pathName, fileName);
    fullFilePath2 = fullfile(pathName2, fileName2);
    load(fullFilePath);
    load(fullFilePath2);
    disp(['User selected "', fullFilePath, '".']);
   
end

%% create savepath and folder directory
% Trim the file name to get the base path.
% We do this by removing everything after the last forward slash.
[base_path, ~, ~] = fileparts(fullFilePath);  % This gives us the path without the filename.

if contains(base_path, 'neck')
    neck = 1;  % Set neck to 1 if "neck" is found
else
    neck = 0;  % Set neck to 0 if "neck" is not found
end
if contains(base_path, 'ves041')
    ves041 = 1;  % Set neck to 1 if "neck" is found
else
    ves041 = 0;  % Set neck to 0 if "neck" is not found
end
if contains(base_path, 'jump')
    jump = 1;  % Set neck to 1 if "neck" is found
else
    jump = 0;  % Set neck to 0 if "neck" is not found
end
if contains(base_path, 'an4')
    an4 = 1;  % Set neck to 1 if "neck" is found
elseif contains(base_path, 'LPLC4')
    an4 = 1;  % Set neck to 1 if "neck" is found
    visual = 1;
elseif contains(base_path, 'LC22')
    an4 = 1;  % Set neck to 1 if "neck" is found
    visual = 1;
elseif contains(base_path, 'cmllp01')
    an4 = 1;  % Set neck to 1 if "neck" is found
else
    an4 = 0;  % Set neck to 0 if "neck" is not found
end
% Define the directory for the stimulus
stim_directory = fullfile(base_path, 'smr_analysis');

% Create the directory if it doesn't exist
if ~exist(stim_directory, 'dir')
    mkdir(stim_directory);
end

% Set the savepath variable
savepath = stim_directory;

% Display the save path
disp(['Save path: ', savepath]);

%% process and pick out jumps
minVel = 2;
forvel_cutoff = 0.5;
daq = process_fictrac_panels_2p(daq, minVel, jump);

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
%% now here i will have a function to extract saccades
[daq, saccade_both] = find_saccades_2p(daq, yaw_information_right,yaw_information_left, forvel_cutoff, 0);
[daq, saccade_both_supp] = find_saccades_2p(daq, yaw_information_right_supp,yaw_information_left_supp, forvel_cutoff, 1);
%% plot basic variables
if ves041
    %% basic velocity and dff plotting
    plotting_ves041(daq,ts, jump, savepath)
    %% test turns - rotational velocity
    % Create a figure for plotting
    figure;
    hold on;
    
    plot(daq.t_supp, daq.smoothedangularVelocity_supp)
    
    % Logical indices for turning states
    turning_indices_on = daq.turning_supp == 1;  % Where turning_supp is 1
    turning_indices_off = daq.turning_supp == 0;  % Where turning_supp is 0
    saccading_indices_on = daq.saccading_supp == 1;  % Where turning_supp is 1
    saccading_indices_off = daq.saccading_supp == 0;  % Where turning_supp is 0
    
    % Extract the data for turning on and off
    time_on = daq.t_supp(turning_indices_on);  % Time values where turning_supp is 1
    angular_velocity_on = daq.smoothedangularVelocity_supp(turning_indices_on);  % Corresponding angular velocity values
    for_velocity_on = daq.smoothedfwdVelocity_supp(turning_indices_on);  % Corresponding angular velocity values
    
    
    time_on_saccade = daq.t_supp(saccading_indices_on);  % Time values where turning_supp is 1
    angular_velocity_saccade = daq.smoothedangularVelocity_supp(saccading_indices_on);  % Corresponding angular velocity values
    for_velocity_saccade = daq.smoothedfwdVelocity_supp(saccading_indices_on);
    
    % Plot angular velocity when turning_supp is 1 (green dots)
    scatter(time_on, angular_velocity_on, 10, 'g', 'filled', 'DisplayName', 'Non saccade turn');  % Green dots
    
    
    % Plot angular velocity when saccading_supp is 1 (red dots)
    scatter(time_on_saccade, angular_velocity_saccade, 10, 'r', 'filled', 'DisplayName', 'Saccade');  % Red dots
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('Smoothed Angular Velocity');
    title('Smoothed Angular Velocity During Turning States (Dots)');
    grid on;
    
    % Add a legend to distinguish between states
    legend('show', 'Location', 'Best');  % Displays the legends based on 'DisplayName'
    
    hold off;
    save_plot_with_title_as_filename('rv', 'turns_v_saccades', savepath)
    
    %% test turns - forward
    figure
    hold on
    
    plot(daq.t_supp, daq.smoothedfwdVelocity_supp)
    % Plot angular velocity when turning_supp is 1 (green dots)
    scatter(time_on, for_velocity_on, 10, 'g', 'filled', 'DisplayName', 'Turning On');  % Green dots
    
    
    % Plot angular velocity when saccading_supp is 1 (red dots)
    scatter(time_on_saccade, for_velocity_saccade, 10, 'r', 'filled', 'DisplayName', 'Turning Off');  % Red dots
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('Smoothed Angular Velocity');
    title('Smoothed Forward Velocity During Turning States (Dots)');
    grid on;
    
    % Add a legend to distinguish between states
    legend('show', 'Location', 'Best');  % Displays the legends based on 'DisplayName'
    
    hold off;
    save_plot_with_title_as_filename('fwdv', 'turns_vs_saccades', savepath)
    
    
    %% upsample for dff with turning / saccades
    % Assume ts{1,1} and daq.bfv are defined
    ts_data = ts{1,1};                       % Your original data (1x2000)
    target_length = length(daq.bfv_supp);        % Length of target data (1x18000)
    
    % Create a time vector for original data
    original_length = length(ts_data);
    original_time = daq.t; % Example normalization to unit length (you can adjust based on actual time)
    
    % Create a time vector for the target length
    target_time = daq.t_supp;   % Corresponding time for target data
    
    % Upsample ts_data using interpolation
    upsampled_ts_data = interp1(original_time, ts_data, target_time, 'linear'); % Linear interpolation
    upsampled_x_data = interp1(original_time, daq.px, target_time, 'linear'); % Linear interpolation
    upsampled_y_data = interp1(original_time, daq.py, target_time, 'linear'); % Linear interpolation
    %% plot data for both all turns and saccades
    plot_dff_saccades(daq, upsampled_ts_data, yaw_both_supp, savepath, 1)
    plot_dff_saccades(daq, upsampled_ts_data, saccade_both_supp, savepath, 0)
    
    plot_flatpath_saccades(daq, jump, savepath, daq.px_supp, daq.py_supp, 1);
    plot_flatpath_saccades(daq, jump, savepath, daq.px_supp, daq.py_supp, 0);


end

if neck
    plotting_neck(daq, ts, savepath)
end

if an4
    plotting_an4(daq, ts, jump, savepath)
    lm_velocity(daq, ts, savepath)
    
end
if visual
    openloop_gratings_process(daq, ts, savepath)
end


%% PILOT ZONE
    %% plot flat path with saccades
    %plot_flatpath_saccades(daq, jump,savepath);
    %%% plot heatmap of forward acceleration and rs and dff
    % ts.ball.foraccel = diff(ts.ball.forvel);
    % ts.ball.yawvel_adjusted = ts.ball.yawdeg(1:end-1); 
    % ts.ball.sidevel_adjusted = ts.ball.sidevel(1:end-1);
    % ts.ball.dffvel_adjusted = ts.resp.i2{1}(1:end-1); 
    % binnedData = dff_v_velocityHeat(ts.ball.foraccel,abs(ts.ball.yawvel_adjusted),abs(ts.ball.sidevel_adjusted),ts.ball.dffvel_adjusted,1);
    %% plot dff split up by saccades 
    %plot_dff_saccades(daq, ts, saccade_both)
