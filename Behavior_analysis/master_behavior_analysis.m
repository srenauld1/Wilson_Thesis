%% This script is an amalgamation of all of the behavior data processing scripts:

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

%% define initial variables
minVel = 3; % minimum total velocity to be count as moving
forvel_cutoff = 1;

%% load in exptData and exptMeta
% Define the base folder path
baseFolder = '/Volumes/Neurobio/wilsonlab/Sophia/behavior/'; % Change this to the desired directory

% Prompt the user to select a file from the base folder
[fileName, pathName] = uigetfile('*.mat', 'Select Processed Behavior Data', baseFolder);

% Check if the user clicked "Cancel"
if isequal(fileName, 0)
    disp('User canceled the file selection.');
else
    % Load the selected file into the workspace
    fullFilePath = fullfile(pathName, fileName);
    load(fullFilePath);
    disp(['User selected "', fullFilePath, '".']);
   
end

%% create savepath and folder directory
% Trim the file name to get the base path.
% We do this by removing everything after the last forward slash.
[base_path, ~, ~] = fileparts(fullFilePath);  % This gives us the path without the filename.

% Extract the stimulus number from the file name.
% Assuming the last part contains the stim part (e.g. "stim2"), we split by spaces:
parts = strsplit(fullFilePath, ' ');

% Get the last element that contains 'stimX' (stimulus info)
stim_part = parts{end-1};  % This gets 'stim2' or similar
stim_number = extractAfter(stim_part, 'stim');  % Extract the number, giving us '2'
% get the trial number and letter
trial_part = parts{1};
trial_num = extractAfter(trial_part, "trial");
% Define the directory for the stimulus
stim_directory = fullfile(base_path, ['trial', trial_num,'stim', stim_number]);

% Create the directory if it doesn't exist
if ~exist(stim_directory, 'dir')
    mkdir(stim_directory);
end

% Set the savepath variable
savepath = stim_directory;

% Display the save path
disp(['Save path: ', savepath]);

%% process experimental data
[exptData, exptMeta] = process_fictrac_panels(exptData,exptMeta, minVel);

%% plot velocity data split up by opto chunks
[exptData,  exptMeta, fwdvelocity_segments_optoon_trimmed, rotvel_segments_optoon_trimmed] = plotExpt_sliced_opto_velocity(exptData, exptMeta, savepath);

%% determine turn direction
[exptData,  exptMeta, fwdvelocity_segments_optoon_trimmed, rotvel_segments_optoon_trimmed] = turn_directions(exptData,  exptMeta, fwdvelocity_segments_optoon_trimmed, rotvel_segments_optoon_trimmed);

%% find turns
yaw_information_right = findYawVelPeaksFT(exptData, 100, [0.2,0.7], exptData.motion, 1, 0);
yaw_information_left = findYawVelPeaksFT(exptData, 100, [0.2,0.7], exptData.motion, 0, 0);

%% now here i will have a function to extract saccades
[exptData, saccade_both] = find_saccades(exptData, yaw_information_right,yaw_information_left, forvel_cutoff);

%% flat path
%plot_flatpath_saccades(exptData, 1, savepath)

%% then maybe plot rasters here
%plot_raster(exptData,saccade_both, savepath)

%% plot flat path opto vs. saccades


