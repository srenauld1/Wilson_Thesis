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

%% load in ts
% Define the base folder path
baseFolder = '/Users/sophiarenauld/stacks/'; % Change this to the desired directory

% Prompt the user to select a file from the base folder
[fileName, pathName] = uigetfile('ts.mat', 'Select Processed Imaging and Behavior Data', baseFolder);

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

%% process experimental data - ves041
if ves041
    plotting_ves041(ts, savepath)
end




