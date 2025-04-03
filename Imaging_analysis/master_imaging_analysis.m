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
jump = 0;

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
if contains(base_path, 'jump')
    jump = 1;  % Set neck to 1 if "neck" is found
else
    jump = 0;  % Set neck to 0 if "neck" is not found
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
minVel = 3;
forvel_cutoff = 1;
ts = process_fictrac_panels_2p(ts, minVel, jump);

%% find turns
yaw_information_right = findYawVelPeaksFT(ts, 50, [0.2,1], ts.motion, 1);
yaw_information_left = findYawVelPeaksFT(ts, 50, [0.2,1], ts.motion, 0);

%% now here i will have a function to extract saccades
[ts, saccade_both] = find_saccades_2p(ts, yaw_information_right,yaw_information_left, forvel_cutoff);

%% plot basic variables
if ves041
    %% basic velocity and dff plotting
    plotting_ves041(ts, savepath)
    %% plot flat path with saccades
    plot_flatpath_saccades(ts, jump,savepath);
    %%% plot heatmap of forward acceleration and rs and dff
    % ts.ball.foraccel = diff(ts.ball.forvel);
    % ts.ball.yawvel_adjusted = ts.ball.yawdeg(1:end-1); 
    % ts.ball.sidevel_adjusted = ts.ball.sidevel(1:end-1);
    % ts.ball.dffvel_adjusted = ts.resp.i2{1}(1:end-1); 
    % binnedData = dff_v_velocityHeat(ts.ball.foraccel,abs(ts.ball.yawvel_adjusted),abs(ts.ball.sidevel_adjusted),ts.ball.dffvel_adjusted,1);
    %% plot dff split up by saccades 
    plot_dff_saccades(ts, saccade_both)

end


if neck
    plotting_neck(ts, savepath)
end


%% PILOT ZONE


