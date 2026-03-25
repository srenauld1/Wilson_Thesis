%% This script is an amalgamation of all of the behavior data processing scripts:

% This first function re-orients the G4 panel
% s to the true front (rather
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

% SMR updated to plot the unsmoothed data only 3/11/2026
% added aligned flat path plot 3/25/2026


close all
clear all

%% define initial variables
minVel = 3; % minimum total velocity to be count as moving

%% load in exptData and exptMeta
% Define the base folder path
baseFolder = '/Volumes/Neurobio/wilsonlab/Sophia/behavior/'; % Change this to the desired directory

% Prompt the user to select a file from the base folder
[fileName, pathName] = uigetfile('*.mat', 'Select Processed Behavior Data', baseFolder);

if isequal(fileName, 0)
    disp('User canceled the file selection.');
    return  % ← add this
end
fullFilePath = fullfile(pathName, fileName);
load(fullFilePath);
disp(['User selected "', fullFilePath, '".']);
if ~exist('exptData','var') || ~exist('exptMeta','var')
    error('Loaded file does not contain exptData and/or exptMeta: %s', fullFilePath);
end

%% create savepath and folder directory
% Trim the file name to get the base path.
% We do this by removing everything after the last forward slash.

[base_path, fname, ~] = fileparts(fullFilePath);  % already defined above
fname_parts = strsplit(fname, ' ');
% fname = '2026_03_04_fly01_cell01_trial01A pattern20 stim2 function1_pro'
% fname_parts{1} = '2026_03_04_fly01_cell01_trial01A'
% fname_parts{2} = 'pattern20'
% fname_parts{3} = 'stim2'
% fname_parts{4} = 'function1_pro'

if numel(fname_parts) < 2
    error('Filename format not recognized (expected spaces between trial/stim info): %s', fname);
end

stim_part   = fname_parts{end-1};           % 'stim2'
stim_number = extractAfter(stim_part, 'stim');   % '2'
trial_part  = fname_parts{1};               % '2026_03_04_fly01_cell01_trial01A'
trial_num   = extractAfter(trial_part, 'trial'); % '01A'

if isempty(stim_number) || isempty(trial_num)
    error('Could not parse stim number or trial number from filename: %s', fname);
end
if isempty(stim_number) || isempty(trial_num)
    error('Could not parse stim number or trial number from filename: %s', fullFilePath);
end
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
% does not smooth data
[exptData, exptMeta] = process_fictrac_panels(exptData,exptMeta, minVel);

%% plot velocity data split up by opto chunks
% using raw data alone
[exptData, exptMeta, fwdvelocity_segments_optoon_extra, rotvel_segments_optoon_extra] = ...
    plotExpt_sliced_opto_velocity(exptData, exptMeta, savepath, minVel);
%% plot flatpath with opto on and off
plot_flatpath_opto(exptData, savepath)

%% now plot flatpath collapsed over whole trial
plot_flatpath_aligned_opto(exptData, savepath)

%% save expt data
%% save processed exptData and exptMeta
save_filename = fullfile(savepath, [fname, '_processed.mat']);

try
    save(save_filename, 'exptData', 'exptMeta');
    disp(['Saved processed data to: ', save_filename]);
catch ME
    warning('Could not save processed data to %s.\nError: %s', ...
            save_filename, ME.message);
end

%% end!

