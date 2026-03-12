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

% SMR updated to plot the unsmoothed data only 3//11/2026

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
% does not smooth data
[exptData, exptMeta] = process_fictrac_panels(exptData,exptMeta, minVel);

%% plot velocity data split up by opto chunks
% using raw data alone
[exptData,  exptMeta, fwdvelocity_segments_optoon_trimmed, rotvel_segments_optoon_trimmed] = plotExpt_sliced_opto_velocity(exptData, exptMeta, savepath);

%% plot flatpath with opto on and off
plot_flatpath_opto(exptData, savepath)

% <<<<<<< HEAD
%% now plot sinuosity
% S = benhamou_sinuosity(exptData.x', exptData.y');




% %% now plot sinuosity
% S = benhamou_sinuosity(exptData.x', exptData.y');
% 
% 
% %% calculate sinuosity
% 
% optostim = exptData.optoStim(:);
% fs = 1/median(diff(exptData.t)); % e.g. 100 Hz
% window_idx = round(1*fs);
% x = exptData.x;  % [Nx1] or [1xN]
% y = exptData.y;
% 
% % --- Find ON windows ---
% on_starts = find([false; diff(optostim)==0.25]);
% on_ends   = find([false; diff(optostim)==-0.25]) - 1;
% if optostim(1)==1, on_starts = [1; on_starts]; end
% if optostim(end)==1, on_ends = [on_ends; length(optostim)]; end
% 
% % --- Find OFF windows ---
% off_starts = find([false; diff(optostim)==-0.25]);
% off_ends   = find([false; diff(optostim)==0.25]) - 1;
% if optostim(1)==0, off_starts = [1; off_starts]; end
% if optostim(end)==0, off_ends = [off_ends; length(optostim)]; end
% 
% % --- 1. Optostim ON windows: take entire periods ---
% sinuosity_on = nan(length(on_starts),1);
% for i = 1:length(on_starts)
%     idxs = on_starts(i):on_ends(i);
%     if length(idxs)<2, continue; end
%     S = benhamou_sinuosity(x(idxs)', y(idxs)');
%     sinuosity_on(i) = S;
% end
% 
% % --- 2. Optostim OFF: pick a random 1s window (between 4–10s in) ---
% rng(2); % fixed seed for reproducibility
% sinuosity_off = nan(length(off_starts),1);
% for i = 1:length(off_starts)
%     idx_start = off_starts(i);
%     idx_end   = off_ends(i);
%     win_len = idx_end - idx_start + 1;
%     min_offset = round(4*fs); 
%     max_offset = round(10*fs); 
%     window_idx = round(1*fs);
% 
%     first_possible = idx_start + min_offset;
%     last_possible  = idx_start + max_offset - window_idx;
%     if last_possible > idx_end - window_idx + 1
%         last_possible = idx_end - window_idx + 1;
%     end
%     possible_starts = first_possible:last_possible;
%     if isempty(possible_starts) || possible_starts(1)+window_idx-1 > idx_end
%         continue
%     end
%     idx0 = possible_starts(randi(length(possible_starts)));
%     idxs = idx0:(idx0+window_idx-1);
%     S = benhamou_sinuosity(x(idxs)', y(idxs)');
%     sinuosity_off(i) = S;    
% end
% 
% % Remove NaNs for averaging and plotting
% sinuosity_on  = sinuosity_on(~isnan(sinuosity_on));
% sinuosity_off = sinuosity_off(~isnan(sinuosity_off));
% 
% % --- Plot or compare
% figure
% boxplot([sinuosity_on; sinuosity_off], ...
%     [repmat({'Opto ON'},length(sinuosity_on),1); ...
%      repmat({'Opto OFF'},length(sinuosity_off),1)]);
% ylabel('Sinuosity (Benhamou)');
% title('Optostim Entire Window (ON) vs. 1s Random in [4,10]s during OFF');
% 
% fprintf('Mean ON: %.3f, Mean OFF: %.3f\n', mean(sinuosity_on), mean(sinuosity_off));
% save_plot_with_title_as_filename("sinuosity_opto_on", "opto_off", savepath)

