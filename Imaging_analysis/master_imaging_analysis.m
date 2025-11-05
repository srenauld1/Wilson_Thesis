%% This script is an amalgamation of all of the calcium data processing scripts:

% This first function re-orients the G4 panels to the true front (rather
% than where the panels decide is 0) 0 is in front of the fly
% This function aclolso downsamples time, opto, and the panels such that they
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
an4_project = 0;
jump = 0;
lplc4=0;
visual = 0;
box=0;
split = 0;
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
    daq = load(fullFilePath);
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
    an4_project = 1;  
elseif contains(base_path, 'LPLC4')
    an4_project = 1;  
    visual = 1;
elseif contains(base_path, 'LC22')
    an4_project = 1;  
    visual = 1;
elseif contains(base_path, 'LT51')
    an4_project = 1;  
    visual = 1;
elseif contains(base_path, 'TMY20')
    an4_project = 1;  
    visual = 1;
elseif contains(base_path, 'tmy20')
    an4_project = 1;  
    visual = 1;
elseif contains(base_path, 'cmllp01')
    an4_project = 1; 
elseif contains(base_path, 'cmllpo1')
    an4_project = 1;  
elseif contains(base_path, 'LC33')
    visual = 1;
else
    an4_project = 0;  
end

if contains(base_path, 'box')
    box=1;
    visual = 1; 
elseif contains(base_path, 'olgrating')
    visual = 1; 
elseif contains(base_path, 'vgrating')
    visual = 1; 
    split=0;
elseif contains(base_path, 'hgrating')
    visual = 1; 
elseif contains(base_path, 'splitgrating')
    visual = 1;
    split=1;
% elseif contains(base_path, 'bbar')
%     visual = 1;
elseif contains(base_path, 'dbar')
    visual = 1;
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

%% adjust variables
fields = fieldnames(daq);
patterns = {'byv','bfv','bsv','byv_supp','bfv_supp','bsv_supp', 'vy', 'vy_supp', 'vyv', 'vyv_supp'};
replacements = {'bvy','bvf','bvs','bvy_supp','bvf_supp','bvs_supp', 'vh', 'vh_supp', 'vvy', 'vvy_supp'};

for i = 1:numel(fields)
    oldname = fields{i};
    newname = oldname;
    for p = 1:numel(patterns)
        % Only do replacement if matches exactly (for _supp and such)
        if strcmp(oldname, patterns{p})
            newname = replacements{p};
            break % stop after first exact match
        end
    end
    if ~strcmp(oldname, newname)
        daq.(newname) = daq.(oldname);
        daq = rmfield(daq, oldname);
    end
end

if exist('ts','var')
    if iscell(ts) && numel(ts)==1
        ts = ts{1};
    end
end

if exist('ts','var')
    dat{1,1}.ts = ts;
end

% Assume dat is a 1x1 cell containing a struct
if iscell(dat) && numel(dat)==1 && isstruct(dat{1})
    dat = dat{1};
end

%% process and pick out jumps
vel_on = 0.5;  % Speed above which fly is moving (enter movement)
vel_off = 0.5; % Speed below which fly is NOT moving (exit movement)
forvel_cutoff = 2;
[daq, yaw_information_right,yaw_information_left,yaw_information_left_supp, yaw_information_right_supp ] = process_fictrac_panels_2p(daq, vel_on, vel_off, jump);
%% now here i will have a function to extract saccades
[daq, saccade_both] = find_saccades_2p(daq, yaw_information_right,yaw_information_left, forvel_cutoff, 0);
[daq, saccade_both_supp] = find_saccades_2p(daq, yaw_information_right_supp,yaw_information_left_supp, forvel_cutoff, 1);
%% plot basic variables
if neck
    plotting_neck(daq, dat, savepath)
end

if an4_project % this includes all other visual cells
    plotting_an4(daq, dat, jump, savepath)
    lm_velocity(daq, dat, savepath)
    
end
if visual
    [daq, dat, dff_motion] = openloop_gratings_process(daq, dat, savepath, box, split);
end

if ves041
    %% basic velocity and dff plotting
    plotting_ves041(daq, dat, jump, savepath)
    lm_velocity(daq, dat, savepath)
    save_fly_data(daq, dat, savepath)
else
    save_fly_data(daq, dat, dff_motion, savepath)
end
