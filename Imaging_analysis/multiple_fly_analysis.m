%% This script is for the purpose of amassing calcium imaging data across flies
% SMR 8/22/2025
parent_folder = '/Users/sophiarenauld/stacks';
pattern = '*LPLC4*splitgrating'; 
dirs = dir(fullfile(parent_folder, pattern));
dirs = dirs([dirs.isdir]);
subfolders = {dirs.name};

fly_collection = load_fly_data_structs(parent_folder, subfolders);

%% Access the data. make structs for forward, pattern, visual etc
all_dff = batch_data.get_all_dff();  % Cell array of dff from all flies
forward_vel = batch_data.get_kinematics('bfv_supp');  % Forward velocity from all flies
all_daq = batch_data.get_all_daq();  % All daq structures

% Find specific fly
% fly_data = batch_data.find_fly('20250821-5_cmllp01_olgrating');
% specific_dff = fly_data.ts{1};
% specific_daq = fly_data.daq;

% Check what you loaded
fprintf('Loaded %d flies: %s\n', batch_data.metadata.num_valid_flies, ...
    strjoin(batch_data.fly_ids, ', '));