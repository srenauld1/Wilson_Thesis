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

%close all
%clear all

%% initialize variables
neck = 0;
ves041 = 0;
an4_project = 0;
jump = 0;
lplc4=0;
visual = 0;
box=0;
split = 0;
all_8=0;
do_pause = 0;
LPLC4=0;
odor = 0;
%% load in daq and ts
% Define the base folder path
baseFolder = '/Users/sophiarenauld/stacks/'; % Change this to the desired directory

% For interactive use, prompt for a processed file as before. For batch use,
% callers may define fullFilePath before running this script.
if ~exist('fullFilePath', 'var') || isempty(fullFilePath)
    [fileName, pathName] = uigetfile('*_s_.mat', 'Select Processed Data', baseFolder);
    if isequal(fileName, 0)
        error('master_imaging_analysis:NoInput', ...
            'User canceled the file selection.');
    end
    fullFilePath = fullfile(pathName, fileName);
end

if ~isfile(fullFilePath)
    error('master_imaging_analysis:InputNotFound', ...
        'Processed input does not exist: %s', fullFilePath);
end

if exist('analysisOmitImageStack','var') && analysisOmitImageStack
    % Behavior/ROI reruns do not need a second copy of the multi-GB movie.
    inputVariables=whos('-file',fullFilePath);
    inputVariables=setdiff({inputVariables.name},{'stack'},'stable');
    a2p_data=load(fullFilePath,inputVariables{:});
else
    a2p_data = load(fullFilePath);
end
a2p_data = import_ves041_protocol(a2p_data,fileparts(fullFilePath));
[a2p_data, compatibility] = harmonize_a2p_data(a2p_data);
if compatibility.synthesizedSupplementalDq
    warning('master_imaging_analysis:SingleDqStream', ...
        ['Only dq(1) is present. dq(2) was copied from dq(1); analyses ' ...
         'will use the imaging-aligned behavioral sampling rate.']);
end
disp(['Analyzing "', fullFilePath, '".']);



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
if isfield(a2p_data,'protocol')
    jump=strcmp(a2p_data.protocol.kind,'bar-jump');
end
if contains(base_path, 'an4')
    an4_project = 1;  
elseif contains(base_path, 'LPLC4')
    LPLC4 = 1;
elseif contains(base_path, 'lplc4')
    LPLC4 = 1;
elseif contains(base_path, 'LC22')
    an4_project = 1;  
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
elseif contains(base_path, 'lc33')
    visual = 1;
else
    an4_project = 0;  
end
if contains(base_path, 'odor')
    odor = 1;  % Set neck to 1 if "neck" is found
    ves041=1;
end

% Jump detection requires visual heading. Keep legacy behavior when vh is
% present, but do not fail newer recordings that omit that channel.
if jump && ~compatibility.hasVisualHeading
    warning('master_imaging_analysis:MissingVisualHeading', ...
        'Disabling jump detection because dq.vh is unavailable.');
    jump = 0;
end


% Define the directory for the stimulus
stim_directory = fullfile(base_path, 'smr_analysis');
if exist('analysisOutputFolder','var') && ~isempty(analysisOutputFolder)
    stim_directory=analysisOutputFolder;
end

% Create the directory if it doesn't exist
if ~exist(stim_directory, 'dir')
    mkdir(stim_directory);
end

% Set the savepath variable
savepath = stim_directory;

% Display the save path
disp(['Save path: ', savepath]);

%% new data separation
daq = a2p_data.dq;


%% process and pick out jumps
vel_on = 0.5;  % Forward Speed above which fly is moving (enter movement)
vel_off = 0.5; % Forward Speed below which fly is NOT moving (exit movement)
forvel_cutoff = 2;
a2p_data = process_fictrac_panels_2p(a2p_data, vel_on, vel_off, jump);
%% now here i will have a function to extract saccades - REMOVED
% [a2p_data, saccade_both] = find_saccades_2p(a2p_data, yaw_information_right,yaw_information_left, forvel_cutoff, 0);
% [a2p_data, saccade_both_supp] = find_saccades_2p(a2p_data, yaw_information_right_supp,yaw_information_left_supp, forvel_cutoff, 1);



%% plot basic variables
if neck
    % NOT ADJUSTED FOR NEW A2P
    plotting_neck(a2p_data, dat, savepath)
end
if an4_project % this includes all other visual cells
    plotting_an4(a2p_data, jump, savepath)  %lm_velocity(a2p_data, dat, savepath)
end
if visual
    plotting_an4(a2p_data, jump, savepath)
    if compatibility.hasVisualHeading
        [a2p_data, dff_motion] = openloop_gratings_process(a2p_data, savepath, box, split, LED);
    else
        warning('master_imaging_analysis:MissingVisualHeading', ...
            'Skipping visual-heading/grating analysis because dq.vh is unavailable.');
        visual = 0;
    end
end
if all_8
      eight_visual_pattern_analysis(a2p_data, savepath, box, split)
end
if do_pause
      a2p_data = process_pause_data(a2p_data, savepath, box, split);
end
if ves041
    %% calculate sinuosity and menotaxis
    % menotaxis code
    % window in seconds to look at heading
    window = 15;
    % minVel, minimum forward velocity to qualify
    minVel = 1;
    % high thresh is rho to count as menotaxis
    highThres = 0.88;
    % low thresh is rho to come out of menotaxis
    lowThres = 0.7;
    [a2p_data, triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(a2p_data, savepath, window, minVel,highThres,lowThres, jump, 0);
    
    if compatibility.hasPosition
        a2p_data = plot_flat_path_colored(a2p_data, 'dff', 'data_only', savepath, odor);
        a2p_data = plot_flat_path_colored(a2p_data, 'fwd', 'data_only', savepath, odor);
    end

    % sinuosity
    % window_sec = 5; % desired window in seconds
    % dt = median(diff(a2p_data.dq(2).t)); % or your timebase
    % window_pts = round(window_sec / dt);
    % a2p_data = add_sinuosity_sliding(a2p_data, window_pts);
    % 
    % a2p_data = plot_flat_path_colored(a2p_data, 'sinuosity', 'data_only', savepath, odor);
    % 
    % % plot dff
    % a2p_data = plot_flat_path_colored(a2p_data, 'dff', 'data_only', savepath, odor);

    %% basic velocity and dff plotting
    a2p_data = plotting_ves041(a2p_data, jump, savepath, odor);
    if isfield(a2p_data,'protocol') && any(strcmp(a2p_data.protocol.kind,{'heat-barrier','heat-stripes'}))
        a2p_data.heat_analysis = analyze_ves041_heat_context(a2p_data,savepath);
    end

    % Event-aligned visual-response plots are generated automatically for
    % discrete bar-jump and moving-grating trials.
    if (jump || contains(lower(base_path),'grating')) && compatibility.hasVisualHeading
        plot_ves041_visual_event_responses(a2p_data,savepath,jump);
    elseif jump || contains(lower(base_path),'grating')
        warning('master_imaging_analysis:MissingVisualHeading', ...
            'Skipping VES041 stimulus-aligned plots because dq.vh is unavailable.');
    end
    % NOT ADJUSTED FOR NEW A2P
    %lm_velocity(a2p_data, dat, savepath)

    %% add in the xcorr stuff for behavior / imaging quantification
    % [results_meno_rho, a2p_data] = moving_window_zero_lag_corr(a2p_data, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'rho', 'dff');
    % [results_meno_fwd, a2p_data] = moving_window_zero_lag_corr(a2p_data, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'fwd', 'dff');
    % [results_meno_sinuosity, a2p_data] = moving_window_zero_lag_corr(a2p_data, [0.25, 0.5, 1, 2, 4, 8], [15, 30, 60, 120], 'sinuosity', 'dff');
    % % plot xcorr for fwd and sinuosity
    % a2p_data = plot_flat_path_colored(a2p_data, 'crosscorr', 'whole_range', savepath, odor);
    % a2p_data = plot_flat_path_colored(a2p_data, 'crosscorr', 'whole_range', savepath, odor);

end
%% odor section
if odor
    % for odor experiments, the odor pulse is shown in the bh and bfv
    % channel. so, px and py are correct, not pxb pyb. If i want rotational
    % speed i should look at vh and vvy
    % segment out odor
    odor_pulses = a2p_data.dq(1).bh;
    if isfield(a2p_data.dq(1),'px'), odor_x = a2p_data.dq(1).px; else, odor_x = []; end
    if isfield(a2p_data.dq(1),'py'), odor_y = a2p_data.dq(1).py; else, odor_y = []; end

    % Define ON/OFF based on odor thresholds
    odor_on  = odor_pulses > -0.5;  % logical vector
    odor_off = odor_pulses < -1;    % logical vector
    % Anything in between (-1 to -0.5) will be treated as OFF here, but you can
    % adjust that if desired.
    
    if ~isempty(odor_x) && ~isempty(odor_y)
    figure; hold on;
    h_on = []; h_off = []; h_start = [];
    
    % Loop once through segments and plot with appropriate color
    for i = 1:length(odor_x)-1
        if odor_on(i)
            h = plot(odor_x(i:i+1), odor_y(i:i+1), '-', 'Color', [1 0 0], 'LineWidth', 1.3); % red
            if isempty(h_on), h_on = h; end
        else
            % treat everything else as odor OFF
            h = plot(odor_x(i:i+1), odor_y(i:i+1), 'k-', 'LineWidth', 1.3); % black
            if isempty(h_off), h_off = h; end
        end
    end
    
    % Start point in green
    h_start = plot(odor_x(1), odor_y(1), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g', ...
                   'DisplayName', 'Start');
    
    xlabel('X Position');
    ylabel('Y Position');
    title('Flat Path: Odor ON (Red), Odor OFF (Black), Start (Green)');
    legend([h_on, h_off, h_start], {'Odor ON', 'Odor OFF', 'Start'}, 'Location', 'best');
    
    axis equal; axis tight; grid on; hold off;
    
    save_plot_with_title_as_filename('flatpath', 'odor', savepath);
    else
        warning('master_imaging_analysis:MissingPosition', ...
            'Skipping odor flat-path plot because dq.px/dq.py are unavailable.');
    end

    %% chop up by odor and plot
    time        = a2p_data.dq(1).t;
    fwd_vel     = a2p_data.dq(1).bvf;
    angular_vel = a2p_data.dq(1).vvy;
    dff         = a2p_data.roi.ts;
    
    odor_pulse_analysis(time, fwd_vel, angular_vel, dff, odor_pulses, savepath);

end
%% LPLC4 specific
if LPLC4
    entry_code = input(['Analyze LPLC4 trial: (0=dark, 1 = replay only, 2 = full closed loop gamut, 3 = bright bar with stars, 4=open loop, ' ...
        '5=all 8 gratings, 6= f2bb2f bar) = ']);
    if entry_code == 0 
        plotting_an4(a2p_data, jump, savepath)
    elseif entry_code ==1
        a2p_data = process_pause_data(a2p_data, savepath, box, split, 1);
    elseif entry_code == 2
        a2p_data = process_pause_data(a2p_data, savepath, box, split, 0);
    elseif entry_code == 3
        a2p_data = process_bar_mp_data(a2p_data, savepath, box, split);
    elseif entry_code==4
        [a2p_data, dff_motion] = openloop_gratings_process(a2p_data, savepath, box, split);
        visual = 1;
    elseif entry_code==5
        a2p_data = eight_visual_pattern_analysis(a2p_data, savepath, box, split);
    elseif entry_code == 6
        a2p_data = openloop_bar_bidirectional_process(a2p_data, savepath, box, split);
    end
end

%% save data

if ~ visual
    save_fly_data(a2p_data, savepath)
else
    save_fly_data(a2p_data, dff_motion, savepath)
end

% Hidden runs must still produce FIG files that open visibly for the user.
if ves041
    repair_saved_figure_visibility(savepath);
end

