function batch_data = load_flies_manual(base_directory, varargin)
% Load fly_data files with manual file selection
%
% Inputs:
%   base_directory - starting directory for file selection dialog
%   varargin - optional parameters:
%     'verbose' - display loading progress (default: true)
%
% Outputs:
%   batch_data - structure containing all loaded fly data
%
% Usage:
%   batch_data = load_flies_manual('/Users/sophiarenauld/stacks/')

    % Parse inputs
    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    parse(p, varargin{:});
    
    verbose = p.Results.verbose;
    
    % Initialize storage
    fly_counter = 0;
    current_path = base_directory;
    
    % Initialize output structure
    batch_data = struct();
    batch_data.flies = {};
    batch_data.fly_ids = {};
    
    % File selection loop
    fprintf('Select fly data files one by one. Click Cancel when finished.\n');
    
    while true
        % Prompt user to select a fly data file
        prompt_text = sprintf('Select Fly Data File #%d (Cancel to finish)', fly_counter + 1);
        [fileName, pathName] = uigetfile('*fly_data*.mat', prompt_text, current_path);
        
        % Check if user canceled
        if isequal(fileName, 0)
            if fly_counter == 0
                fprintf('No files selected. Exiting.\n');
                batch_data = [];
                return;
            else
                fprintf('File selection complete. Selected %d files.\n', fly_counter);
                break;
            end
        end
        
        fly_counter = fly_counter + 1;
        current_path = pathName; % Remember the path for next selection
        
        % Load the selected file
        filepath = fullfile(pathName, fileName);
        
        if verbose
            fprintf('Loading file %d: %s\n', fly_counter, fileName);
        end
        
        try
            % Load the fly_data structure
            loaded = load(filepath);
            
            % Store the fly data
            batch_data.flies{fly_counter} = loaded.fly_data;
            batch_data.fly_ids{fly_counter} = loaded.fly_data.fly_id;
            
            if verbose
                fprintf('  Loaded fly: %s\n', loaded.fly_data.fly_id);
            end
            
        catch ME
            warning('Failed to load file %s: %s', fileName, ME.message);
            batch_data.flies{fly_counter} = [];
            batch_data.fly_ids{fly_counter} = sprintf('ERROR_%d', fly_counter);
        end
    end
    
    % Remove any failed loads
    valid_flies = ~cellfun(@isempty, batch_data.flies);
    batch_data.flies = batch_data.flies(valid_flies);
    batch_data.fly_ids = batch_data.fly_ids(valid_flies);
    
    % Add metadata
    batch_data.metadata.num_flies = length(batch_data.flies);
    batch_data.metadata.num_valid_flies = sum(valid_flies);
    batch_data.metadata.base_directory = base_directory;
    batch_data.metadata.load_date = datestr(now);
    
    if verbose
        fprintf('Successfully loaded %d flies\n', sum(valid_flies));
        fprintf('Fly IDs: %s\n', strjoin(batch_data.fly_ids, ', '));
    end
    
    % Add data access functions
    batch_data = add_data_access_functions(batch_data);
end

function batch_data = add_data_access_functions(batch_data)
% Add helper functions for accessing the loaded data

    % Function to get all dff data
    batch_data.get_all_dff = @() cellfun(@(x) x.ts{1}, batch_data.flies, 'UniformOutput', false);
    
    % Function to get specific kinematic variable from all flies
    batch_data.get_kinematics = @(var_name) cellfun(@(x) x.daq.(var_name), batch_data.flies, 'UniformOutput', false);
    
    % Function to find fly by ID
    batch_data.find_fly = @(fly_id) batch_data.flies{strcmp(batch_data.fly_ids, fly_id)};
    
    % Function to get all daq structures
    batch_data.get_all_daq = @() cellfun(@(x) x.daq, batch_data.flies, 'UniformOutput', false);
    
    % Function to get all ts data
    batch_data.get_all_ts = @() cellfun(@(x) x.ts, batch_data.flies, 'UniformOutput', false);
    
    if batch_data.metadata.num_valid_flies > 0
        fprintf('\nData access functions added:\n');
        fprintf('  batch_data.get_all_dff() - get all dff data\n');
        fprintf('  batch_data.get_kinematics(''var_name'') - get kinematic variable from all flies\n');
        fprintf('  batch_data.find_fly(''fly_id'') - find specific fly data\n');
        fprintf('  batch_data.get_all_daq() - get all daq structures\n');
        fprintf('  batch_data.get_all_ts() - get all ts data\n');
    end
end