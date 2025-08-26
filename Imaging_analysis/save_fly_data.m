function save_fly_data(daq, ts, dff_motion, fly_identifier, savepath, varargin)
% Save fly data with identifier for later batch analysis
%
% Inputs:
%   daq - daq structure with kinematic data
%   ts - time series data (dff)
%   fly_identifier - unique identifier for this fly (string or char)
%   savepath - directory to save data
%   varargin - optional additional data to save (e.g., analysis results)
%
% Usage:
%   save_fly_data(daq, ts, 'fly_001', savepath)
%   save_fly_data(daq, ts, 'fly_001', savepath, 'analysis_results', results)

    % Create filename with fly identifier
    filename = sprintf('fly_data_%s.mat', fly_identifier);
    full_filepath = fullfile(savepath, filename);
    
    % Prepare data structure
    fly_data = struct();
    fly_data.daq = daq;
    fly_data.ts = ts;
    fly_data.dff_motion = dff_motion;
    fly_data.fly_id = fly_identifier;
    fly_data.save_date = datestr(now);
    
    % Add any additional data passed as arguments
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
            if i+1 <= length(varargin)
                field_name = varargin{i};
                field_data = varargin{i+1};
                fly_data.(field_name) = field_data;
            end
        end
    end
    
    % Save the data
    save(full_filepath, 'fly_data');
    fprintf('Fly data saved: %s\n', full_filepath);
end