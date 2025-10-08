function save_fly_data(daq, dat, varargin)
% Save fly data with identifier for later batch analysis
%
% Required:
%   daq - kinematic data structure
%   ts - time series data (dff)
% Required, but see below:
%   dff_motion (optional, can be omitted)
%   savepath - dir to save data
%   varargin - more fieldname,data pairs

    % Parse inputs (detect whether dff_motion was provided)
    if nargin < 3
        error('Not enough input arguments.');
    end
    if ischar(varargin{1}) || isstring(varargin{1})
        % dff_motion omitted; arguments are daq, ts, savepath, ...
        dff_motion = [];
        savepath = varargin{1};
        extra_args = varargin(2:end);
    else
        % dff_motion present; arguments are daq, ts, dff_motion, savepath, ...
        dff_motion = varargin{1};
        savepath = varargin{2};
        extra_args = varargin(3:end);
    end

    % Create filename with fly identifier
    filename = 'fly_data.mat';
    full_filepath = fullfile(savepath, filename);

    % Prepare data structure
    fly_data = struct();
    fly_data.daq = daq;
    fly_data.dat = dat;
    if ~isempty(dff_motion)
        fly_data.dff_motion = dff_motion;
    end
    fly_data.save_date = datestr(now);

    % Add any additional data passed as arguments
    if ~isempty(extra_args)
        for i = 1:2:length(extra_args)
            if i+1 <= length(extra_args)
                field_name = extra_args{i};
                field_data = extra_args{i+1};
                fly_data.(field_name) = field_data;
            end
        end
    end

    % Save the data
    save(full_filepath, 'fly_data');
    fprintf('Fly data saved: %s\n', full_filepath);
end