function fly_collection = load_fly_data_structs(parent_folder, subfolders)
% load_fly_data_structs loads fly_data structs from smr_analysis/fly_data.mat in specified subfolders

if nargin < 2 || isempty(subfolders)
    % Scan all direct subfolders (not recursive)
    d = dir(parent_folder);
    is_dir = [d.isdir] & ~ismember({d.name},{'.','..'});
    subfolders = {d(is_dir).name};
end

fly_collection = struct();

for k = 1:numel(subfolders)
    this_analysis_folder = fullfile(parent_folder, subfolders{k}, 'smr_analysis');
    if ~isfolder(this_analysis_folder)
        warning('No smr_analysis folder in %s, skipping.', fullfile(parent_folder, subfolders{k}));
        continue
    end
    all_files = dir(fullfile(this_analysis_folder, '*.mat'));
    for i = 1:numel(all_files)
        mat_path = fullfile(all_files(i).folder, all_files(i).name);
        S = load(mat_path);

        % Try to identify fly_data struct automatically
        if isfield(S, 'fly_data')
            fly_data = S.fly_data;
        else
            % Try to find variable that matches your structure
            varnames = fieldnames(S);
            found = false;
            for j = 1:length(varnames)
                candidate = S.(varnames{j});
                if isstruct(candidate) && isfield(candidate,'daq')
                    fly_data = candidate;
                    found = true;
                    break;
                end
            end
            if ~found
                warning('No fly_data struct found in %s', mat_path);
                continue
            end
        end

        % --- fly_id or fallback ---
        if isfield(fly_data, 'fly_id')
            fly_id = fly_data.fly_id;
        elseif isfield(fly_data, 'id') && isfield(fly_data.id, 'recid')
            fly_id = fly_data.id.recid;
        elseif isfield(fly_data, 'daq') && isfield(fly_data.daq, 'recid')
            fly_id = fly_data.daq.recid;
        else
            warning('No fly_id or id.recid found—skipping file %s', mat_path);
            continue
        end
        fly_id_field = matlab.lang.makeValidName(fly_id);

        % If duplicate, warn & skip (or append a suffix if you prefer)
        if isfield(fly_collection, fly_id_field)
            warning('Duplicate fly_id: %s; skipping this entry (%s)', fly_id, mat_path)
            continue
        end

        fly_collection.(fly_id_field) = fly_data;
        fprintf('Loaded %s --> %s\n', mat_path, fly_id_field);
    end
end
end