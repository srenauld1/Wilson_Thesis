% Simple, direct approach:
function data_out = remove_timepoints(data_in, keep_idx)
    % Initialize output
    data_out = data_in;
    
    % Process top-level fields
    if isfield(data_in, 't')
        data_out.t = data_in.t(keep_idx);
    end
    if isfield(data_in, 'epochinds')
        data_out.epochinds = data_in.epochinds(keep_idx);
    end
    if isfield(data_in, 'dff')
        data_out.dff = data_in.dff(keep_idx);
    end
    
    % Handle nested structs
    if isfield(data_in, 'flypos')
        fields = fieldnames(data_in.flypos);
        for i = 1:length(fields)
            if size(data_in.flypos.(fields{i}), 2) == length(data_in.t)
                data_out.flypos.(fields{i}) = data_in.flypos.(fields{i})(:, keep_idx);
            end
        end
    end
    
    if isfield(data_in, 'ball')
        fields = fieldnames(data_in.ball);
        for i = 1:length(fields)
            if size(data_in.ball.(fields{i}), 2) == length(data_in.t)
                data_out.ball.(fields{i}) = data_in.ball.(fields{i})(:, keep_idx);
            end
        end
    end
    
    if isfield(data_in, 'vis')
        fields = fieldnames(data_in.vis);
        for i = 1:length(fields)
            if size(data_in.vis.(fields{i}), 2) == length(data_in.t)
                data_out.vis.(fields{i}) = data_in.vis.(fields{i})(:, keep_idx);
            end
        end
    end
end