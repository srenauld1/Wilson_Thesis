function a2p_data = compute_absolute_circular_diff_2p(a2p_data)
    % Computes absolute circular differences in ts.vis.yaw
    %
    % Parameters:
    %   daq: structure containing daq.vy (heading in radians)
    %
    % Returns:
    %   ts: modified structure with ts.absolute_circular_diff added
    
    % Check if heading field exists
    if ~isfield(a2p_data.dq(1), 'vh')
        error('Field daq.vy not found in data structure.');
    end
    
    % Get heading data
    heading =a2p_data.dq(1).vh;
    
    % Check if row or column vector and handle accordingly
    if isrow(heading)
        % For row vector
        heading_diff = [heading(1), diff(heading)];
    else
        % For column vector
        heading_diff = [heading(1); diff(heading)];
    end
    
    % Compute absolute circular difference
    circular_diff = abs(atan2(...
        sin(heading_diff), ...
        cos(heading_diff)));
    
    % Add results to data structure
    a2p_data.dq(1).absolute_circular_diff = circular_diff;
end