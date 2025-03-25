function exptData = compute_absolute_circular_diff(exptData)
    % Computes absolute circular differences in ts.vis.yaw
    %
    % Parameters:
    %   ts: structure containing ts.vis.yaw (heading in radians)
    %
    % Returns:
    %   ts: modified structure with ts.vis.absolute_circular_diff added
    
    % Check if heading field exists
    if ~isfield(exptData, 'g4displayXPos')
        error('Field g4displayXPos not found in data structure.');
    end
    
    % Get heading data
    heading = exptData.g4displayXPos;
    
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
    exptData.absolute_circular_diff = circular_diff;
end