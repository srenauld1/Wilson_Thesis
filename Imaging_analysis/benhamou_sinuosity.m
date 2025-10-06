function S = benhamou_sinuosity(x, y, quarters)
% Compute ecological sinuosity (Benhamou, 2004) for 2D paths
% Input:  x, y = matrices of coordinates (rows = time steps, cols = trajectories)
%         quarters (optional logical): if true, compute sinuosity in 4 segments per trajectory
% Output: S = row vector (#traj) if quarters not specified, or (#traj x 4) if quarters==true

    % Check quarters argument
    if nargin < 3
        quarters = false;
    end

    % Ensure inputs are the same size
    if ~isequal(size(x), size(y))
        error('x and y must have the same size');
    end

    nTraj = size(x,2);       % number of trajectories

    if ~quarters
        S = nan(1, nTraj);   % preallocate output
    else
        S = nan(nTraj, 4);   % preallocate for 4 quarters
    end
    
    for k = 1:nTraj
        % Extract trajectory
        xk_full = x(:,k);
        yk_full = y(:,k);

        % Remove NaNs (if any)
        valid = ~(isnan(xk_full) | isnan(yk_full));
        xk_full = xk_full(valid);
        yk_full = yk_full(valid);

        nPts = numel(xk_full);

        if nPts < 3
            if ~quarters
                S(k) = NaN;
            else
                S(k,:) = NaN;
            end
            continue;
        end

        if ~quarters
            % --------- Whole track ------------
            S(k) = benhamou_core(xk_full, yk_full);
        else
            % -------- Quartered --------------
            qedges = round(linspace(1, nPts+1, 5)); % define four intervals
            for q = 1:4
                idx = qedges(q):(qedges(q+1)-1);
                if numel(idx) < 3
                    S(k,q) = NaN;
                else
                    S(k,q) = benhamou_core(xk_full(idx), yk_full(idx));
                end
            end
        end
    end
end

% --------- internal subfunction for single segment sinuosity-------------
function Sval = benhamou_core(xk, yk)
    dx = diff(xk);
    dy = diff(yk);
    steps = sqrt(dx.^2 + dy.^2);
    p = mean(steps);

    b = std(steps) / p;

    ang = atan2(dy, dx);
    dtheta = diff(ang);
    dtheta = mod(dtheta + pi, 2*pi) - pi;
    c = mean(cos(dtheta));

    denominator = p * ((1 + c)/(1 - c) + b^2);
    Sval = 2 / sqrt(denominator);
end