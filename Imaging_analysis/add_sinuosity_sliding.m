function a2p_data = add_sinuosity_sliding(a2p_data, window)
% Add pointwise sinuosity to a2p_data using a sliding window.
% Input:
%   a2p_data - struct with fields x and y (both nx1 or 1xn)
%   window   - number of POINTS in the window (should be odd for symmetric)
%
% Output:
%   a2p_data.sinuosity - [n x 1] vector with sliding window sinuosity

    x = a2p_data.dq(2).pxb(:); % or a2p_data.x(:)
    y = a2p_data.dq(2).pyb(:); % or a2p_data.y(:)
    N = numel(x);

    halfwin = floor(window / 2);
    sinuosity = nan(N, 1);

    for i = 1:N
        w_start = max(1, i - halfwin);
        w_end   = min(N, i + halfwin);

        xw = x(w_start:w_end);
        yw = y(w_start:w_end);

        % Must be vector; re-shape as nx1 for "benhamou_sinuosity"
        % This uses benhamou_sinuosity with one column as trajectory.
        try
            s_val = benhamou_sinuosity(xw, yw);
            sinuosity(i) = s_val;
        catch
            sinuosity(i) = NaN;
        end
    end

    a2p_data.sinuosity = sinuosity;
end