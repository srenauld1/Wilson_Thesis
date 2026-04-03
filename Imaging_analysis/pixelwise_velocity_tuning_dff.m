function [pref_state_vol, pref_amp_vol] = pixelwise_velocity_tuning_dff( ...
            stack, vel, vel_thresh, diff_thresh, savepath, smooth_sigma)
% PIXELWISE_VELOCITY_TUNING_DFF
%
%   Compute voxel-wise preferred velocity state from volumetric imaging and
%   a 1D velocity time series, using dF/F (baseline = 5th percentile).
%
%   Velocity states:
%       1 = still           (abs(vel) <= vel_thresh)
%       2 = moving positive ( vel  > vel_thresh)
%       3 = moving negative ( vel  < -vel_thresh)
%
%   For each voxel:
%       - Compute mean dF/F in each of the 3 states.
%       - Pick the state with highest mean as "preferred".
%       - Keep only voxels where (max - min) across states >= diff_thresh.
%       - Also compute a contrast measure: (p - n) / (p + n), capped at ±0.5.
%
%   Optionally, the contrast volume is spatially smoothed (3-D Gaussian)
%   before plotting.
%
% Inputs:
%   stack       : [nx x ny x nz x nt] imaging data (raw fluorescence)
%   vel         : [nt x 1] or [1 x nt] velocity trace
%   vel_thresh  : scalar, threshold for still vs moving
%   diff_thresh : scalar, min (max-min) dF/F range across states
%                 to be considered tuned
%   savepath    : directory to save the tuning maps figure (optional, [] to skip)
%   smooth_sigma: (optional) std dev of 3-D Gaussian kernel in voxels for
%                 smoothing the contrast volume. 0 or [] = no smoothing.
%
% Outputs:
%   pref_state_vol : [nx x ny x nz] preferred state index
%                    (1=still,2=pos,3=neg, NaN=untuned)
%   pref_amp_vol   : [nx x ny x nz] mean dF/F in that preferred state

if nargin < 5
    savepath = [];
end
if nargin < 6 || isempty(smooth_sigma)
    smooth_sigma = 0;   % no spatial smoothing by default
end

% ------------ Basic checks ------------
sz = size(stack);
if numel(sz) ~= 4
    error('stack must be 4D: [x y z t].');
end
[nx, ny, nz, nt] = deal(sz(1), sz(2), sz(3), sz(4));

vel = vel(:);
if numel(vel) ~= nt
    error('Length of vel (%d) must match # of frames in stack (%d).', ...
          numel(vel), nt);
end

% ------------ Compute dF/F per voxel (F0 = 5th percentile) ------------
stack_d = double(stack);              % [nx x ny x nz x nt]
F0 = prctile(stack_d, 5, 4);          % [nx x ny x nz]
F0(F0 == 0 | isnan(F0)) = eps;        % avoid /0 and NaNs

% dF/F = (F - F0) / F0
dff = bsxfun(@rdivide, bsxfun(@minus, stack_d, F0), F0);  % [nx x ny x nz x nt]

% ------------ Define velocity states ------------
state_still = abs(vel) <= vel_thresh;
state_pos   = vel >  vel_thresh;
state_neg   = vel < -vel_thresh;

% Ensure mutual exclusivity (if boundaries overlap due to ==, still wins)
state_pos(state_still) = false;
state_neg(state_still) = false;

% ------------ Reshape dF/F to [nVox x nt] ------------
V = reshape(dff, [], nt);            % [nVox x nt]
nVox = size(V,1);

% ------------ Mean dF/F per voxel for each state ------------
nStates = 3;
tuning = nan(nVox, nStates);        % col 1: still, 2: pos, 3: neg

% State 1: still
idx = state_still;
if any(idx)
    tuning(:,1) = mean(V(:,idx), 2, 'omitnan');
end

% State 2: positive
idx = state_pos;
if any(idx)
    tuning(:,2) = mean(V(:,idx), 2, 'omitnan');
end

% State 3: negative
idx = state_neg;
if any(idx)
    tuning(:,3) = mean(V(:,idx), 2, 'omitnan');
end

% ------------ Preferred state and amplitude ------------
[pref_amp, pref_state] = max(tuning, [], 2);  % over columns = states

% Tuning strength across states
tuning_range = max(tuning, [], 2, 'omitnan') - ...
               min(tuning, [], 2, 'omitnan');
allNaN = all(isnan(tuning), 2);

% Mask: must have data & range >= diff_thresh
tuned_mask = ~allNaN & (tuning_range >= diff_thresh);

pref_amp(~tuned_mask)   = NaN;
pref_state(~tuned_mask) = NaN;

% Reshape volumes
pref_state_vol = reshape(pref_state, [nx ny nz]);
pref_amp_vol   = reshape(pref_amp,   [nx ny nz]);

% ------------ Compute (p - n) / (p + n) per voxel ------------
p = tuning(:,2);      % positive state mean
n = tuning(:,3);      % negative state mean

denom = p + n;
num   = p - n;

contrast = num ./ denom;          % (p - n) / (p + n)
contrast(denom == 0) = NaN;       % avoid division by zero
bad = isnan(p) | isnan(n);
contrast(bad) = NaN;

% ---- Cap at [-0.5, 0.5] ----
contrast = max(min(contrast, 0.5), -0.5);

% Reshape to volume
contrast_vol = reshape(contrast, [nx ny nz]);

% ------------ Optional spatial smoothing of contrast ------------
if smooth_sigma > 0
    % If you have Image Processing Toolbox:
    try
        contrast_vol = imgaussfilt3(contrast_vol, smooth_sigma);
    catch
        % Fallback: slice-by-slice 2-D smoothing (no smoothing across z)
        warning('imgaussfilt3 not found; using 2-D imgaussfilt per slice.');
        for z = 1:nz
            contrast_vol(:,:,z) = imgaussfilt(contrast_vol(:,:,z), smooth_sigma);
        end
    end
end

% ------------ Plot: contrast per z-plane ------------
figure('Name','(pos - neg) / (pos + neg) contrast','Color','w');
nCols = ceil(sqrt(nz));
nRows = ceil(nz / nCols);

for z = 1:nz
    subplot(nRows, nCols, z);
    imagesc(contrast_vol(:,:,z), [-0.2 0.2]);   % fixed to [-0.5, 0.5]
    axis image off;
    colormap(gca, jet);                         % or another diverging colormap
    colorbar;
    title(sprintf('z = %d', z));
end

sgtitle(sprintf('(pos - neg) / (pos + neg) dF/F contrast (capped at ±0.25, smoothing σ=%.2f)', ...
        smooth_sigma));

% ------------ Save figure if requested ------------
if ~isempty(savepath)
    save_plot_with_title_as_filename('vel_tuning_dff', 'volume', savepath);
end
end