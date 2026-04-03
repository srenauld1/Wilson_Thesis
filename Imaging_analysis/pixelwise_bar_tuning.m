function [pref_bin_vol, pref_amp_vol] = pixelwise_bar_tuning( ...
            stack, bar_phase, nBins, savepath, phase_thresh, tuning_thresh)
% PIXELWISE_BAR_TUNING_DFF
%   Compute voxel-wise preferred bar position using dF/F (baseline = 5th
%   percentile of each voxel's fluorescence over time).
%
%   Bar phase is continuous in [-3,3], where -3 and +3 both represent
%   0/360 degrees. Only frames where the bar is MOVING are used.
%
% Inputs:
%   stack         : [nx x ny x nz x nt] imaging data (raw F)
%   bar_phase     : [nt x 1] or [1 x nt] bar phase in [-3, 3]
%   nBins         : number of spatial bins (e.g. 16, 24)
%   savepath      : directory to save tuning maps figure (optional, [] to skip)
%   phase_thresh  : minimum |Δphase| between frames for "moving" (default 1e-4)
%   tuning_thresh : minimum tuning range (max-min across bins) in dF/F
%                   required to consider voxel as tuned (default 0.05)
%
% Outputs:
%   pref_bin_vol : [nx x ny x nz] preferred bin index (1..nBins, NaN if untuned)
%   pref_amp_vol : [nx x ny x nz] mean dF/F at that preferred bin

if nargin < 6 || isempty(tuning_thresh)
    tuning_thresh = 0.05;   % default dF/F range
end
if nargin < 5 || isempty(phase_thresh)
    phase_thresh = 1e-4;
end
if nargin < 4
    savepath = [];
end

% ------------- Basic checks -------------
sz = size(stack);
if numel(sz) ~= 4
    error('stack must be 4D: [x y z t].');
end
[nx, ny, nz, nt] = deal(sz(1), sz(2), sz(3), sz(4));

bar_phase = bar_phase(:);
if numel(bar_phase) ~= nt
    error('Length of bar_phase (%d) must match # of frames in stack (%d).', ...
          numel(bar_phase), nt);
end

% ------------- Compute dF/F per voxel -------------
stack_d = double(stack);            % [nx x ny x nz x nt]

% F0 = 5th percentile over time for each voxel
F0 = prctile(stack_d, 5, 4);        % [nx x ny x nz]
F0(F0 == 0 | isnan(F0)) = eps;     % avoid divide-by-zero

% dFF: [nx x ny x nz x nt]
dff = bsxfun(@rdivide, bsxfun(@minus, stack_d, F0), F0);

% ------------- Map phase [-3,3] -> angle [0,360) -------------
angle_deg = ((bar_phase + 3) / 6) * 360;
angle_deg = mod(angle_deg, 360);   % [0,360)

% ------------- Determine when bar is moving -------------
dPhase = [0; diff(bar_phase)];     % length nt
moving = abs(dPhase) > phase_thresh;

% ------------- Bin angles into nBins -------------
edges  = linspace(0, 360, nBins+1);
binIdx = discretize(angle_deg, edges);   % 1..nBins or NaN
binIdx(~moving) = NaN;                   % ignore static frames

% ------------- Reshape dF/F: [nVox x nt] -------------
V = reshape(dff, [], nt);                % [nVox x nt]
nVox = size(V,1);

% ------------- Compute tuning: mean dF/F per voxel per bin -------------
tuning = nan(nVox, nBins);               % [nVox x nBins]

for b = 1:nBins
    idx = (binIdx == b);
    if ~any(idx)
        continue;
    end
    tuning(:,b) = mean(V(:,idx), 2, 'omitnan');   % mean dF/F in bin b
end

% ------------- Tuning strength & significance mask -------------
tuning_range = max(tuning, [], 2, 'omitnan') - min(tuning, [], 2, 'omitnan');
tuned_mask   = tuning_range >= tuning_thresh;

% ------------- Preferred bin and amplitude per voxel -------------
[pref_amp, pref_bin] = max(tuning, [], 2);    % over bins

% voxels with all-NaN tuning or below threshold → NaN
allNaN = all(isnan(tuning), 2);
pref_amp(~tuned_mask | allNaN) = NaN;
pref_bin(~tuned_mask | allNaN) = NaN;

% Reshape volumes
pref_bin_vol = reshape(pref_bin, [nx ny nz]);
pref_amp_vol = reshape(pref_amp, [nx ny nz]);

% ------------- Plot preferred bar bin per z-plane -------------
figure('Name','Pixelwise preferred bar position (dF/F, tuned only)','Color','w');
nCols = ceil(sqrt(nz));
nRows = ceil(nz / nCols);

cmap = hsv(nBins);  % one color per bin

for z = 1:nz
    subplot(nRows, nCols, z);

    plane_bins = pref_bin_vol(:,:,z);   % 1..nBins or NaN

    % Background black for untuned voxels
    set(gca,'Color','k');   % axes background is black

    % AlphaData: only plot where plane_bins is finite (tuned)
    alphaMask = ~isnan(plane_bins);

    im = imagesc(plane_bins, [1 nBins]);   % show bins 1..nBins
    im.AlphaData = alphaMask;              % NaN => alpha=0 => background
    axis image off;
    colormap(gca, cmap);

    c = colorbar;
    c.Ticks = 1:nBins;
    c.TickLabels = arrayfun(@(k) sprintf('%d',k), 1:nBins, 'UniformOutput', false);

    title(sprintf('z = %d (preferred bar bin)', z));
end

sgtitle(sprintf('Pixelwise preferred bar position (moving; dF/F, range >= %.3f)', ...
                tuning_thresh));
if ~isempty(savepath)
    save_plot_with_title_as_filename('bar_tuning_dff_tuned', 'volume', savepath);
end
end