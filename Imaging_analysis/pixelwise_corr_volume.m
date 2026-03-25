function corr_vol = pixelwise_corr_volume(stack, var, savepath, var_name)
% PIXELWISE_CORR_VOLUME
%   corr_vol = pixelwise_corr_volume(stack, var, savepath, var_name)
%
% Inputs:
%   stack    : [nx x ny x nz x nt] imaging data
%   var      : [nt x 1] (or 1 x nt) time series to correlate with (e.g. fwd velocity)
%   savepath : folder where figure will be saved (optional, can be [])
%   var_name : string/name of var for labels (optional)
%
% Output:
%   corr_vol : [nx x ny x nz] correlation coefficient for each voxel

if nargin < 4 || isempty(var_name)
    var_name = 'variable';
end
if nargin < 3
    savepath = [];
end

% ----- Basic checks -----
sz = size(stack);
if numel(sz) ~= 4
    error('stack must be 4D: [x y z t].');
end
nx = sz(1); ny = sz(2); nz = sz(3); nt = sz(4);

var = var(:);  % column vector
if numel(var) ~= nt
    error('Length of var (%d) must match number of frames/volumes in stack (%d).', numel(var), nt);
end

% ----- Reshape and convert to double -----
stack_d = double(stack);                 % avoid integer overflow
V = reshape(stack_d, [], nt);            % [nVoxels x nt], nVoxels = nx*ny*nz

% Optionally handle NaNs: ignore voxels with all NaNs
allNaN = all(isnan(V), 2);
V(allNaN,:) = 0;   % temporary; will set r to NaN later for these voxels

% ----- Demean -----
V = V - mean(V, 2, 'omitnan');           % subtract mean over time for each voxel
v = var - mean(var, 'omitnan');          % demean variable

% ----- Compute Pearson correlation for each voxel -----
% r = sum(V.*v,2) / sqrt(sum(V.^2,2)*sum(v.^2))
numerator   = V * v;                     % [nVoxels x 1]
denom_vox   = sqrt(sum(V.^2, 2));        % [nVoxels x 1]
denom_var   = sqrt(sum(v.^2));           % scalar
denominator = denom_vox * denom_var;     % [nVoxels x 1]
r = numerator ./ denominator;           % [nVoxels x 1]

% Handle degenerate cases (zero variance)
r(denominator == 0) = NaN;
r(allNaN) = NaN;

% ----- Reshape back to [nx x ny x nz] -----
corr_vol = reshape(r, [nx ny nz]);

% ----- Plot each z-plane -----
figure('Name',['Pixelwise corr with ' var_name], 'Color','w');
nCols = ceil(sqrt(nz));
nRows = ceil(nz / nCols);

for z = 1:nz
    subplot(nRows, nCols, z);
    imagesc(corr_vol(:,:,z), [-0.2 0.2]);  % fix scale to [-1 1]
    axis image off;
    title(sprintf('z = %d', z));
    colormap(gca, jet);               % per-axes colormap
    colorbar;
end

sgtitle(sprintf('Pixelwise correlation with %s', var_name));

% ----- Save figure if requested -----
if ~isempty(savepath)
    save_plot_with_title_as_filename(['corr_with_' var_name], 'volume', savepath);
end
end