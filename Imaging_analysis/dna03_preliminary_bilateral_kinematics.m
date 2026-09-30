function results = dna03_preliminary_bilateral_kinematics(dataDir, outputDir)
% DNA03_PRELIMINARY_BILATERAL_KINEMATICS
% Preliminary analysis of bilateral DNa03 calcium activity and behavior.
%
% This pipeline is intentionally self-contained and heavily commented so it
% can serve as a readable starting point for subsequent DNa03 analyses.
%
% For every *_pData.mat file, the pipeline:
%   1. Loads the preprocessed left and right dF/F traces.
%   2. Interpolates FicTrac kinematics onto the imaging time base.
%   3. Computes zero-lag Pearson correlations between each fluorescence
%      representation and each kinematic variable.
%   4. Computes lagged correlations for forward velocity and yaw velocity.
%   5. Masks brief, extreme fluorescence spikes from analysis copies while
%      leaving the source pData files unchanged.
%   6. Compares kinematics when both copies are jointly high versus jointly
%      low and aligns behavior to entries into those common-mode states.
%   7. Finds windows with sliding left-right correlation r >= 0.6 and
%      asks which kinematics accompany unusually coordinated activity.
%   8. Plots the same flat path colored by left dF/F, right dF/F, and yaw
%      speed for direct spatial comparison.
%   9. Regresses signed yaw and yaw speed out of each side, then correlates
%      the turn-residualized fluorescence with every kinematic variable.
%  10. Plots 4-by-4 dF/F-versus-kinematic scatter matrices pooled across
%      all flies and separately for each biological fly, with linear fits.
%  11. Saves per-trial diagnostic figures and dataset-level CSV summaries.
%
% Fluorescence representations:
%   left, right          Individual DNa03 copies.
%   bilateral_mean      Mean of left and right; common-mode activity.
%   right_minus_left    Right minus left; lateralized activity.
%
% Kinematic variables:
%   forward_velocity, lateral_velocity, yaw_velocity, yaw_speed,
%   total_speed.
%
% USAGE
%   results = dna03_preliminary_bilateral_kinematics;
%   results = dna03_preliminary_bilateral_kinematics(dataDir, outputDir);
%
% OUTPUTS
%   <outputDir>/trial_correlations.csv
%   <outputDir>/trial_quality_control.csv
%   <outputDir>/dataset_correlation_summary.csv
%   <outputDir>/common_mode_state_metrics.csv
%   <outputDir>/rolling_correlation_metrics.csv
%   <outputDir>/turn_residualized_correlations.csv
%   <outputDir>/turn_residualized_correlation_summary.csv
%   <outputDir>/figures/<trial>_overview.png
%   <outputDir>/figures/<trial>_lag_correlations.png
%   <outputDir>/figures/<trial>_common_mode_states.png
%   <outputDir>/figures/<trial>_rolling_bilateral_correlation.png
%   <outputDir>/figures/<trial>_flat_paths_left_right_yawspeed.png
%   <outputDir>/figures/<trial>_flat_path_left_right_correlation.png
%
% Notes:
% - The script uses img.filtDFF rather than raw img.dFF because these files
%   already provide a standardized filtered signal.
% - Pearson correlations are descriptive. Serial autocorrelation means the
%   nominal sample count should not be used directly for significance tests.
% - Positive lag means that dF/F is shifted later than behavior; therefore,
%   a peak at positive lag is consistent with neural activity following the
%   behavioral variable.
% - Wilson-lab helper code in findYawVelPeaksFT.m defines positive yaw as a
%   right turn and negative yaw as a left turn. Verify this convention with
%   a calibration trial before making a final anatomical interpretation.

if nargin < 1 || isempty(dataDir)
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'pData_DNa03');
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(dataDir, 'preliminary_bilateral_analysis');
end

figureDir = fullfile(outputDir, 'figures');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
if ~exist(figureDir, 'dir'), mkdir(figureDir); end

files = dir(fullfile(dataDir, '*_pData.mat'));
if isempty(files)
    error('No *_pData.mat files found in %s', dataDir);
end

% Analysis settings are collected here so they are easy to find and change.
maxLagSeconds = 5;
minimumValidSamples = 100;
overviewDurationSeconds = 120;
highCorrelationThreshold = 0.6;
rollingCorrelationWindowSeconds = 1;
minimumStateDurationSeconds = 0.2;

correlationRows = table();
qcRows = table();
stateRows = table();
rollingCorrelationRows = table();
turnResidualRows = table();
scatterTrialData = struct('fly_id',{},'trial_id',{},'dff',{},'kinematics',{});
errorLog = strings(0, 2);

fprintf('DNa03 preliminary analysis: %d files\n', numel(files));
for fileIndex = 1:numel(files)
    filePath = fullfile(files(fileIndex).folder, files(fileIndex).name);
    [~, trialID] = fileparts(filePath);
    trialID = erase(trialID, '_pData');
    fprintf('[%d/%d] %s\n', fileIndex, numel(files), trialID);

    try
        loaded = load(filePath, 'img', 'fictrac', 'name', 'exptCond');
        validate_trial(loaded, filePath);

        % Put all signals into column-vector form.
        imageTime = double(loaded.img.t(:));
        leftDFF = double(loaded.img.filtDFF.left(:));
        rightDFF = double(loaded.img.filtDFF.right(:));
        fictracTime = double(loaded.fictrac.t(:));

        % Crop imaging variables conservatively if a file contains a minor
        % endpoint mismatch. The removed count is recorded in the QC table.
        nImage = min([numel(imageTime), numel(leftDFF), numel(rightDFF)]);
        croppedImageSamples = max([numel(imageTime), numel(leftDFF), ...
            numel(rightDFF)]) - nImage;
        imageTime = imageTime(1:nImage);
        leftDFF = leftDFF(1:nImage);
        rightDFF = rightDFF(1:nImage);

        % Remove only brief, extreme transients from analysis copies. This
        % does not alter img.filtDFF in the source MAT file. A point is
% flagged when dF/F exceeds 10 or either its high-pass residual or
% frame-to-frame jump exceeds a conservative 25 robust SD; the mask is expanded by
        % 30 ms on either side to remove the immediate spike footprint.
        medianImagingDt = median(diff(imageTime), 'omitnan');
        [leftDFF, leftArtifactMask] = mask_fluorescence_spikes( ...
            leftDFF, medianImagingDt);
        [rightDFF, rightArtifactMask] = mask_fluorescence_spikes( ...
            rightDFF, medianImagingDt);
        artifactMask = leftArtifactMask | rightArtifactMask;
        leftDFF(artifactMask) = NaN;
        rightDFF(artifactMask) = NaN;

        % Interpolate each FicTrac variable onto the imaging timestamps.
        % NaN is used outside the measured FicTrac interval so extrapolated
        % behavior can never contribute to a correlation.
        behaviorNames = {'forward_velocity', 'lateral_velocity', ...
            'yaw_velocity', 'yaw_speed', 'total_speed'};
        behaviorFields = {'fwdVel', 'slideVel', 'yawAngVel', ...
            'yawAngSpd', 'totSpd'};
        behavior = nan(nImage, numel(behaviorFields));
        for behaviorIndex = 1:numel(behaviorFields)
            rawBehavior = double(loaded.fictrac.(behaviorFields{behaviorIndex})(:));
            nFictrac = min(numel(fictracTime), numel(rawBehavior));
            behavior(:, behaviorIndex) = interp1( ...
                fictracTime(1:nFictrac), rawBehavior(1:nFictrac), ...
                imageTime, 'linear', NaN);
        end

        % Position uses the same interpolation and endpoint rules as the
        % kinematic variables, ensuring all three flat-path panels share
        % exactly the same samples and geometry.
        xPosition = interpolate_fictrac_field(loaded.fictrac, ...
            'xPos', fictracTime, imageTime);
        yPosition = interpolate_fictrac_field(loaded.fictrac, ...
            'yPos', fictracTime, imageTime);

        fluorescenceNames = {'left', 'right', 'bilateral_mean', ...
            'right_minus_left'};
        fluorescence = [leftDFF, rightDFF, ...
            mean([leftDFF rightDFF], 2, 'omitnan'), rightDFF-leftDFF];

        % Retain aligned samples for the requested scatter matrices. Here
        % R+L is the literal sum (not the bilateral mean), and rotation is
        % signed yaw velocity. Date is included in fly_id so fly numbers
        % reused on different experimental dates are not merged.
        scatterTrialData(end+1).fly_id = string(loaded.name.dateName) + "_" + ...
            string(loaded.name.flyName); %#ok<AGROW>
        scatterTrialData(end).trial_id = string(trialID);
        scatterTrialData(end).dff = [leftDFF rightDFF ...
            rightDFF-leftDFF rightDFF+leftDFF];
        scatterTrialData(end).kinematics = behavior(:,[1 3 2 5]);

        % Remove the contemporaneous turning component separately from the
        % left and right traces. Signed yaw captures turn direction, while
        % yaw speed captures turn magnitude. The residual therefore means
        % "fluorescence not linearly explained by turning"; source traces
        % are never modified. Correlation with yaw variables should be near
        % zero by construction and is retained as a useful sanity check.
        turnPredictors = behavior(:,[3 4]);
        [leftTurnResidual,leftTurnR2] = residualize_turning(leftDFF,turnPredictors);
        [rightTurnResidual,rightTurnR2] = residualize_turning(rightDFF,turnPredictors);
        turnResidualRows = [turnResidualRows; make_turn_residual_rows( ...
            trialID,loaded.name,leftTurnResidual,rightTurnResidual, ...
            leftTurnR2,rightTurnR2,behavior,behaviorNames)]; %#ok<AGROW>

        % Bilateral correlation is a central QC and biological measurement.
        bilateralR = pairwise_corr(leftDFF, rightDFF, minimumValidSamples);

        % Store long-format zero-lag correlations. Long format is convenient
        % for grouping by fly, FOV, date, or variable in MATLAB/R/Python.
        for fluorescenceIndex = 1:numel(fluorescenceNames)
            for behaviorIndex = 1:numel(behaviorNames)
                r = pairwise_corr(fluorescence(:, fluorescenceIndex), ...
                    behavior(:, behaviorIndex), minimumValidSamples);
                newRow = table(string(trialID), string(loaded.name.dateName), ...
                    string(loaded.name.flyName), string(loaded.name.fovName), ...
                    string(loaded.name.trialName), ...
                    string(fluorescenceNames{fluorescenceIndex}), ...
                    string(behaviorNames{behaviorIndex}), r, ...
                    'VariableNames', {'trial_id','date','fly','fov','trial', ...
                    'fluorescence_signal','kinematic_variable','pearson_r'});
                correlationRows = [correlationRows; newRow]; %#ok<AGROW>
            end
        end

        validFraction = mean(all(isfinite([fluorescence behavior]), 2));
        qcRow = table(string(trialID), nImage, numel(fictracTime), ...
            croppedImageSamples, mean(artifactMask), validFraction, ...
            medianImagingDt, bilateralR, ...
            'VariableNames', {'trial_id','n_imaging_samples', ...
            'n_fictrac_samples','cropped_imaging_samples', ...
            'artifact_fraction','valid_fraction','median_imaging_dt_s', ...
            'left_right_pearson_r'});
        qcRows = [qcRows; qcRow]; %#ok<AGROW>

        % Common-mode states ask what the animal is doing when both copies
        % are simultaneously elevated or suppressed. Each side is robustly
        % z-scored separately, preventing a brighter ROI from dominating.
        leftZ = robust_zscore(leftDFF);
        rightZ = robust_zscore(rightDFF);
        jointHigh = clean_binary_state(leftZ >= 1 & rightZ >= 1, ...
            medianImagingDt, minimumStateDurationSeconds, 0.1);
        jointLow = clean_binary_state(leftZ <= -1 & rightZ <= -1, ...
            medianImagingDt, minimumStateDurationSeconds, 0.1);

        stateRows = [stateRows; make_state_rows(trialID, loaded.name, ...
            jointHigh, jointLow, behavior)]; %#ok<AGROW>

        % A one-second sliding Pearson correlation captures coordination
        % independent of whether both signals are high or low. "Especially
        % correlated" uses the same r >= 0.6 criterion in every
        % trial, then removes fragments shorter than 200 ms.
        rollingR = rolling_pairwise_corr(leftDFF, rightDFF, ...
            medianImagingDt, rollingCorrelationWindowSeconds);
        highCorrelationState = clean_binary_state( ...
            rollingR >= highCorrelationThreshold, medianImagingDt, ...
            minimumStateDurationSeconds, 0.1);
        rollingCorrelationRows = [rollingCorrelationRows; ...
            make_correlation_state_rows(trialID, loaded.name, rollingR, ...
            highCorrelationState, highCorrelationThreshold, behavior)]; %#ok<AGROW>

        % Lag correlations are calculated for signed left/right activity,
        % their common mode, and their difference against forward and yaw.
        lagBehaviorIndices = [1 3];
        [lagsSeconds, lagR] = lag_correlations(fluorescence, ...
            behavior(:, lagBehaviorIndices), imageTime, maxLagSeconds, ...
            minimumValidSamples);

        save_overview_figure(trialID, imageTime, fluorescence, behavior, ...
            bilateralR, overviewDurationSeconds, figureDir);
        save_lag_figure(trialID, lagsSeconds, lagR, fluorescenceNames, ...
            behaviorNames(lagBehaviorIndices), figureDir);
        save_common_mode_figure(trialID, imageTime, leftZ, rightZ, ...
            jointHigh, jointLow, behavior, figureDir);
        save_rolling_correlation_figure(trialID, imageTime, rollingR, ...
            highCorrelationState, highCorrelationThreshold, behavior, figureDir);
        save_flat_path_figure(trialID, xPosition, yPosition, leftDFF, ...
            rightDFF, behavior(:,4), figureDir);
        save_correlation_flat_path_figure(trialID, xPosition, yPosition, ...
            rollingR, highCorrelationState, highCorrelationThreshold, ...
            leftZ, rightZ, figureDir);

    catch ME
        warning('Could not analyze %s: %s', files(fileIndex).name, ME.message);
        errorLog(end+1, :) = [string(files(fileIndex).name), string(ME.message)]; %#ok<AGROW>
    end
end

if isempty(correlationRows)
    error('No trials were successfully analyzed.');
end

% Dataset-level mean and between-trial SD. Individual trials remain the
% units here; this avoids weighting longer trials more heavily.
summaryRows = groupsummary(correlationRows, ...
    {'fluorescence_signal','kinematic_variable'}, {'mean','std'}, 'pearson_r');

writetable(correlationRows, fullfile(outputDir, 'trial_correlations.csv'));
writetable(qcRows, fullfile(outputDir, 'trial_quality_control.csv'));
writetable(summaryRows, fullfile(outputDir, 'dataset_correlation_summary.csv'));
writetable(stateRows, fullfile(outputDir, 'common_mode_state_metrics.csv'));
writetable(rollingCorrelationRows, ...
    fullfile(outputDir, 'rolling_correlation_metrics.csv'));
writetable(turnResidualRows, ...
    fullfile(outputDir, 'turn_residualized_correlations.csv'));
turnResidualSummary = groupsummary(turnResidualRows, ...
    {'fluorescence_side','kinematic_variable'}, {'mean','std'}, 'pearson_r');
writetable(turnResidualSummary, ...
    fullfile(outputDir, 'turn_residualized_correlation_summary.csv'));

% Scatter plots are descriptive frame-level views. Trend fits use every
% finite aligned sample; plotting is deterministically thinned only to keep
% PNG/FIG files responsive. Per-fly panels pool that fly's trials/FOVs.
scatterRegressionRows = save_dff_kinematic_scatter_matrices( ...
    scatterTrialData, figureDir);
writetable(scatterRegressionRows, ...
    fullfile(outputDir, 'dff_kinematic_scatter_regressions.csv'));

if ~isempty(errorLog)
    errorTable = array2table(errorLog, 'VariableNames', {'file','error'});
    writetable(errorTable, fullfile(outputDir, 'processing_errors.csv'));
else
    errorTable = table();
end

results = struct();
results.correlations = correlationRows;
results.quality_control = qcRows;
results.summary = summaryRows;
results.common_mode_states = stateRows;
results.rolling_correlation = rollingCorrelationRows;
results.turn_residualized_correlations = turnResidualRows;
results.turn_residualized_summary = turnResidualSummary;
results.dff_kinematic_scatter_regressions = scatterRegressionRows;
results.errors = errorTable;
results.output_directory = outputDir;

fprintf('Finished: %d/%d trials analyzed successfully.\n', ...
    height(qcRows), numel(files));
fprintf('Results saved to %s\n', outputDir);
end

function validate_trial(loaded, filePath)
% Fail early with a readable message if a required field is absent.
requiredTop = {'img','fictrac','name'};
for i = 1:numel(requiredTop)
    if ~isfield(loaded, requiredTop{i})
        error('%s is missing top-level field %s.', filePath, requiredTop{i});
    end
end
requiredImage = {'t','filtDFF'};
requiredFictrac = {'t','fwdVel','slideVel','yawAngVel','yawAngSpd','totSpd'};
requiredFictrac = [requiredFictrac {'xPos','yPos'}];
for i = 1:numel(requiredImage)
    if ~isfield(loaded.img, requiredImage{i})
        error('%s is missing img.%s.', filePath, requiredImage{i});
    end
end

if ~isfield(loaded.img.filtDFF,'left') || ~isfield(loaded.img.filtDFF,'right')
    error('%s does not contain both img.filtDFF.left and .right.', filePath);
end
for i = 1:numel(requiredFictrac)
    if ~isfield(loaded.fictrac, requiredFictrac{i})
        error('%s is missing fictrac.%s.', filePath, requiredFictrac{i});
    end
end
end

function aligned = interpolate_fictrac_field(fictrac, fieldName, ...
    fictracTime, imageTime)
% Interpolate one FicTrac field without extrapolating beyond recorded data.
raw = double(fictrac.(fieldName)(:));
n = min(numel(fictracTime),numel(raw));
aligned = interp1(fictracTime(1:n),raw(1:n),imageTime,'linear',NaN);
end

function r = pairwise_corr(x, y, minimumValidSamples)
% Pearson correlation after pairwise removal of NaN/Inf samples.
valid = isfinite(x) & isfinite(y);
if nnz(valid) < minimumValidSamples || std(x(valid)) == 0 || std(y(valid)) == 0
    r = NaN;
    return;
end
c = corrcoef(x(valid), y(valid));
r = c(1,2);
end

function [residual, modelR2] = residualize_turning(fluorescence, turnPredictors)
% Regress dF/F on [signed yaw velocity, yaw speed] with an intercept.
% Fitting and residual calculation use only rows where all three signals
% are finite. NaNs remain NaN, preventing interpolation across artifacts.
fluorescence=fluorescence(:);
valid=isfinite(fluorescence)&all(isfinite(turnPredictors),2);
residual=nan(size(fluorescence)); modelR2=NaN;
if nnz(valid)<3, return; end
X=[ones(nnz(valid),1),turnPredictors(valid,:)];
y=fluorescence(valid);
coefficients=X\y;
fitted=X*coefficients;
residual(valid)=y-fitted;
totalSS=sum((y-mean(y)).^2);
if totalSS>0, modelR2=1-sum((y-fitted).^2)/totalSS; end
end

function rows = make_turn_residual_rows(trialID,nameStruct,leftResidual, ...
    rightResidual,leftR2,rightR2,behavior,behaviorNames)
% Store correlations after removing each side's linear turning component.
residuals={leftResidual,rightResidual}; sides={'left','right'};
modelR2=[leftR2 rightR2]; rows=table();
for sideIndex=1:2
    for behaviorIndex=1:numel(behaviorNames)
        r=pairwise_corr(residuals{sideIndex},behavior(:,behaviorIndex),100);
        newRow=table(string(trialID),string(nameStruct.dateName), ...
            string(nameStruct.flyName),string(nameStruct.fovName), ...
            string(nameStruct.trialName),string(sides{sideIndex}), ...
            modelR2(sideIndex),string(behaviorNames{behaviorIndex}),r, ...
            'VariableNames',{'trial_id','date','fly','fov','trial', ...
            'fluorescence_side','turn_model_r_squared', ...
            'kinematic_variable','pearson_r'});
        rows=[rows;newRow]; %#ok<AGROW>
    end
end
end

function [cleanTrace, artifactMask] = mask_fluorescence_spikes(trace, dt)
% Identify short fluorescence artifacts using robust amplitude and slope.
% The local median removes slow biological fluctuations before thresholding.
trace = trace(:);
windowSamples = max(3, round(2/dt));
if mod(windowSamples,2) == 0, windowSamples = windowSamples + 1; end
localBaseline = movmedian(trace, windowSamples, 'omitnan');
residual = trace - localBaseline;

residualScale = robust_scale(residual);
traceDerivative = [0; diff(trace)];
derivativeScale = robust_scale(traceDerivative);

artifactCore = trace > 10 | abs(residual) > 25*residualScale | ...
    abs(traceDerivative) > 25*derivativeScale;

% Expand the core detection so the rising and falling edges of each spike
% are removed together. Convolution avoids requiring an extra toolbox.
paddingSamples = max(1, round(0.03/dt));
artifactMask = conv(double(artifactCore), ...
    ones(2*paddingSamples+1,1), 'same') > 0;
artifactMask = artifactMask | ~isfinite(trace);

cleanTrace = trace;
cleanTrace(artifactMask) = NaN;
end

function scale = robust_scale(x)
% Median absolute deviation converted to a Gaussian-equivalent SD.
x = x(isfinite(x));
if isempty(x)
    scale = NaN;
    return;
end
center = median(x);
scale = 1.4826 * median(abs(x-center));
if ~isfinite(scale) || scale <= eps
    scale = std(x);
end
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
end

function z = robust_zscore(x)
% Robust z-score used only to define simultaneous high and low states.
valid = isfinite(x);
z = nan(size(x));
if ~any(valid), return; end
center = median(x(valid));
scale = robust_scale(x(valid));
z(valid) = (x(valid)-center)/scale;
end

function cleaned = clean_binary_state(mask, dt, minimumDuration, maximumGap)
% Merge brief gaps and reject threshold crossings that are too short to be
% interpreted as a calcium-activity state.
cleaned = logical(mask(:));
maximumGapSamples = max(0, round(maximumGap/dt));
minimumStateSamples = max(1, round(minimumDuration/dt));

% Fill short false gaps only when bounded by true samples on both sides.
falseMask = ~cleaned;
[gapStarts,gapEnds] = logical_runs(falseMask);
for i = 1:numel(gapStarts)
    gapLength = gapEnds(i)-gapStarts(i)+1;
    bounded = gapStarts(i)>1 && gapEnds(i)<numel(cleaned) && ...
        cleaned(gapStarts(i)-1) && cleaned(gapEnds(i)+1);
    if bounded && gapLength <= maximumGapSamples
        cleaned(gapStarts(i):gapEnds(i)) = true;
    end
end

% Remove true bouts that do not meet the minimum duration.
[boutStarts,boutEnds] = logical_runs(cleaned);
for i = 1:numel(boutStarts)
    if boutEnds(i)-boutStarts(i)+1 < minimumStateSamples
        cleaned(boutStarts(i):boutEnds(i)) = false;
    end
end
end

function [runStarts,runEnds] = logical_runs(mask)
% Return inclusive start/end indices for all true runs in a logical vector.
edges = diff([false; logical(mask(:)); false]);
runStarts = find(edges==1);
runEnds = find(edges==-1)-1;
end

function rows = make_state_rows(trialID, nameStruct, jointHigh, jointLow, behavior)
% Summarize occupancy and kinematics in the two common-mode states.
stateNames = {'joint_high','joint_low'};
stateMasks = {jointHigh,jointLow};
kinematicNames = {'forward_velocity','lateral_velocity','yaw_velocity', ...
    'yaw_speed','total_speed'};
rows = table();
for stateIndex = 1:numel(stateNames)
    mask = stateMasks{stateIndex};
    for behaviorIndex = 1:numel(kinematicNames)
        values = behavior(mask,behaviorIndex);
        values = values(isfinite(values));
        if isempty(values)
            meanValue = NaN; medianValue = NaN;
        else
            meanValue = mean(values); medianValue = median(values);
        end
        newRow = table(string(trialID), string(nameStruct.dateName), ...
            string(nameStruct.flyName), string(nameStruct.fovName), ...
            string(nameStruct.trialName), string(stateNames{stateIndex}), ...
            string(kinematicNames{behaviorIndex}), mean(mask), nnz(mask), ...
            meanValue, medianValue, ...
            'VariableNames', {'trial_id','date','fly','fov','trial','state', ...
            'kinematic_variable','state_fraction','n_samples', ...
            'mean_kinematic_value','median_kinematic_value'});
        rows = [rows; newRow]; %#ok<AGROW>
    end
end
end

function rollingR = rolling_pairwise_corr(x, y, dt, windowSeconds)
% Sliding Pearson correlation using a centered window and pairwise-valid data.
x = x(:); y = y(:);
windowSamples = max(5,round(windowSeconds/dt));
if mod(windowSamples,2)==0, windowSamples=windowSamples+1; end
halfWindow = floor(windowSamples/2);
window = [halfWindow halfWindow];
valid = isfinite(x) & isfinite(y);
x0 = x; y0 = y; x0(~valid)=0; y0(~valid)=0;
n = movsum(double(valid),window,'Endpoints','shrink');
sumX = movsum(x0,window,'Endpoints','shrink');
sumY = movsum(y0,window,'Endpoints','shrink');
sumXX = movsum(x0.^2,window,'Endpoints','shrink');
sumYY = movsum(y0.^2,window,'Endpoints','shrink');
sumXY = movsum(x0.*y0,window,'Endpoints','shrink');
covXY = sumXY - (sumX.*sumY)./max(n,1);
varX = sumXX - (sumX.^2)./max(n,1);
varY = sumYY - (sumY.^2)./max(n,1);
rollingR = covXY ./ sqrt(varX.*varY);
rollingR(n < 0.5*windowSamples | varX<=0 | varY<=0) = NaN;
rollingR = max(-1,min(1,rollingR));
end

function rows = make_correlation_state_rows(trialID, nameStruct, rollingR, ...
    highState, threshold, behavior)
% Compare kinematics in the most strongly correlated windows with all other
% windows having a valid rolling correlation estimate.
kinematicNames = {'forward_velocity','lateral_velocity','yaw_velocity', ...
    'yaw_speed','total_speed'};
stateNames = {'r_ge_0_8','other_valid'};
stateMasks = {highState, isfinite(rollingR) & ~highState};
rows = table();
for stateIndex = 1:numel(stateNames)
    mask = stateMasks{stateIndex};
    for behaviorIndex = 1:numel(kinematicNames)
        values = behavior(mask,behaviorIndex);
        values = values(isfinite(values));
        if isempty(values)
            meanValue=NaN; medianValue=NaN;
        else
            meanValue=mean(values); medianValue=median(values);
        end
        newRow = table(string(trialID),string(nameStruct.dateName), ...
            string(nameStruct.flyName),string(nameStruct.fovName), ...
            string(nameStruct.trialName),string(stateNames{stateIndex}), ...
            threshold,mean(mask),nnz(mask),mean(rollingR(mask),'omitnan'), ...
            string(kinematicNames{behaviorIndex}),meanValue,medianValue, ...
            'VariableNames',{'trial_id','date','fly','fov','trial', ...
            'coordination_state','correlation_threshold','state_fraction', ...
            'n_samples','mean_rolling_r','kinematic_variable', ...
            'mean_kinematic_value','median_kinematic_value'});
        rows=[rows;newRow]; %#ok<AGROW>
    end
end
end

function [lagsSeconds, lagR] = lag_correlations(fluorescence, behavior, ...
    imageTime, maxLagSeconds, minimumValidSamples)
% Compute Pearson r at integer imaging-frame offsets.
dt = median(diff(imageTime), 'omitnan');
maxLagSamples = max(1, round(maxLagSeconds/dt));
lagSamples = (-maxLagSamples:maxLagSamples)';
lagsSeconds = lagSamples * dt;
lagR = nan(numel(lagSamples), size(fluorescence,2), size(behavior,2));

for lagIndex = 1:numel(lagSamples)
    lag = lagSamples(lagIndex);
    if lag >= 0
        fluorescenceIndices = (1+lag):size(fluorescence,1);
        behaviorIndices = 1:(size(fluorescence,1)-lag);
    else
        fluorescenceIndices = 1:(size(fluorescence,1)+lag);
        behaviorIndices = (1-lag):size(fluorescence,1);
    end
    for fluorescenceIndex = 1:size(fluorescence,2)
        for behaviorIndex = 1:size(behavior,2)
            lagR(lagIndex, fluorescenceIndex, behaviorIndex) = pairwise_corr( ...
                fluorescence(fluorescenceIndices, fluorescenceIndex), ...
                behavior(behaviorIndices, behaviorIndex), minimumValidSamples);
        end
    end
end
end

function save_overview_figure(trialID, imageTime, fluorescence, behavior, ...
    bilateralR, overviewDurationSeconds, figureDir)
% Save a compact trace overview to catch alignment or signal-quality issues.
startTime = imageTime(1);
keep = imageTime <= min(imageTime(end), startTime + overviewDurationSeconds);
t = imageTime(keep) - startTime;

fig = figure('Visible','off','Color','white','Position',[100 100 1300 900]);
tiledlayout(fig, 4, 1, 'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(t, fluorescence(keep,1), 'Color',[0 0.45 0.74], 'DisplayName','left');
plot(t, fluorescence(keep,2), 'Color',[0.85 0.33 0.10], 'DisplayName','right');
ylabel('dF/F'); title(sprintf('%s | left vs right (r = %.3f)', trialID, bilateralR), ...
    'Interpreter','none'); legend('Location','best'); grid on;

nexttile; plot(t, fluorescence(keep,4), 'Color',[0.49 0.18 0.56]);
ylabel('R - L dF/F'); title('Lateralized fluorescence'); grid on;

nexttile; plot(t, behavior(keep,1), 'k');
ylabel('Forward velocity'); title('Forward velocity'); grid on;

nexttile; plot(t, behavior(keep,3), 'Color',[0.64 0.08 0.18]);
ylabel('Yaw velocity'); xlabel('Time from trial start (s)');
title('Signed yaw velocity'); grid on;

exportgraphics(fig, fullfile(figureDir, [trialID '_overview.png']), ...
    'Resolution',150,'BackgroundColor','white');
savefig(fig,fullfile(figureDir,[trialID '_overview.fig']),'compact');
close(fig);
end

function save_lag_figure(trialID, lagsSeconds, lagR, fluorescenceNames, ...
    behaviorNames, figureDir)
% Plot lag-dependent correlation curves for forward and yaw velocity.
colors = lines(numel(fluorescenceNames));
fig = figure('Visible','off','Color','white','Position',[100 100 1200 500]);
tiledlayout(fig, 1, numel(behaviorNames), 'TileSpacing','compact','Padding','compact');
for behaviorIndex = 1:numel(behaviorNames)
    nexttile; hold on;
    for fluorescenceIndex = 1:numel(fluorescenceNames)
        plot(lagsSeconds, lagR(:,fluorescenceIndex,behaviorIndex), ...
            'LineWidth',1.5,'Color',colors(fluorescenceIndex,:), ...
            'DisplayName',strrep(fluorescenceNames{fluorescenceIndex},'_',' '));
    end
    xline(0,':k','HandleVisibility','off');
    yline(0,':k','HandleVisibility','off');
    xlabel('Lag: dF/F relative to behavior (s)'); ylabel('Pearson r');
    title(strrep(behaviorNames{behaviorIndex},'_',' '));
    if behaviorIndex == 1, legend('Location','best'); end
end
sgtitle([trialID ' | lagged correlations'], 'Interpreter','none');
exportgraphics(fig, fullfile(figureDir, [trialID '_lag_correlations.png']), ...
    'Resolution',150,'BackgroundColor','white');
savefig(fig,fullfile(figureDir,[trialID '_lag_correlations.fig']),'compact');
close(fig);
end

function save_common_mode_figure(trialID, imageTime, leftZ, rightZ, ...
    jointHigh, jointLow, behavior, figureDir)
% Visualize simultaneous high/low states and behavior around state entries.
% Entry-triggered means are descriptive; the number of entries is shown so
% sparse states are immediately recognizable.
dt = median(diff(imageTime), 'omitnan');
halfWindow = max(1, round(3/dt));
relativeTime = (-halfWindow:halfWindow)' * dt;

highEntries = find(jointHigh & ~[false; jointHigh(1:end-1)]);
lowEntries = find(jointLow & ~[false; jointLow(1:end-1)]);
highEntries = highEntries(highEntries>halfWindow & ...
    highEntries<=numel(imageTime)-halfWindow);
lowEntries = lowEntries(lowEntries>halfWindow & ...
    lowEntries<=numel(imageTime)-halfWindow);

% Use forward velocity, signed yaw, yaw speed, and total speed. Lateral
% velocity remains available in the CSV but is omitted to keep this figure
% readable.
behaviorIndices = [1 3 4 5];
behaviorLabels = {'Forward velocity','Yaw velocity (+ = right)', ...
    'Yaw speed','Total speed'};
highTriggered = triggered_mean(behavior(:,behaviorIndices), highEntries, halfWindow);
lowTriggered = triggered_mean(behavior(:,behaviorIndices), lowEntries, halfWindow);

fig = figure('Visible','off','Color','white','Position',[100 100 1500 850]);
tiledlayout(fig, 2, 4, 'TileSpacing','compact','Padding','compact');

% A short trace segment confirms that the state definitions behave as
% intended and makes residual artifacts easy to spot.
nexttile([1 3]); hold on;
duration = min(120, imageTime(end)-imageTime(1));
keep = imageTime <= imageTime(1)+duration;
t = imageTime(keep)-imageTime(1);
plot(t,leftZ(keep),'Color',[0 0.45 0.74],'DisplayName','left z');
plot(t,rightZ(keep),'Color',[0.85 0.33 0.10],'DisplayName','right z');
plot(t,2.5*double(jointHigh(keep)),'Color',[0.47 0.67 0.19], ...
    'DisplayName','both increased to high');
plot(t,-2.5*double(jointLow(keep)),'Color',[0.49 0.18 0.56], ...
    'DisplayName','both dropped to low');
yline(1,':','HandleVisibility','off'); yline(-1,':','HandleVisibility','off');
xlabel('Time (s)'); ylabel('Robust z / state'); grid on; legend('Location','best');
title('Simultaneous high and low bilateral activity');

nexttile;
bar([mean(jointHigh),mean(jointLow)]*100, 'FaceColor',[0.3 0.5 0.7]);
set(gca,'XTickLabel',{'both high','both low'}); ylabel('Trial occupancy (%)');
title(sprintf('high n=%d entries; low n=%d',numel(highEntries),numel(lowEntries)));

for behaviorIndex = 1:numel(behaviorIndices)
    nexttile; hold on;
    plot(relativeTime, highTriggered(:,behaviorIndex), ...
        'Color',[0.47 0.67 0.19],'LineWidth',2,'DisplayName','both increase');
    plot(relativeTime, lowTriggered(:,behaviorIndex), ...
        'Color',[0.49 0.18 0.56],'LineWidth',2,'DisplayName','both decrease');
    xline(0,':k','HandleVisibility','off'); grid on;
    xlabel('Time from state entry (s)'); ylabel(behaviorLabels{behaviorIndex});
    if behaviorIndex == 1, legend('Location','best'); end
end
sgtitle([trialID ' | common-mode activity states'], 'Interpreter','none');
exportgraphics(fig, fullfile(figureDir, [trialID '_common_mode_states.png']), ...
    'Resolution',150,'BackgroundColor','white');
savefig(fig,fullfile(figureDir,[trialID '_common_mode_states.fig']),'compact');
close(fig);
end

function meanTrace = triggered_mean(signals, eventIndices, halfWindow)
% Extract fixed windows and return an event-triggered mean for each signal.
windowLength = 2*halfWindow+1;
meanTrace = nan(windowLength,size(signals,2));
if isempty(eventIndices), return; end
stack = nan(windowLength,size(signals,2),numel(eventIndices));
for eventIndex = 1:numel(eventIndices)
    idx = eventIndices(eventIndex)-halfWindow:eventIndices(eventIndex)+halfWindow;
    stack(:,:,eventIndex) = signals(idx,:);
end
meanTrace = mean(stack,3,'omitnan');
end

function save_rolling_correlation_figure(trialID, imageTime, rollingR, ...
    highState, threshold, behavior, figureDir)
% Show unusually coordinated bilateral periods and behavior at their onset.
dt = median(diff(imageTime),'omitnan');
halfWindow = max(1,round(3/dt));
relativeTime = (-halfWindow:halfWindow)'*dt;
entries = find(highState & ~[false; highState(1:end-1)]);
entries = entries(entries>halfWindow & entries<=numel(imageTime)-halfWindow);
behaviorIndices = [1 3 4 5];
labels = {'Forward velocity','Yaw velocity (+ = right)','Yaw speed','Total speed'};
eventMean = triggered_mean(behavior(:,behaviorIndices),entries,halfWindow);

fig=figure('Visible','off','Color','white','Position',[100 100 1500 850]);
tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
nexttile([1 3]); hold on;
keep=imageTime<=imageTime(1)+min(120,imageTime(end)-imageTime(1));
t=imageTime(keep)-imageTime(1);
plot(t,rollingR(keep),'Color',[0.2 0.2 0.2],'LineWidth',1);
yline(threshold,'--','Color',[0.85 0.33 0.10], ...
    'Label',sprintf('r = %.1f cutoff',threshold));
area(t,double(highState(keep)),'FaceColor',[0.47 0.67 0.19], ...
    'FaceAlpha',0.22,'EdgeColor','none');
yline(0,':k'); ylim([-1.05 1.05]); grid on;
xlabel('Time (s)'); ylabel('1-s rolling Pearson r');
title('Local left-right coordination');
nexttile;
bar(mean(highState)*100,'FaceColor',[0.47 0.67 0.19]);
xlim([0.4 1.6]); set(gca,'XTick',[]); ylabel('Trial occupancy (%)');
title(sprintf('r >= %.1f; n=%d entries',threshold,numel(entries)));
for k=1:numel(behaviorIndices)
    nexttile; plot(relativeTime,eventMean(:,k),'LineWidth',2, ...
        'Color',[0.47 0.67 0.19]); hold on; xline(0,':k'); grid on;
    xlabel('Time from state entry (s)'); ylabel(labels{k});
end
sgtitle([trialID ' | unusually correlated left-right activity'],'Interpreter','none');
exportgraphics(fig,fullfile(figureDir,[trialID '_rolling_bilateral_correlation.png']), ...
    'Resolution',150,'BackgroundColor','white');
savefig(fig,fullfile(figureDir,[trialID '_rolling_bilateral_correlation.fig']),'compact');
close(fig);
end

function save_flat_path_figure(trialID,x,y,leftDFF,rightDFF,yawSpeed,figureDir)
% Plot identical path samples in every panel. The common-mode panel sums
% separately robust-z-scored sides, preventing the brighter/more variable
% ROI from dominating R+L. Display clipping does not change measurements.
valid=isfinite(x)&isfinite(y)&isfinite(leftDFF)&isfinite(rightDFF)&isfinite(yawSpeed);
x=x(valid); y=y(valid);
leftValues=leftDFF(valid); rightValues=rightDFF(valid);
commonModeZ=robust_zscore(leftValues)+robust_zscore(rightValues);
signals={leftValues,rightValues,commonModeZ,yawSpeed(valid)};
labels={'Left dF/F','Right dF/F','z(Right) + z(Left)','Yaw speed'};
if isempty(x), return; end
xLimits=[min(x) max(x)]; yLimits=[min(y) max(y)];
if diff(xLimits)==0, xLimits=xLimits+[-1 1]; end
if diff(yLimits)==0, yLimits=yLimits+[-1 1]; end
fig=figure('Visible','off','Color','white','Position',[100 100 1900 500]);
tiledlayout(fig,1,4,'TileSpacing','compact','Padding','compact');
for k=1:4
    ax=nexttile; values=signals{k};
    if k==3
        maxAbs=prctile(abs(values),98);
        limits=[-maxAbs maxAbs];
    else
        limits=prctile(values,[2 98]);
    end
    if ~all(isfinite(limits)) || limits(1)==limits(2)
        limits=[min(values) max(values)];
        if limits(1)==limits(2), limits=limits(1)+[-1 1]; end
    end
    displayed=max(limits(1),min(limits(2),values));
    scatter(x,y,8,displayed,'filled'); hold on;
    plot(x(1),y(1),'kp','MarkerFaceColor','w','MarkerSize',10);
    axis equal; xlim(xLimits); ylim(yLimits); grid on;
    xlabel('x position'); ylabel('y position'); title(labels{k});
    if k==3
        blueWhiteRed=[linspace(0.1,1,128)' linspace(0.3,1,128)' ones(128,1); ...
            ones(128,1) linspace(1,0.15,128)' linspace(1,0.1,128)'];
        colormap(ax,blueWhiteRed);
    else
        colormap(ax,turbo);
    end
    clim(limits); colorbar;
end
sgtitle([trialID ' | matched flat paths (star = start)'],'Interpreter','none');
exportgraphics(fig,fullfile(figureDir,[trialID '_flat_paths_left_right_yawspeed.png']), ...
    'Resolution',180,'BackgroundColor','white');
savefig(fig,fullfile(figureDir,[trialID '_flat_paths_left_right_yawspeed.fig']),'compact');
close(fig);
end

function save_correlation_flat_path_figure(trialID,x,y,rollingR,highState,threshold,leftZ,rightZ,figureDir)
% Plot local bilateral coordination and its activity level side by side.
% rollingR is the centered one-second left-right Pearson correlation used
% elsewhere in the analysis. The second panel retains the full path in gray
% and colors high-correlation samples by mean(left z-score,right z-score).
% The threshold is shared with the state analysis (currently r >= 0.6). Red
% means jointly elevated, blue means jointly suppressed, and white means
% highly correlated but close to baseline. Correlation color limits span
% each trial's observed minimum and maximum; bilateral z limits stay fixed.
valid=isfinite(x)&isfinite(y)&isfinite(rollingR)&isfinite(leftZ)&isfinite(rightZ);
x=x(valid); y=y(valid); rollingR=rollingR(valid);
highState=highState(valid); bilateralZ=(leftZ(valid)+rightZ(valid))/2;
if isempty(x), return; end
fig=figure('Visible','off','Color','white','Position',[100 100 1500 700]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1=nexttile; scatter(x,y,10,rollingR,'filled'); hold on;
plot(x(1),y(1),'kp','MarkerFaceColor','w','MarkerSize',11);
correlationLimits=[min(rollingR) max(rollingR)];
if correlationLimits(1)==correlationLimits(2)
    correlationLimits=correlationLimits(1)+[-0.01 0.01];
end
axis equal tight; grid on; colormap(ax1,turbo); clim(ax1,correlationLimits);
cb=colorbar; cb.Label.String='1-s rolling left-right Pearson r';
xlabel('x position'); ylabel('y position'); title('Local left-right correlation');

ax2=nexttile; plot(x,y,'Color',[0.82 0.82 0.82],'LineWidth',1); hold on;
displayedZ=max(-2,min(2,bilateralZ));
scatter(x(highState),y(highState),12,displayedZ(highState),'filled');
plot(x(1),y(1),'kp','MarkerFaceColor','w','MarkerSize',11);
axis equal tight; grid on;
blueWhiteRed=[linspace(0.1,1,128)' linspace(0.3,1,128)' ones(128,1); ...
    ones(128,1) linspace(1,0.15,128)' linspace(1,0.1,128)'];
colormap(ax2,blueWhiteRed); clim(ax2,[-2 2]);
cb=colorbar; cb.Label.String='Bilateral mean activity (z-score)';
xlabel('x position'); ylabel('y position');
title(sprintf('r >= %.1f: blue = jointly low, red = jointly high',threshold));

sgtitle([trialID ' | bilateral correlation and activity level (star = start)'], ...
    'Interpreter','none');
exportgraphics(fig,fullfile(figureDir, ...
    [trialID '_flat_path_left_right_correlation.png']), ...
    'Resolution',180,'BackgroundColor','white');
savefig(fig,fullfile(figureDir, ...
    [trialID '_flat_path_left_right_correlation.fig']),'compact');
close(fig);
end

function regressionRows = save_dff_kinematic_scatter_matrices(trialData,figureDir)
% Plot fluorescence against four kinematic variables, pooled and per fly.
% Each biological fly is date + fly ID. Every linear model is y = dF/F and
% x = kinematics. Pearson r and the fit use all pairwise-finite frames.
% Scatter rendering is thinned to at most 6000 points per panel, and axes
% show the 0.5th-99.5th percentiles so isolated behavior spikes do not make
% the central relationship unreadable. Neither fitting nor CSV values are
% clipped or downsampled.
regressionRows=table();
if isempty(trialData), return; end

dffNames={'left_dff','right_dff','right_minus_left_dff','right_plus_left_dff'};
dffLabels={'Left dF/F','Right dF/F','R - L dF/F','R + L dF/F'};
kinematicNames={'forward_velocity','signed_yaw_velocity', ...
    'lateral_velocity','total_velocity'};
kinematicLabels={'Forward velocity','Signed yaw velocity', ...
    'Lateral velocity','Total velocity'};

flyIDs=unique([trialData.fly_id],'stable');
groupIDs=["all_flies" flyIDs];
for groupIndex=1:numel(groupIDs)
    groupID=groupIDs(groupIndex);
    if groupID=="all_flies"
        use=true(1,numel(trialData)); groupType="pooled";
        displayTitle='All flies pooled';
    else
        use=[trialData.fly_id]==groupID; groupType="per_fly";
        displayTitle=strrep(char(groupID),'_',' ');
    end
    groupTrials=trialData(use);
    allDFF=vertcat(groupTrials.dff);
    allKinematics=vertcat(groupTrials.kinematics);

    fig=figure('Visible','off','Color','white','Position',[50 50 1800 1500]);
    tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for dffIndex=1:4
        for kinematicIndex=1:4
            ax=nexttile; hold on;
            x=allKinematics(:,kinematicIndex);
            y=allDFF(:,dffIndex);
            valid=isfinite(x)&isfinite(y);
            xValid=x(valid); yValid=y(valid); nValid=numel(xValid);

            if nValid>=3 && std(xValid)>0 && std(yValid)>0
                coefficients=polyfit(xValid,yValid,1);
                r=pairwise_corr(xValid,yValid,3);
                xLimits=scatter_axis_limits(xValid);
                yLimits=scatter_axis_limits(yValid);
                sampleCount=min(6000,nValid);
                sampleIndices=unique(round(linspace(1,nValid,sampleCount)));
                scatter(xValid(sampleIndices),yValid(sampleIndices),7, ...
                    'filled','MarkerFaceColor',[0.25 0.35 0.48], ...
                    'MarkerFaceAlpha',0.16,'MarkerEdgeAlpha',0.08);
                fitX=linspace(xLimits(1),xLimits(2),100);
                plot(fitX,polyval(coefficients,fitX),'Color',[0.78 0.12 0.12], ...
                    'LineWidth',2.2);
                xlim(xLimits); ylim(yLimits);
                text(0.04,0.93,sprintf('r = %.3f | n = %d',r,nValid), ...
                    'Units','normalized','VerticalAlignment','top', ...
                    'FontSize',10,'FontWeight','bold');
                slope=coefficients(1); intercept=coefficients(2);
            else
                r=NaN; slope=NaN; intercept=NaN;
                text(0.5,0.5,'Insufficient variable data','Units','normalized', ...
                    'HorizontalAlignment','center');
            end
            grid on; box off;
            if dffIndex==1, title(kinematicLabels{kinematicIndex}); end
            if dffIndex==4, xlabel(kinematicLabels{kinematicIndex}); end
            if kinematicIndex==1, ylabel(dffLabels{dffIndex}); end

            newRow=table(groupType,string(groupID),string(dffNames{dffIndex}), ...
                string(kinematicNames{kinematicIndex}),nValid,r,slope,intercept, ...
                'VariableNames',{'group_type','group_id','dff_signal', ...
                'kinematic_variable','n_samples','pearson_r','linear_slope', ...
                'linear_intercept'});
            regressionRows=[regressionRows;newRow]; %#ok<AGROW>
        end
    end
    sgtitle(sprintf('%s | dF/F versus kinematics',displayTitle), ...
        'FontWeight','bold','Interpreter','none');
    fileStem=regexprep(char(groupID),'[^A-Za-z0-9_-]','_');
    exportgraphics(fig,fullfile(figureDir, ...
        [fileStem '_dff_kinematic_scatter_matrix.png']), ...
        'Resolution',180,'BackgroundColor','white');
    savefig(fig,fullfile(figureDir, ...
        [fileStem '_dff_kinematic_scatter_matrix.fig']),'compact');
    close(fig);
end
end

function limits = scatter_axis_limits(values)
% Robust display limits only; statistical fits continue to use all values.
limits=prctile(values,[0.5 99.5]);
if ~all(isfinite(limits)) || limits(1)==limits(2)
    limits=[min(values) max(values)];
end
if limits(1)==limits(2), limits=limits(1)+[-1 1]; end
padding=0.04*diff(limits);
limits=limits+[-padding padding];
end
