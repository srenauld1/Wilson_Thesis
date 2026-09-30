%% VES041 spontaneous behavior correlation audit
% Quantifies VES041 relationships to forward speed, slowing, and pivot turns
% in closed-loop/dark trials without odor or heat. This script is read-only:
% it loads existing fly_data.mat files and writes only a CSV summary.

trialFolders = {
    '/Users/sophiarenauld/stacks/20251001-21_ves041_fly3_clbbar'
    '/Users/sophiarenauld/stacks/20260213-2_ves041_clbar'
    '/Users/sophiarenauld/stacks/20260213-4_ves041_fly2_clbar'
    '/Users/sophiarenauld/stacks/20260217-2_ves041_vt_clbbar1'
    '/Users/sophiarenauld/stacks/20260217-3_ves041_vt_clbbar2_actuallydark'
    '/Users/sophiarenauld/stacks/20260217-4_ves041_vt_clbbar3'
    };

rows = cell(numel(trialFolders), 14);

for k = 1:numel(trialFolders)
    trialFolder = trialFolders{k};
    loaded = load(fullfile(trialFolder, 'smr_analysis', 'fly_data.mat'), 'fly_data');
    data = loaded.fly_data;

    % The processed behavior and fluorescence vectors share the imaging time
    % base after process_fictrac_panels_2p.
    t = data.dq(2).t(:);
    fwd = data.dq(2).bvf(:);
    yaw = data.dq(2).bvy_deg(:);
    dff = data.roi.ts(:);
    n = min([numel(t), numel(fwd), numel(yaw), numel(dff)]);
    t = t(1:n); fwd = fwd(1:n); yaw = yaw(1:n); dff = dff(1:n);

    % Retain the complete finite dF/F range, including values above 10.
    valid = isfinite(t) & isfinite(fwd) & isfinite(yaw) & isfinite(dff);
    dt = median(diff(t(valid)));
    fs = 1 / dt;
    smoothN = max(3, round(0.25 * fs));
    fwdSm = smoothdata(fwd, 'movmedian', smoothN, 'omitmissing');
    yawSm = smoothdata(yaw, 'movmedian', smoothN, 'omitmissing');
    dffSm = smoothdata(dff, 'movmean', smoothN, 'omitmissing');

    % Robust within-trial dF/F z-score prevents high-amplitude flies from
    % dominating event comparisons.
    dffMedian = median(dffSm(valid), 'omitnan');
    dffMAD = median(abs(dffSm(valid) - dffMedian), 'omitnan');
    dffZ = (dffSm - dffMedian) ./ max(1.4826 * dffMAD, eps);

    % Pointwise monotonic associations. Negative acceleration is coded so
    % larger values mean stronger slowing.
    accel = gradient(fwdSm, dt);
    slowing = -accel;
    rhoFwd = corr(dffSm(valid), fwdSm(valid), 'Type', 'Spearman');
    rhoAbsYaw = corr(dffSm(valid), abs(yawSm(valid)), 'Type', 'Spearman');
    rhoSlowing = corr(dffSm(valid), slowing(valid), 'Type', 'Spearman');

    % Movement offsets: sustained movement (>1 mm/s for 0.5 s) crossing below
    % 0.5 mm/s, with at least 1 s between events. Compare the second before
    % offset with the second after offset.
    moving = false(n,1);
    moving(1) = fwdSm(1) >= 1;
    for j = 2:n
        if ~moving(j-1) && fwdSm(j) >= 1
            moving(j) = true;
        elseif moving(j-1) && fwdSm(j) > 0.5
            moving(j) = true;
        end
    end
    offsets = find(moving(1:end-1) & ~moving(2:end)) + 1;
    offsets = keep_refractory(offsets, round(fs));
    preN = round(fs); postN = round(fs);
    offsets = offsets(offsets > preN & offsets + postN <= n);
    offsets = offsets(arrayfun(@(i) mean(moving(i-round(0.5*fs):i-1)) >= 0.8, offsets));
    offsetPre = event_mean(dffZ, offsets, -preN, -1);
    offsetPost = event_mean(dffZ, offsets, 0, postN);

    % Rapid slowdowns that need not end in a full stop: lowest 5% of
    % acceleration while the fly was moving >1 mm/s immediately beforehand.
    decelThreshold = prctile(accel(valid & fwdSm > 1), 5);
    slowdownEvents = find(accel <= decelThreshold & fwdSm > 0.5);
    slowdownEvents = keep_refractory(slowdownEvents, round(fs));
    slowdownEvents = slowdownEvents(slowdownEvents > preN & slowdownEvents + postN <= n);
    slowdownDff = event_mean(dffZ, slowdownEvents, -round(0.5*fs), round(0.5*fs));

    % Pivot turns: top 5% absolute yaw while translational speed is low
    % (<=2 mm/s). Refractory period prevents one turn from being counted many
    % times. Compare event dF/F against the trial's median (= robust z of 0).
    pivotThreshold = prctile(abs(yawSm(valid & fwdSm <= 2)), 95);
    pivotEvents = find(abs(yawSm) >= pivotThreshold & fwdSm <= 2);
    pivotEvents = keep_refractory(pivotEvents, round(fs));
    pivotEvents = pivotEvents(pivotEvents > preN & pivotEvents + postN <= n);
    pivotDff = event_mean(dffZ, pivotEvents, -round(0.5*fs), round(0.5*fs));

    [~, trialName] = fileparts(trialFolder);
    rows(k,:) = {trialName, n, rhoFwd, rhoAbsYaw, rhoSlowing, ...
        numel(offsets), mean(offsetPre,'omitnan'), mean(offsetPost,'omitnan'), ...
        numel(slowdownEvents), mean(slowdownDff,'omitnan'), ...
        numel(pivotEvents), mean(pivotDff,'omitnan'), ...
        decelThreshold, pivotThreshold};
end

summary = cell2table(rows, 'VariableNames', {
    'trial','nSamples','rhoDffFwd','rhoDffAbsYaw','rhoDffSlowing', ...
    'nOffsets','offsetPreDffZ','offsetPostDffZ','nSlowdowns', ...
    'slowdownDffZ','nPivots','pivotDffZ','decelThreshold','pivotYawThreshold'});

outFile = '/tmp/ves041_spontaneous_correlation_audit.csv';
writetable(summary, outFile);
disp(summary);
fprintf('Saved %s\n', outFile);

function kept = keep_refractory(events, refractoryN)
    if isempty(events), kept = events; return; end
    kept = events(1);
    for j = 2:numel(events)
        if events(j) - kept(end) >= refractoryN
            kept(end+1,1) = events(j); %#ok<AGROW>
        end
    end
end

function values = event_mean(signal, events, firstOffset, lastOffset)
    values = nan(numel(events),1);
    for j = 1:numel(events)
        idx = events(j) + (firstOffset:lastOffset);
        % Median preserves every finite dF/F sample while limiting the
        % leverage of isolated fluorescence spikes.
        values(j) = median(signal(idx), 'omitnan');
    end
end
