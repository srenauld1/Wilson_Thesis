%% Test whether VES041 peaks mark low-speed reorientation in 20250812-22
% Starts from the exact dF/F peaks used by plotting_ves041. Behavioral
% association is tested against circularly shifted peak times, preserving
% the number and temporal spacing of fluorescence events while breaking
% their alignment to behavior.

clearvars
rng(544);

trialFolder = '/Users/sophiarenauld/stacks/20250812-22_ves041_fly3_bbarjumps2';
outputFolder = fullfile(trialFolder, 'smr_analysis', 'peak_reorientation');
if ~isfolder(outputFolder), mkdir(outputFolder); end

loaded = load(fullfile(trialFolder, 'smr_analysis', 'fly_data.mat'), 'fly_data');
data = loaded.fly_data;
t = data.dq(1).t(:);
dff = data.roi.ts(:);
fwd = data.dq(1).bvf(:);
yaw = data.dq(1).bvy_deg(:);
n = min([numel(t),numel(dff),numel(fwd),numel(yaw)]);
t=t(1:n); dff=dff(1:n); fwd=fwd(1:n); yaw=yaw(1:n);
fs = 1/(max(t)/length(t));

% Exact fluorescence peak definition used in plotting_ves041.m.
[peakHeight, peakIdx] = findpeaks(dff, 'MinPeakDistance',30, 'MinPeakHeight',2.5);
analysisWindowN = round(2*fs);
keep = peakIdx > analysisWindowN & peakIdx+analysisWindowN <= n;
peakIdx=peakIdx(keep); peakHeight=peakHeight(keep);

% Operational definitions. These are intentionally explicit and easy to
% alter for sensitivity checks.
eventN = round(0.5*fs);        % behavior measured +/-0.5 s around peak
speedN = round(0.5*fs);        % pre/post speed windows
pivotThreshold = 50;           % deg/s
lowSpeedThreshold = 2;         % mm/s, median around peak
speedChangeThreshold = 1;      % mm/s

metrics = measure_behavior(peakIdx,t,fwd,yaw,eventN,speedN, ...
    pivotThreshold,lowSpeedThreshold,speedChangeThreshold);

% Bar-related means the peak lies in the existing +/-2 s jump buffer.
jumpMask = false(n,1);
rawJumpMask = logical(data.omit_jump_blocks(:));
jumpMask(1:min(n,numel(rawJumpMask))) = rawJumpMask(1:min(n,numel(rawJumpMask)));
barRelated = jumpMask(peakIdx);

% Circular-shift null: shift the complete set of peak indices by one random
% offset each iteration. This preserves inter-peak timing and behavioral
% autocorrelation better than independently sampled random points.
nPerm = 10000;
minShift = round(5*fs);
allowedLength = n-2*analysisWindowN;
null = nan(nPerm,8);
for p=1:nPerm
    shift = randi([minShift,allowedLength-minShift]);
    shifted = analysisWindowN + mod((peakIdx-analysisWindowN-1)+shift,allowedLength)+1;
    m = measure_behavior(shifted,t,fwd,yaw,eventN,speedN, ...
        pivotThreshold,lowSpeedThreshold,speedChangeThreshold);
    null(p,:) = summarize_metrics(m);
end
observed = summarize_metrics(metrics);

metricNames = ["meanPeakAbsYaw","medianSpeedAtPeak","meanAbsHeadingChange", ...
    "fractionPivot","fractionLowSpeed","fractionLowSpeedPivot", ...
    "fractionSlowdown","fractionSpeedup"];

% Directional permutation p values: yaw, heading change, and fractions are
% tested for enrichment; speed is tested for unusually low values.
pValue = nan(size(observed));
for q=1:numel(observed)
    if q==2
        pValue(q)=(1+sum(null(:,q)<=observed(q)))/(nPerm+1);
    else
        pValue(q)=(1+sum(null(:,q)>=observed(q)))/(nPerm+1);
    end
end
nullLow = prctile(null,2.5,1);
nullHigh = prctile(null,97.5,1);
statsTable = table(metricNames',observed',mean(null,1)',nullLow',nullHigh',pValue', ...
    'VariableNames',{'metric','observed','nullMean','nullCI2_5','nullCI97_5','permutationP'});

peakTable = table(t(peakIdx),peakHeight,metrics.speedBefore,metrics.speedAfter, ...
    metrics.speedChange,metrics.medianSpeed,metrics.peakAbsYaw, ...
    metrics.absHeadingChange,metrics.isPivot,metrics.isLowSpeed, ...
    metrics.isLowSpeedPivot,metrics.isSlowdown,metrics.isSpeedup,barRelated, ...
    'VariableNames',{'peakTime','peakDff','speedBefore','speedAfter','speedChange', ...
    'medianSpeedAtPeak','peakAbsYaw','absoluteHeadingChange','isPivot','isLowSpeed', ...
    'isLowSpeedPivot','isSlowdown','isSpeedup','barJumpRelated'});

writetable(statsTable,fullfile(outputFolder,'peak_reorientation_permutation_stats.csv'));
writetable(peakTable,fullfile(outputFolder,'peak_reorientation_events.csv'));
disp(statsTable);
fprintf(['n=%d; pivot=%d; low speed=%d; joint=%d; slowdown=%d; ' ...
    'speedup=%d; stable=%d; bar-related=%d\n'],height(peakTable), ...
    sum(peakTable.isPivot),sum(peakTable.isLowSpeed),sum(peakTable.isLowSpeedPivot), ...
    sum(peakTable.isSlowdown),sum(peakTable.isSpeedup), ...
    sum(~peakTable.isSlowdown & ~peakTable.isSpeedup),sum(peakTable.barJumpRelated));

function m=measure_behavior(indices,t,fwd,yaw,eventN,speedN,pivotThr,lowThr,changeThr)
nEvents=numel(indices);
m.speedBefore=nan(nEvents,1); m.speedAfter=nan(nEvents,1);
m.speedChange=nan(nEvents,1); m.medianSpeed=nan(nEvents,1);
m.peakAbsYaw=nan(nEvents,1); m.absHeadingChange=nan(nEvents,1);
for k=1:nEvents
    i=indices(k);
    m.speedBefore(k)=median(fwd(i-speedN:i-1),'omitnan');
    m.speedAfter(k)=median(fwd(i:i+speedN),'omitnan');
    m.speedChange(k)=m.speedAfter(k)-m.speedBefore(k);
    win=i-eventN:i+eventN;
    m.medianSpeed(k)=median(fwd(win),'omitnan');
    m.peakAbsYaw(k)=max(abs(yaw(win)),[],'omitnan');
    % Absolute net angular displacement over the event window.
    m.absHeadingChange(k)=abs(trapz(t(win),yaw(win)));
end
m.isPivot=m.peakAbsYaw>=pivotThr;
m.isLowSpeed=m.medianSpeed<=lowThr;
m.isLowSpeedPivot=m.isPivot&m.isLowSpeed;
m.isSlowdown=m.speedChange<=-changeThr;
m.isSpeedup=m.speedChange>=changeThr;
end

function values=summarize_metrics(m)
values=[mean(m.peakAbsYaw,'omitnan'),median(m.medianSpeed,'omitnan'), ...
    mean(m.absHeadingChange,'omitnan'),mean(m.isPivot),mean(m.isLowSpeed), ...
    mean(m.isLowSpeedPivot),mean(m.isSlowdown),mean(m.isSpeedup)];
end
