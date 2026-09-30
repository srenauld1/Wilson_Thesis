function [eventTable, summaryTable] = compare_ves041_slowdown_contexts(baseFolder, outputFolder, opts)
%COMPARE_VES041_SLOWDOWN_CONTEXTS Compare VES041 activity at slowdowns.
%
% Separates slowdown events along two independent dimensions:
%   1. spontaneous versus stimulus-associated
%   2. with a pivot turn versus without a pivot turn
%
% Only folders carrying the macOS Finder tag "good behavior" are included.
% Supported stimulus contexts are heat, odor, bar jumps, and open-loop
% grating motion. Closed-loop and dark recordings have no stimulus mask.
%
% Outputs:
%   eventTable   - one row per slowdown
%   summaryTable - one row per context/pivot combination
%   CSV, PNG, and editable MATLAB FIG files in outputFolder
%
% No finite dF/F values are removed based on amplitude. Robust medians and
% Spearman-style/event-based summaries limit isolated-point leverage without
% deleting high dF/F observations.
%
% Example:
%   [events, summary] = compare_ves041_slowdown_contexts();

arguments
    baseFolder (1,1) string = "/Users/sophiarenauld/stacks"
    outputFolder (1,1) string = fullfile(baseFolder, "ves041_behavior_quantification", "slowdown_contexts")
    opts.SmoothWindowSec (1,1) double {mustBePositive} = 0.25
    opts.DecelerationPercentile (1,1) double {mustBeGreaterThan(opts.DecelerationPercentile,0),mustBeLessThan(opts.DecelerationPercentile,50)} = 10
    opts.MinimumSpeedBefore (1,1) double {mustBeNonnegative} = 1
    opts.MinimumSpeedDrop (1,1) double {mustBeNonnegative} = 1
    opts.RefractorySec (1,1) double {mustBePositive} = 0.75
    opts.EventWindowSec (1,2) double = [-2 2]
    opts.PreWindowSec (1,2) double = [-1 -0.25]
    opts.ResponseWindowSec (1,2) double = [-0.25 0.5]
    opts.SpeedWindowSec (1,1) double {mustBePositive} = 0.5
    opts.PivotWindowSec (1,2) double = [-0.5 0.5]
    opts.PivotYawThreshold (1,1) double {mustBeNonnegative} = 50
    opts.StimulusAssociationWindowSec (1,2) double = [-0.25 2]
    opts.GratingVelocityThreshold (1,1) double {mustBeNonnegative} = 0.02
    opts.HeatThreshold (1,1) double {mustBeNonnegative} = 0.0001
    opts.MakePlots (1,1) logical = true
end

if ~isfolder(outputFolder), mkdir(outputFolder); end
trialFolders = find_good_behavior_trials(baseFolder);
fprintf('Found %d folders tagged "good behavior".\n', numel(trialFolders));

allRows = cell(0,15);
skipped = strings(0,1);

for k = 1:numel(trialFolders)
    folder = trialFolders(k);
    fprintf('Processing %d/%d: %s\n', k, numel(trialFolders), folder);
    try
        trial = load_trial_context(folder, opts);
    catch ME
        warning('Skipping %s: %s', folder, ME.message);
        skipped(end+1,1) = folder + ": " + string(ME.message); %#ok<AGROW>
        continue
    end

    valid = isfinite(trial.t) & isfinite(trial.dff) & ...
        isfinite(trial.fwd) & isfinite(trial.yaw);
    if nnz(valid) < 20
        skipped(end+1,1) = folder + ": insufficient finite samples"; %#ok<AGROW>
        continue
    end

    dt = median(diff(trial.t(valid)));
    fs = 1 / dt;
    smoothN = max(3, round(opts.SmoothWindowSec * fs));
    fwd = smoothdata(trial.fwd, 'movmedian', smoothN, 'omitmissing');
    yaw = smoothdata(trial.yaw, 'movmedian', smoothN, 'omitmissing');
    dff = smoothdata(trial.dff, 'movmean', smoothN, 'omitmissing');
    accel = gradient(fwd, dt);

    movingAccel = accel(valid & fwd >= opts.MinimumSpeedBefore);
    if isempty(movingAccel), continue; end
    accelThreshold = prctile(movingAccel, opts.DecelerationPercentile);
    candidate = accel <= accelThreshold & fwd >= 0;
    events = representative_events(candidate, accel);
    events = keep_refractory(events, round(opts.RefractorySec * fs));

    preSpeedN = round(opts.SpeedWindowSec * fs);
    postSpeedN = preSpeedN;
    boundsN = round(max(abs([opts.EventWindowSec opts.PreWindowSec ...
        opts.ResponseWindowSec opts.PivotWindowSec opts.StimulusAssociationWindowSec])) * fs);
    events = events(events > max(boundsN,preSpeedN) & ...
        events + max(boundsN,postSpeedN) <= numel(fwd));

    for j = 1:numel(events)
        idx = events(j);
        speedBefore = median(fwd(idx-preSpeedN:idx-1), 'omitnan');
        speedAfter = median(fwd(idx:idx+postSpeedN), 'omitnan');
        speedDrop = speedBefore - speedAfter;
        if speedBefore < opts.MinimumSpeedBefore || speedDrop < opts.MinimumSpeedDrop
            continue
        end

        preIdx = idx + seconds_to_indices(opts.PreWindowSec, fs);
        responseIdx = idx + seconds_to_indices(opts.ResponseWindowSec, fs);
        pivotIdx = idx + seconds_to_indices(opts.PivotWindowSec, fs);
        stimIdx = idx + seconds_to_indices(opts.StimulusAssociationWindowSec, fs);
        traceIdx = idx + seconds_to_indices(opts.EventWindowSec, fs);

        dffPre = median(dff(preIdx), 'omitnan');
        dffEvent = median(dff(responseIdx), 'omitnan');
        dffDelta = dffEvent - dffPre;
        peakAbsYaw = max(abs(yaw(pivotIdx)), [], 'omitnan');
        hasPivot = peakAbsYaw >= opts.PivotYawThreshold;
        stimulusAssociated = any(trial.stimulusMask(stimIdx));

        if stimulusAssociated
            context = "stimulus-associated";
            stimulusType = trial.stimulusType;
        else
            context = "spontaneous";
            stimulusType = "none";
        end
        if hasPivot, pivotLabel = "pivot"; else, pivotLabel = "no pivot"; end
        group = context + " | " + pivotLabel;

        allRows(end+1,:) = {trial.name, trial.paradigm, context, ... %#ok<AGROW>
            stimulusType, pivotLabel, group, trial.t(idx), speedBefore, ...
            speedAfter, speedDrop, peakAbsYaw, dffPre, dffEvent, dffDelta, ...
            struct('time',trial.t(traceIdx)-trial.t(idx), ...
                    'dff',dff(traceIdx),'fwd',fwd(traceIdx),'yaw',abs(yaw(traceIdx)))};
    end
end

eventTable = cell2table(allRows, 'VariableNames', { ...
    'trial','paradigm','context','stimulusType','pivotStatus','group', ...
    'eventTime','speedBefore','speedAfter','speedDrop','peakAbsYaw', ...
    'dffPre','dffEvent','dffDelta','trace'});
eventTable.trial = string(eventTable.trial);
eventTable.paradigm = string(eventTable.paradigm);
eventTable.context = string(eventTable.context);
eventTable.stimulusType = string(eventTable.stimulusType);
eventTable.pivotStatus = string(eventTable.pivotStatus);
eventTable.group = string(eventTable.group);

if isempty(eventTable)
    error('compare_ves041_slowdown_contexts:NoEvents', ...
        'No slowdowns passed the configured criteria.');
end

eventCsv = removevars(eventTable, 'trace');
writetable(eventCsv, fullfile(outputFolder, 'ves041_slowdown_events.csv'));
summaryTable = summarize_events(eventTable);
writetable(summaryTable, fullfile(outputFolder, 'ves041_slowdown_context_summary.csv'));
if ~isempty(skipped)
    writelines(skipped, fullfile(outputFolder, 'skipped_trials.txt'));
end

if opts.MakePlots
    plot_event_traces(eventTable, outputFolder, opts.EventWindowSec);
    plot_event_effects(eventTable, outputFolder);
end
disp(summaryTable);
end

function folders = find_good_behavior_trials(baseFolder)
entries = dir(fullfile(baseFolder, '*ves041*'));
entries = entries([entries.isdir]);
folders = strings(0,1);
for k = 1:numel(entries)
    folder = string(fullfile(entries(k).folder, entries(k).name));
    command = sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"', ...
        strrep(folder, '"', '\"'));
    [status, tags] = system(command);
    if status == 0 && contains(tags, 'good behavior')
        folders(end+1,1) = folder; %#ok<AGROW>
    end
end
folders = sort(folders);
end

function trial = load_trial_context(folder, opts)
[~, name] = fileparts(folder);
trial.name = string(name);
trial.paradigm = infer_paradigm(lower(string(name)));

flyFile = fullfile(folder, 'smr_analysis', 'fly_data.mat');
if isfile(flyFile)
    loaded = load(flyFile, 'fly_data');
    data = loaded.fly_data;
    t = data.dq(2).t(:);
    dff = data.roi.ts(:);
    fwd = data.dq(2).bvf(:);
    if isfield(data.dq(2), 'bvy_deg')
        yaw = data.dq(2).bvy_deg(:);
    else
        yaw = rad2deg(data.dq(2).bvy(:));
    end
    n = min([numel(t),numel(dff),numel(fwd),numel(yaw)]);
    trial.t = t(1:n); trial.dff = dff(1:n);
    trial.fwd = fwd(1:n); trial.yaw = yaw(1:n);
    trial.stimulusMask = false(n,1);
    trial.stimulusType = "none";

    if contains(trial.paradigm, "bar jump") && isfield(data, 'jump_detected')
        jump = logical(data.jump_detected(:));
        jump = jump(1:min(n,numel(jump)));
        jumpMask = false(n,1); jumpMask(1:numel(jump)) = jump;
        trial.stimulusMask = expand_onsets(jumpMask, trial.t, [0 2]);
        trial.stimulusType = "bar jump";
    elseif trial.paradigm == "grating"
        trial.stimulusMask = grating_mask(data, trial.t, opts.GratingVelocityThreshold);
        trial.stimulusType = "grating";
    elseif trial.paradigm == "odor"
        rawT = data.dq(1).t(:);
        odor = data.dq(1).bh(:) > -0.5;
        trial.stimulusMask = interp1(rawT, double(odor), trial.t, 'nearest', 0) > 0.5;
        trial.stimulusType = "odor";
    end
    return
end

% Heat trials use Virmen time and kinematics rather than fly_data.mat.
vmFile = dir(fullfile(folder, '*virmenData.mat'));
imFile = dir(fullfile(folder, '*_or_z0z1_s_.mat'));
if numel(vmFile) == 1 && numel(imFile) == 1
    vm = load(fullfile(vmFile.folder, vmFile.name));
    imaging = load(fullfile(imFile.folder, imFile.name));
    tKin = seconds(vm.vmData.Seconds{1,1}); tKin = tKin(:);
    fwdKin = vm.vmData.vel_for{1,1}; fwdKin = fwdKin(:);
    yawKin = vm.vmData.vel_yaw{1,1}; yawKin = yawKin(:);
    heat = vm.vmData.heat; if iscell(heat), heat = cell2mat(heat); end
    heat = heat(:);
    t = vm.daqVolumeTime.daqVolumeTime{1,1}(:);
    dff = imaging.roi.ts(:);
    n = min(numel(t),numel(dff)); t=t(1:n); dff=dff(1:n);
    within = t >= min(tKin) & t <= max(tKin);
    t=t(within); dff=dff(within);
    fwd = interp1(tKin,fwdKin,t,'linear',NaN);
    yaw = interp1(tKin,yawKin,t,'linear',NaN);
    heatIm = interp1(tKin,double(heat),t,'nearest',0);
    trial.t=t; trial.dff=dff; trial.fwd=fwd; trial.yaw=yaw;
    trial.stimulusMask = heatIm > opts.HeatThreshold;
    trial.stimulusType = "heat";
    return
end

error('No usable fly_data.mat or complete heat/Virmen input set.');
end

function paradigm = infer_paradigm(name)
if contains(name,'heat'), paradigm="heat";
elseif contains(name,'odor'), paradigm="odor";
elseif contains(name,'grating'), paradigm="grating";
elseif contains(name,'jump'), paradigm="bar jump";
elseif contains(name,'dark'), paradigm="dark";
else, paradigm="closed loop";
end
end

function mask = grating_mask(data, imagingTime, threshold)
rawT = data.dq(1).t(:);
if isfield(data.dq(1), 'vvyo'), velocity = data.dq(1).vvyo(:);
else, velocity = data.dq(1).vvy(:);
end
moving = abs(velocity) >= threshold;
mask = interp1(rawT, double(moving), imagingTime, 'nearest', 0) > 0.5;
end

function expanded = expand_onsets(mask, t, windowSec)
expanded = false(size(mask));
onsets = find(mask & ~[false;mask(1:end-1)]);
for k = 1:numel(onsets)
    expanded(t >= t(onsets(k))+windowSec(1) & t <= t(onsets(k))+windowSec(2)) = true;
end
end

function events = representative_events(candidate, score)
d = diff([false;candidate(:);false]);
starts = find(d==1); stops = find(d==-1)-1;
events = zeros(numel(starts),1);
for k = 1:numel(starts)
    [~,rel] = min(score(starts(k):stops(k)));
    events(k) = starts(k)+rel-1;
end
end

function kept = keep_refractory(events, refractoryN)
if isempty(events), kept=events; return; end
kept=events(1);
for k=2:numel(events)
    if events(k)-kept(end)>=refractoryN, kept(end+1,1)=events(k); end %#ok<AGROW>
end
end

function idx = seconds_to_indices(windowSec, fs)
idx = round(windowSec(1)*fs):round(windowSec(2)*fs);
end

function summary = summarize_events(events)
groups = ["spontaneous | no pivot";"spontaneous | pivot"; ...
          "stimulus-associated | no pivot";"stimulus-associated | pivot"];
rows = cell(numel(groups),8);
for k=1:numel(groups)
    subset = events(events.group==groups(k),:);
    rows(k,:)={groups(k),height(subset),numel(unique(subset.trial)), ...
        mean(subset.dffDelta,'omitnan'),median(subset.dffDelta,'omitnan'), ...
        mean(subset.speedDrop,'omitnan'),median(subset.peakAbsYaw,'omitnan'), ...
        strjoin(unique(string(subset.stimulusType)),', ')};
end
summary=cell2table(rows,'VariableNames',{'group','nEvents','nTrials', ...
    'meanDffDelta','medianDffDelta','meanSpeedDrop','medianPeakAbsYaw','stimuli'});
end

function plot_event_traces(events, outputFolder, eventWindow)
groups = ["spontaneous | no pivot","spontaneous | pivot", ...
          "stimulus-associated | no pivot","stimulus-associated | pivot"];
colors = lines(4); fields={'dff','fwd','yaw'};
labels={'dF/F','Forward velocity (mm/s)','Absolute yaw speed (deg/s)'};
fig=figure('Visible','off','Color','w','Position',[100 100 1100 850]);
for row=1:3
    ax=subplot(3,1,row); hold(ax,'on');
    for g=1:numel(groups)
        subset=events(events.group==groups(g),:);
        if isempty(subset), continue; end
        commonT=linspace(eventWindow(1),eventWindow(2),201)';
        mat=nan(numel(commonT),height(subset));
        for j=1:height(subset)
            tr=subset.trace{j};
            mat(:,j)=interp1(tr.time,tr.(fields{row}),commonT,'linear',NaN);
        end
        plot(ax,commonT,median(mat,2,'omitnan'),'LineWidth',2,'Color',colors(g,:), ...
            'DisplayName',sprintf('%s (n=%d)',groups(g),height(subset)));
    end
    xline(ax,0,':k'); ylabel(ax,labels{row}); grid(ax,'on');
    if row==1, legend(ax,'Location','best'); title(ax,'VES041 slowdown context comparison'); end
end
xlabel('Time from maximal deceleration (s)');
savefig(fig,fullfile(outputFolder,'ves041_slowdown_context_traces.fig'));
print(fig,fullfile(outputFolder,'ves041_slowdown_context_traces.png'),'-dpng','-r200');
close(fig);
end

function plot_event_effects(events, outputFolder)
fig=figure('Visible','off','Color','w','Position',[100 100 1100 500]);
subplot(1,2,1); boxchart(categorical(events.group),events.dffDelta); hold on;
scatter(categorical(events.group),events.dffDelta,18,'filled','MarkerFaceAlpha',0.35);
yline(0,':k'); ylabel('Event dF/F - pre-event dF/F'); title('All slowdown events');
xtickangle(20); grid on;
subplot(1,2,2);
trialMeans=groupsummary(events,{'trial','group'},'mean','dffDelta');
boxchart(categorical(trialMeans.group),trialMeans.mean_dffDelta); hold on;
scatter(categorical(trialMeans.group),trialMeans.mean_dffDelta,35,'filled');
yline(0,':k'); ylabel('Mean dF/F change per trial'); title('Trials as replicates');
xtickangle(20); grid on;
savefig(fig,fullfile(outputFolder,'ves041_slowdown_context_effects.fig'));
print(fig,fullfile(outputFolder,'ves041_slowdown_context_effects.png'),'-dpng','-r200');
close(fig);
end
