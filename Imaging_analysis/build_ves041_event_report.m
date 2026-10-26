function status = build_ves041_event_report(baseFolder, outputFolder, opts)
%BUILD_VES041_EVENT_REPORT Extract good-behavior calcium peaks and slowdowns.
% Run with MATLAB R2023b and Signal Processing / Statistics toolboxes.
% Inputs are read only. All derived files are written inside outputFolder.
% SourceReport reuses saved events, retaining all heat trials and only the
% no-heat trials whose mean plotted peak-aligned forward curve exceeds the
% specified threshold. Mirroring is applied to event paths before averaging.
arguments
    baseFolder (1,1) string = "/Users/sophiarenauld/stacks"
    outputFolder (1,1) string = fullfile(baseFolder,"ves041_peak_slowdown_analysis_2026-10-02")
    opts.SourceReport (1,1) string = ""
    opts.MinimumPeakAlignedForwardSpeed (1,1) double {mustBeNonnegative,mustBeFinite} = 2
    opts.MirrorLateralDisplacement (1,1) logical = true
    opts.NormalizeWithinFly (1,1) logical = false
end
if opts.NormalizeWithinFly
    assert(opts.SourceReport~="",'Normalization requires a saved, selected source report.');
    status=normalize_saved_report(opts.SourceReport,outputFolder);
    return
end
if opts.SourceReport~=""
    status=collapse_saved_report(opts.SourceReport,outputFolder, ...
        opts.MinimumPeakAlignedForwardSpeed,opts.MirrorLateralDisplacement);
    return
end
p = struct('window',[-2 2], 'outputHz',20, 'calciumSmoothSec',.25, ...
    'baselineSec',30, 'noiseMultiplier',4, 'minimumAmplitude',1, ...
    'mergeGapSec',.5, 'minimumEventSec',.25, 'kinematicSmoothSec',.25, ...
    'decelerationPercentile',10, 'minimumSpeedBefore',1, 'minimumSpeedDrop',1, ...
    'speedWindowSec',.5, 'slowdownSeparationSec',.75, ...
    'turnWindowSec',[.25 .75], 'turnThresholdDeg',15, 'maximumGapFactor',1.5);
if ~isfolder(outputFolder), mkdir(outputFolder); end
copyfile([mfilename('fullpath') '.m'],fullfile(outputFolder,'build_ves041_event_report.m'));
write_text(fullfile(outputFolder,'parameters.json'),jsonencode(p,PrettyPrint=true));
inventory = find_trials(baseFolder);
writetable(inventory,fullfile(outputFolder,'finder_tag_inventory.csv'));
selected = inventory(inventory.goodBehavior,:);
assert(~isempty(selected),'No folders currently tagged good behavior were found.');
records = cell(height(selected),1); results = cell(height(selected),1);
allPeaks = table(); allSlowdowns = table();
for k=1:height(selected)
    name=selected.trial(k); fprintf('[%d/%d] %s\n',k,height(selected),name);
    r=struct('trial',name,'condition',string('no heat'),'status',string('unavailable'), ...
        'detail',string(''),'durationSec',NaN,'imagingHz',NaN,'behaviorHz',NaN, ...
        'nPeaksDetected',NaN,'nPeaksIncluded',NaN,'nPeaksExcludedWindow',NaN, ...
        'nSlowdownCandidates',NaN,'nSlowdownsIncluded',NaN, ...
        'nSlowdownsWithReorientation',NaN,'nSlowdownsWithoutReorientation',NaN, ...
        'nSlowdownsExcludedWindow',NaN,'nSlowdownsExcludedSpeed',NaN, ...
        'nSlowdownsExcludedSeparation',NaN,'calciumThreshold',NaN, ...
        'noiseSigma',NaN,'decelerationThreshold',NaN,'upstreamVelocityWindowSec',NaN);
    if contains(lower(name),'heat'), r.condition="heat"; end
    try
        [d,meta]=load_trial(fullfile(baseFolder,name),r.condition);
        d=prepare_trial(d,p);
        [peakIndex,peakInfo]=detect_peaks(d,p);
        [peakEvents,peakTraces,nBadPeak]=extract_events(d,peakIndex,name,r.condition,"peak",p);
        [slowIndex,slowInfo]=detect_slowdowns(d,p);
        [slowEvents,slowTraces,nBadSlow]=extract_events(d,slowIndex,name,r.condition,"slowdown",p);
        r.status="included";r.detail="";r.durationSec=d.t(end)-d.t(1);
        r.imagingHz=d.fs;r.behaviorHz=meta.behaviorHz;
        r.nPeaksDetected=numel(peakIndex);r.nPeaksIncluded=height(peakEvents);
        r.nPeaksExcludedWindow=nBadPeak;r.nSlowdownCandidates=slowInfo.candidates;
        r.nSlowdownsIncluded=height(slowEvents);
        r.nSlowdownsWithReorientation=sum(slowEvents.reorientation);
        r.nSlowdownsWithoutReorientation=sum(~slowEvents.reorientation);
        r.nSlowdownsExcludedWindow=slowInfo.badWindow+nBadSlow;
        r.nSlowdownsExcludedSpeed=slowInfo.badSpeed;
        r.nSlowdownsExcludedSeparation=slowInfo.tooClose;
        r.calciumThreshold=peakInfo.threshold;r.noiseSigma=peakInfo.noiseSigma;
        r.decelerationThreshold=slowInfo.threshold;
        r.upstreamVelocityWindowSec=meta.velocityWindowSec;
        trialFolder=fullfile(outputFolder,'trials',name);
        if ~isfolder(trialFolder),mkdir(trialFolder);end
        writetable(peakEvents,fullfile(trialFolder,'peak_events.csv'));
        writetable(slowEvents,fullfile(trialFolder,'slowdown_events.csv'));
        write_text(fullfile(trialFolder,'source_metadata.json'),jsonencode(meta,PrettyPrint=true));
        save(fullfile(trialFolder,'extracted_events.mat'),'d','meta','p', ...
            'peakIndex','peakInfo','peakEvents','peakTraces','slowInfo','slowEvents','slowTraces','-v7');
        write_trace_table(peakEvents,peakTraces,fullfile(trialFolder,'peak_aligned_traces.csv'));
        write_trace_table(slowEvents,slowTraces,fullfile(trialFolder,'slowdown_aligned_traces.csv'));
        R=struct('trial',name,'condition',r.condition,'peaks',peakTraces, ...
            'slowdowns',slowTraces,'peakEvents',peakEvents,'slowEvents',slowEvents);
        plot_trial_peaks(R,trialFolder);plot_trial_slowdowns(R,trialFolder);
        plot_detection(d,peakIndex,peakEvents,slowEvents,peakInfo,name,trialFolder);
        results{k}=R;
        allPeaks=[allPeaks;peakEvents];allSlowdowns=[allSlowdowns;slowEvents]; %#ok<AGROW>
        fprintf('  included: %d peaks; %d slowdowns (%d with reorientation)\n', ...
            r.nPeaksIncluded,r.nSlowdownsIncluded,r.nSlowdownsWithReorientation);
    catch ME
        r.status="unavailable";r.detail=string(ME.message);
        fprintf(2,'  UNAVAILABLE: %s\n',getReport(ME,'extended','hyperlinks','off'));
    end
    records{k}=r;
    status=struct2table(vertcat(records{1:k}));
    writetable(status,fullfile(outputFolder,'trial_status.csv'));
end
results=vertcat(results{~cellfun(@isempty,results)});
assert(~isempty(results),'No trials were successfully analyzed.');
writetable(allPeaks,fullfile(outputFolder,'all_peak_events.csv'));
writetable(allSlowdowns,fullfile(outputFolder,'all_slowdown_events.csv'));
plot_combined_peaks(results,outputFolder);
plot_combined_slowdowns(results,outputFolder);
save(fullfile(outputFolder,'combined_results.mat'),'results','status','p','-v7');
write_status(status,outputFolder);
write_methods(status,p,outputFolder);
validate_report(results,status,p,outputFolder);
fprintf('Complete: %s\n',outputFolder);
end

function status=collapse_saved_report(source,out,threshold,mirror)
assert(string(java.io.File(char(source)).getCanonicalPath())~= ...
    string(java.io.File(char(out)).getCanonicalPath()),'Use a new output folder to preserve the source report.');
S=load(fullfile(source,'combined_results.mat'),'results','status','p');
status=S.status;status.sourceStatus=status.status;
status.meanPeakAlignedForwardMmSec=nan(height(status),1);
p=S.p;p.sourceReport=source;p.minimumNoHeatPeakAlignedForwardSpeed=threshold;
p.selectionDefinition="Retain all heat trials. For no-heat trials, retain only mean(mean(peaks.fwd,2)) > threshold over the 81 plotted time samples from -2 to +2 seconds.";
p.mirrorLateralDisplacement=mirror;
p.mirrorDefinition="Multiply the whole event x path by -1 if x(+2 s) is negative; leave y and all other traces unchanged.";
colors=turbo(numel(S.results));keep=false(numel(S.results),1);
for k=1:numel(S.results)
    R=S.results(k);row=find(status.trial==R.trial);
    assert(numel(row)==1 && status.sourceStatus(row)=="included",'Source result/status mismatch.');
    mu=mean(mean(R.peaks.fwd,2));status.meanPeakAlignedForwardMmSec(row)=mu;
    keep(k)=R.condition=="heat" || (isfinite(mu) && mu>threshold);
    if R.condition=="heat"
        status.detail(row)="Heat trial retained; cutoff applies only to no-heat trials.";
    elseif keep(k)
        status.detail(row)=sprintf('Mean plotted peak-aligned forward curve %.6f mm/s > %.6f mm/s.',mu,threshold);
    else
        status.status(row)="excluded_speed";
        status.detail(row)=sprintf('Mean plotted peak-aligned forward curve %.6f mm/s does not exceed %.6f mm/s.',mu,threshold);
    end
    S.results(k).plotColor=colors(k,:);
    fprintf('%s | %.6f mm/s | %s\n',R.trial,mu,status.status(row));
end
assert(any(keep),'No trials retained.');results=S.results(keep);
if ~isfolder(out),mkdir(out);end
copyfile([mfilename('fullpath') '.m'],fullfile(out,'build_ves041_event_report.m'));
copyfile(fullfile(source,'finder_tag_inventory.csv'),fullfile(out,'source_finder_tag_inventory.csv'));
copyfile(fullfile(source,'methods.txt'),fullfile(out,'source_methods.txt'));
writetable(status,fullfile(out,'trial_status.csv'));
writetable(status(:,{'trial','condition','sourceStatus','status','meanPeakAlignedForwardMmSec','detail'}),fullfile(out,'trial_speed_selection.csv'));
write_text(fullfile(out,'parameters.json'),jsonencode(p,PrettyPrint=true));
allPeaks=table();allSlowdowns=table();
for k=1:numel(results)
    R=results(k);X=load(fullfile(source,'trials',R.trial,'extracted_events.mat'));
    if mirror
        [R.peaks,R.peakEvents]=mirror_paths(R.peaks,R.peakEvents);
        [R.slowdowns,R.slowEvents]=mirror_paths(R.slowdowns,R.slowEvents);
    end
    results(k)=R;X.peakTraces=R.peaks;X.peakEvents=R.peakEvents;
    X.slowTraces=R.slowdowns;X.slowEvents=R.slowEvents;X.p=p;
    trialFolder=fullfile(out,'trials',R.trial);if ~isfolder(trialFolder),mkdir(trialFolder);end
    save(fullfile(trialFolder,'extracted_events.mat'),'-struct','X','-v7');
    copyfile(fullfile(source,'trials',R.trial,'source_metadata.json'),fullfile(trialFolder,'source_metadata.json'));
    writetable(R.peakEvents,fullfile(trialFolder,'peak_events.csv'));
    writetable(R.slowEvents,fullfile(trialFolder,'slowdown_events.csv'));
    write_trace_table(R.peakEvents,R.peaks,fullfile(trialFolder,'peak_aligned_traces.csv'));
    write_trace_table(R.slowEvents,R.slowdowns,fullfile(trialFolder,'slowdown_aligned_traces.csv'));
    plot_trial_peaks(R,trialFolder);
    for stem=["slowdowns_with_without_reorientation","event_detection_overview"]
        for ext=[".png",".fig"]
            copyfile(fullfile(source,'trials',R.trial,stem+ext),fullfile(trialFolder,stem+ext));
        end
    end
    allPeaks=[allPeaks;R.peakEvents];allSlowdowns=[allSlowdowns;R.slowEvents]; %#ok<AGROW>
end
writetable(allPeaks,fullfile(out,'all_peak_events.csv'));
writetable(allSlowdowns,fullfile(out,'all_slowdown_events.csv'));
plot_combined_peaks(results,out,p);plot_combined_slowdowns(results,out,p);
save(fullfile(out,'combined_results.mat'),'results','status','p','-v7');
write_collapsed_notes(status,p,out);
validate_report(results,status,p,out);
assert(all(status.meanPeakAlignedForwardMmSec(status.status=="included" & status.condition=="no heat")>threshold),'Speed selection failed.');
assert(all(ismember([S.results([S.results.condition]=="heat").trial],[results.trial])),'A heat trial was removed.');
for k=1:numel(results)
    original=S.results(find([S.results.trial]==results(k).trial,1));
    for type=["peaks","slowdowns"]
        C=results(k).(type);O=original.(type);
        fields=fieldnames(O);fields(strcmp(fields,'x'))=[];
        for j=1:numel(fields),assert(isequaln(C.(fields{j}),O.(fields{j})),'A non-lateral trace changed.');end
        if mirror
            assert(isequaln(C.xUnmirrored,O.x),'Unmirrored path differs from source.');
            if type=="peaks",E=results(k).peakEvents;else,E=results(k).slowEvents;end
            assert(isequaln(C.x,O.x.*E.lateralMirrorSign'),'Reflection mismatch.');
            assert(all(C.x(end,:)>=0),'Mirrored endpoint is negative.');
            assert(all(abs(hypot(diff(O.x),diff(O.y))-hypot(diff(C.x),diff(C.y)))<1e-10,'all'),'Reflection changed path length.');
        end
    end
end
fid=fopen(fullfile(out,'validation.txt'),'a');cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'PASS: no-heat trials selected using the plotted peak-aligned curve mean; all heat trials retained.\n');
fprintf(fid,'PASS: retained events and non-lateral traces equal source values; whole-path reflections preserve y and path lengths.\n');
fprintf('Complete: %s\n',out);
end

function [C,E]=mirror_paths(C,E)
direction=ones(size(C.x,2),1);direction(C.x(end,:)<0)=-1;
C.xUnmirrored=C.x;C.x=C.x.*direction';E.lateralMirrorSign=direction;
end

function status=normalize_saved_report(source,out)
assert(string(java.io.File(char(source)).getCanonicalPath())~= ...
    string(java.io.File(char(out)).getCanonicalPath()),'Use a separate normalized output folder.');
S=load(fullfile(source,'combined_results.mat'),'results','status','p');
assert(~is_normalized(S.results),'Source must contain physical-unit traces.');
originalResults=S.results;results=S.results;status=S.status;p=S.p;
p.normalizationSourceReport=source;p.normalization="within-fly z-score; origin-preserving RMS for geometry";
p.normalizationReference="All finite native imaging samples from the retained full recordings, pooled within fly; sample SD (N-1).";
p.geometryReference="All retained peak and slowdown windows pooled within fly. Heading: RMS heading change. Paths: RMS radial displacement, shared by x and y.";
n=numel(results);trial=[results.trial]';flyId=strings(n,1);identityBasis=strings(n,1);X=cell(n,1);
for k=1:n
    token=regexp(char(trial(k)),'^(\d{8})-.*_fly(\d+)_','tokens','once');
    if ~isempty(token)
        flyId(k)=string(token{1})+"_fly"+string(token{2});
        identityBasis(k)="Date and explicit fly number in trial folder name.";
    else
        assert(any(trial(k)==["20260217-2_ves041_vt_clbbar1","20260217-3_ves041_vt_clbbar2_actuallydark","20260217-4_ves041_vt_clbbar3"]), ...
            'No fly mapping established for %s.',trial(k));
        flyId(k)="20260217_fly1_inferred";
        identityBasis(k)="Inferred shared fly: sequential clbbar1-3 recordings at 09:58, 10:12, 10:19; later same-day folders explicitly identify fly2. No unique fly ID in acquisition CSV.";
    end
    X{k}=load(fullfile(source,'trials',trial(k),'extracted_events.mat'));
end
flyMap=table(trial,flyId,[results.condition]',identityBasis, ...
    'VariableNames',{'trial','flyId','condition','identityBasis'});
flies=unique(flyId,'stable');factors=table();
signalNames=["dffRaw","dff","localDff","fwdRaw","fwd","sideRaw","side","yawRaw","yaw","absYaw"];
for q=1:numel(flies)
    ix=find(flyId==flies(q));
    for field=signalNames
        values=[];
        for k=ix'
            v=normalization_reference(X{k}.d,field);values=[values;v(isfinite(v))]; %#ok<AGROW>
        end
        center=mean(values);scale=std(values,0);assert(isfinite(scale)&&scale>0,'Degenerate normalization reference: %s %s.',flies(q),field);
        z=(values-center)/scale;
        assert(abs(mean(z))<1e-10 && abs(std(z,0)-1)<1e-10,'Reference is not centered with unit SD.');
        one=table(flies(q),field,center,scale,numel(values),"z-score", ...
            'VariableNames',{'flyId','signal','center','scale','nReferenceSamples','method'});
        factors=[factors;one]; %#ok<AGROW>
    end
    angles=[];rSquared=[];
    for k=ix'
        for type=["peaks","slowdowns"]
            C=results(k).(type);angles=[angles;C.headingChange(:)]; %#ok<AGROW>
            rSquared=[rSquared;C.x(:).^2+C.y(:).^2]; %#ok<AGROW>
        end
    end
    headingScale=sqrt(mean(angles.^2));pathScale=sqrt(mean(rSquared));
    assert(all(isfinite([headingScale,pathScale]))&&min(headingScale,pathScale)>0,'Degenerate geometric normalization reference.');
    for field=["headingChange","absHeadingChange","x","y","xUnmirrored"]
        if startsWith(field,"x")||field=="y",scale=pathScale;method="radial RMS";else,scale=headingScale;method="heading RMS";end
        one=table(flies(q),field,0,scale,numel(angles),method, ...
            'VariableNames',{'flyId','signal','center','scale','nReferenceSamples','method'});
        factors=[factors;one]; %#ok<AGROW>
    end
end
if ~isfolder(out),mkdir(out);end
copyfile([mfilename('fullpath') '.m'],fullfile(out,'build_ves041_event_report.m'));
copyfile(fullfile(source,'methods.txt'),fullfile(out,'source_selection_methods.txt'));
copyfile(fullfile(source,'source_methods.txt'),fullfile(out,'source_detection_methods.txt'));
for file=["trial_speed_selection.csv","source_finder_tag_inventory.csv"]
    copyfile(fullfile(source,file),fullfile(out,file));
end
writetable(flyMap,fullfile(out,'fly_trial_mapping.csv'));
writetable(factors,fullfile(out,'fly_normalization_parameters.csv'));
write_text(fullfile(out,'parameters.json'),jsonencode(p,PrettyPrint=true));
allPeaks=table();allSlowdowns=table();
for k=1:n
    F=factors(factors.flyId==flyId(k),:);results(k).flyId=flyId(k);results(k).normalizedByFly=true;
    for type=["peaks","slowdowns"]
        C=results(k).(type);
        for j=1:height(F)
            field=F.signal(j);if ~isfield(C,field),continue,end
            C.(field)=(C.(field)-F.center(j))/F.scale(j);
            assert(all(isfinite(C.(field)),'all'),'Nonfinite normalized trace.');
        end
        assert(all(C.x(41,:)==0 & C.y(41,:)==0),'Normalization displaced the origin.');
        assert(all(C.x(end,:)>=0),'Normalization reversed mirrored endpoints.');
        results(k).(type)=C;
    end
    R=results(k);Q=X{k};Q.p=p;Q.normalization=F;Q.flyId=flyId(k);
    Q.peakTracesPhysical=Q.peakTraces;Q.slowTracesPhysical=Q.slowTraces;
    Q.peakTraces=R.peaks;Q.slowTraces=R.slowdowns;
    % Keep source d and event measurements in physical units for reproducibility.
    D=Q.d;
    for field=signalNames
        row=F.signal==field;v=normalization_reference(D,field);
        if field=="localDff",target="signal";else,target=field;end
        Q.normalizedData.(target)=(v-F.center(row))/F.scale(row);
    end
    Q.normalizedData.t=D.t;Q.normalizedData.normalizedByFly=true;
    trialFolder=fullfile(out,'trials',R.trial);if ~isfolder(trialFolder),mkdir(trialFolder);end
    save(fullfile(trialFolder,'extracted_events.mat'),'-struct','Q','-v7');
    copyfile(fullfile(source,'trials',R.trial,'source_metadata.json'),fullfile(trialFolder,'source_metadata.json'));
    writetable(R.peakEvents,fullfile(trialFolder,'peak_events.csv'));
    writetable(R.slowEvents,fullfile(trialFolder,'slowdown_events.csv'));
    write_trace_table(R.peakEvents,R.peaks,fullfile(trialFolder,'peak_aligned_traces.csv'));
    write_trace_table(R.slowEvents,R.slowdowns,fullfile(trialFolder,'slowdown_aligned_traces.csv'));
    plot_trial_peaks(R,trialFolder);plot_trial_slowdowns(R,trialFolder);
    E=R.peakEvents;row=F.signal=="dff";
    E.smoothedDffAtEvent=(E.smoothedDffAtEvent-F.center(row))/F.scale(row);
    info=Q.peakInfo;row=F.signal=="localDff";info.threshold=(info.threshold-F.center(row))/F.scale(row);
    plot_detection(Q.normalizedData,Q.peakIndex,E,R.slowEvents,info,R.trial,trialFolder);
    allPeaks=[allPeaks;R.peakEvents];allSlowdowns=[allSlowdowns;R.slowEvents]; %#ok<AGROW>
    fprintf('[%d/%d] %s | %s | %d peaks, %d slowdowns\n',k,n,R.trial,flyId(k),height(R.peakEvents),height(R.slowEvents));
end
writetable(allPeaks,fullfile(out,'all_peak_events.csv'));writetable(allSlowdowns,fullfile(out,'all_slowdown_events.csv'));
status.flyId=strings(height(status),1);
for k=1:n,status.flyId(status.trial==trial(k))=flyId(k);end
writetable(status,fullfile(out,'trial_status.csv'));write_status(status,out);
fid=fopen(fullfile(out,'included_trials_and_status.txt'),'a');
fprintf(fid,'\nNormalization uses %d fly groups. See fly_trial_mapping.csv for identities and the inferred February 17 group.\n',numel(flies));
fprintf(fid,'Cohort, event counts, thresholds and event-table measurements remain in original physical units. Aligned signal files and figures use fly-normalized values.\n');fclose(fid);
plot_combined_peaks(results,out,p);plot_combined_slowdowns(results,out,p);
save(fullfile(out,'combined_results.mat'),'results','originalResults','status','p','flyMap','factors','-v7');
write_normalization_methods(p,out,numel(flies));
for k=1:n
    assert(isequaln(results(k).peakEvents,originalResults(k).peakEvents)&&isequaln(results(k).slowEvents,originalResults(k).slowEvents),'Event membership or measurements changed.');
end
write_text(fullfile(out,'validation.txt'),sprintf('PASS: %d trials, %d fly groups, %d peaks and %d slowdowns retained exactly.\nPASS: full-recording references have mean zero and sample SD one after z-scoring.\nPASS: all event traces finite, mirrored endpoints nonnegative and path origins unchanged.\nPASS: event times, measurements and reorientation groups unchanged.\n',n,numel(flies),height(allPeaks),height(allSlowdowns)));
fprintf('Complete: %s\n',out);
end

function v=normalization_reference(d,field)
if field=="localDff",v=d.signal;elseif field=="absYaw",v=abs(d.yaw);else,v=d.(field);end
v=v(:);
end

function write_normalization_methods(p,out,nFlies)
lines=["VES041 fly-normalized analysis";""; ...
    "Source: "+p.normalizationSourceReport; ...
    "Same 8 no-heat trials above the original 2 mm/s peak-aligned forward-speed cutoff and all 3 heat trials. Trial selection precedes normalization. All 352 peaks and 993 slowdowns retain their original times, physical measurements and reorientation classifications."; ...
    string(nFlies)+" fly groups. Explicit fly numbers are paired with recording date, so fly3 on different days is not pooled. February 17 trials 2-4 are treated as one fly, inferred from their consecutive clbbar1-3 names and acquisition times; their acquisition metadata does not give a unique fly ID. This mapping is recorded in fly_trial_mapping.csv."; ...
    "Activity and velocity: z = (value - fly mean) / fly sample SD, with sample SD denominator N-1. Estimate each fly's mean and SD from finite native imaging samples pooled over all of its retained full recordings, including stationary periods and heat-off periods. Each native sample has equal weight; no additional interpolation is used to estimate these parameters. Excluded trials do not contribute. The same constants apply to every event and both event types from that fly."; ...
    "Normalize dffRaw, dff, localDff, fwdRaw, fwd, sideRaw, side, yawRaw and yaw separately using their corresponding full-recording signals. localDff is the existing calcium signal after subtraction of the 30 s moving-median baseline. For absYaw, calculate the fly mean and SD of abs(smoothed yaw), then z-score each event's already extracted absolute yaw. Negative normalized absolute-yaw values therefore mean below that fly's usual absolute yaw, not rotation in the opposite direction."; ...
    "Aligned geometric quantities keep their zero origin. Heading scale = sqrt(mean(headingChange^2)) across all 81 samples of all retained peak and slowdown windows for that fly. Divide both signed and absolute heading change by this common RMS angle without subtracting a mean."; ...
    "Path scale = sqrt(mean(x^2 + y^2)) across those same event-window samples. Divide mirrored x, y and xUnmirrored by that one positive distance scale. Using the same scale for both axes preserves path shape and heading; the time-zero origin and the whole-path left/right reflection are preserved. Heading and path scales use event windows, so overlapping windows contribute more than once. These panels show fly RMS units, not z-scores. RMS means root mean square."; ...
    "Time stays in seconds; heat remains a 0/1 indicator. Normalized activity, velocity, heading and displacement axes are dimensionless. fly_normalization_parameters.csv records the exact center, divisor, sample count and method for every fly/signal. There is no clipping, min-max scaling or event-specific baseline adjustment."; ...
    "Average normalized events within each trial first. Colored curves show trial means; black curves give equal weight to each contributing trial. Repeated trials from a fly still contribute separately to the overall mean, matching the source report. Normalization is within fly; aggregation remains across trials. Figures share limits for the same signal across comparison panels."; ...
    "No new smoothing, event detection, exclusion or reorientation threshold is introduced. The source 0.25 s Gaussian calcium smoothing, 0.25 s moving-median kinematic smoothing, 30 s local calcium baseline and upstream velocity processing remain in effect. See source_detection_methods.txt for all original detection, filtering, gap and source-specific details, and source_selection_methods.txt for the speed cutoff and mirroring."; ...
    "Root/per-trial event CSVs and trial_status.csv retain source physical measurements and thresholds. Aligned-trace CSVs, result trace structs, figures and mean-curve CSVs contain the normalized values. In each extracted_events.mat, d, peakEvents and slowEvents remain physical; normalizedData stores full normalized signals; peakTraces/slowTraces are normalized and peakTracesPhysical/slowTracesPhysical preserve their physical versions. combined_results.mat also includes originalResults. This report can be checked without the October 2 report."; ...
    "Per-trial detection overviews show the source detections at their original times on normalized full-recording traces. Raw and smoothed calcium each use their own saved normalization constants. The local calcium threshold is transformed using the localDff center and scale."; ...
    "Rerun: build_ves041_event_report('/Users/sophiarenauld/stacks','/path/to/new/output',SourceReport='"+p.normalizationSourceReport+"',NormalizeWithinFly=true);"];
write_text(fullfile(out,'methods.txt'),join(lines,newline));
end

function write_collapsed_notes(T,p,out)
ok=T.status=="included";
fid=fopen(fullfile(out,'included_trials_and_status.txt'),'w');cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'VES041 collapsed analysis: no-heat mean peak-aligned forward curve > %.6g mm/s\nGenerated %s\n\n',p.minimumNoHeatPeakAlignedForwardSpeed,char(datetime('now')));
fprintf(fid,'Source report: %s\nOriginal tagged cohort: %d\nRetained: %d (%d heat, %d no heat)\nNo-heat trials excluded by speed: %d\nPreviously unavailable: %d\n', ...
    p.sourceReport,height(T),sum(ok),sum(ok&T.condition=="heat"),sum(ok&T.condition=="no heat"),sum(T.status=="excluded_speed"),sum(T.status=="unavailable"));
fprintf(fid,'Retained peaks: %d\nRetained slowdowns: %d (%d with reorientation, %d without)\n\n', ...
    sum(T.nPeaksIncluded(ok)),sum(T.nSlowdownsIncluded(ok)),sum(T.nSlowdownsWithReorientation(ok)),sum(T.nSlowdownsWithoutReorientation(ok)));
fprintf(fid,'Trial | condition | plotted mean forward (mm/s) | status | peaks | slowdowns | with reorientation | without\n');
for i=1:height(T)
    fprintf(fid,'%s | %s | %.6f | %s | %g | %g | %g | %g\n',T.trial(i),T.condition(i),T.meanPeakAlignedForwardMmSec(i),T.status(i),T.nPeaksIncluded(i),T.nSlowdownsIncluded(i),T.nSlowdownsWithReorientation(i),T.nSlowdownsWithoutReorientation(i));
end
fprintf(fid,'\nExcluded rows show their original event counts, not retained events. NaN means unavailable. All heat trials are retained regardless of speed.\n');
clear cleaner
lines=["VES041 collapsed analysis: methods";""; ...
    "Source report: "+p.sourceReport; ...
    "The exact source cohort and saved event detections are reused. Current Finder tags are not reread. Original processing and detection methods are preserved in source_methods.txt."; ...
    "Selection applies only to no-heat trials. For each trial, first average forward velocity across all included calcium-peak events at each of the 81 plotted time samples from -2 to +2 seconds. Then take the arithmetic mean of that trial-mean curve across those 81 samples. Retain the trial only when this value is strictly greater than "+string(p.minimumNoHeatPeakAlignedForwardSpeed)+" mm/s. This is the curve in the right-hand forward-velocity panel of all_trials_peaks_heat_no_heat, not the full-recording mean or a mean restricted to walking."; ...
    "All heat trials are retained. The same selected trial set is used for both the peak and slowdown combined figures. All previously accepted events from each retained trial remain; no new within-trial event exclusion is introduced. The root event tables and combined_results.mat contain only retained trials. trial_status.csv and trial_speed_selection.csv also document excluded and unavailable source trials."; ...
    "Trial colors are preserved from the source plots. Thin colored curves are trial means; black curves are equal-weight averages of the contributing trial means. Event count does not determine a trial's weight. Each trial with zero events in a subgroup is omitted from that subgroup mean."; ...
    "Lateral mirroring enabled: "+string(p.mirrorLateralDisplacement)+". Each event path is already centered at (0,0) and oriented along +y at time zero. If its x displacement at +2 seconds is negative, multiply its entire x trajectory by -1; otherwise use +1. Leave y unchanged. An exactly zero endpoint is not reflected. This is a whole-path reflection rather than a pointwise absolute value, so a path may still cross x=0 earlier in the window."; ...
    "Apply reflections to individual events before averaging within and across trials. Both peak and slowdown paths are reflected in the saved data. xUnmirrored preserves the original trajectory, and lateralMirrorSign in each event table records +1 or -1. The mean aligned path in the figures uses the reflected coordinates."; ...
    "Calcium activity, forward velocity, side velocity, signed yaw, heading, event times, and reorientation classifications retain their source values. Mirroring affects only path coordinates. No new smoothing, baseline subtraction, amplitude normalization, or significance test is added. The prior odor-channel correction remains in effect."; ...
    "Per-trial peak/path figures are regenerated for retained trials. Their slowdown figures and detection overviews are copied unchanged because those signals and event lists have not changed. The original report and all its per-trial analyses remain available in the source folder."; ...
    "Rerun: build_ves041_event_report('/Users/sophiarenauld/stacks', '/path/to/new/output', SourceReport='"+p.sourceReport+"', MinimumPeakAlignedForwardSpeed="+string(p.minimumNoHeatPeakAlignedForwardSpeed)+", MirrorLateralDisplacement="+string(p.mirrorLateralDisplacement)+");"];
write_text(fullfile(out,'methods.txt'),join(lines,newline));
end

function T=find_trials(base)
entries=dir(fullfile(base,'*ves041*'));entries=entries([entries.isdir]);
rows=cell(0,3);
for i=1:numel(entries)
    % Hex output avoids the extra newline xattr adds to binary plist output.
    path=fullfile(entries(i).folder,entries(i).name);
    cmd=sprintf('/usr/bin/xattr -px com.apple.metadata:_kMDItemUserTags %s 2>/dev/null | /usr/bin/xxd -r -p | /usr/bin/plutil -convert json -o - - 2>/dev/null',shellquote(path));
    [code,out]=system(cmd);tags=strings(0,1);
    if code==0
        tags=string(jsondecode(out));tags=regexprep(tags,'\n[0-9]+$','');
    end
    rows(end+1,:)={string(entries(i).name),any(tags=="good behavior"),join(tags,'; ')}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'trial','goodBehavior','finderTags'});
T=sortrows(T,'trial');
end

function s=shellquote(s)
s=['''' strrep(char(s),'''','''"''"''') ''''];
end

function [d,m]=load_trial(folder,condition)
f=dir(fullfile(folder,'*_s_.mat'));
assert(numel(f)==1,'Expected one processed *_s_.mat with an extracted ROI; found %d.',numel(f));
source=fullfile(f.folder,f.name);a=load(source,'dq','roi');
assert(numel(a.roi)==1 && isvector(a.roi.ts),'Expected exactly one VES041 ROI trace.');
m=struct('imagingSource',source,'sourceBytes',f.bytes,'sourceModified',f.date, ...
    'roiOptions',a.roi.opt,'velocityWindowSec',NaN,'behaviorHz',NaN);
d.dffRaw=double(a.roi.ts(:));
if condition=="heat"
    f=dir(fullfile(folder,'*virmenData.mat'));
    assert(numel(f)==1,'Expected one Virmen data file; found %d.',numel(f));
    m.behaviorSource=fullfile(f.folder,f.name);v=load(m.behaviorSource,'vmData','daqVolumeTime');
    d.t=double(v.daqVolumeTime.daqVolumeTime{1,1}(:));
    vm=v.vmData;tk=column(vm,'Seconds');if isduration(tk),tk=seconds(tk);end
    fw=column(vm,'vel_for');sd=column(vm,'vel_side');yw=rad2deg(column(vm,'vel_yaw'));
    hd=unwrap(column(vm,'heading'));heat=column(vm,'heat');
    m.yawGain=vm.yawGain(1);m.translationGain=vm.transGain(1);
    m.kinematicSource='Virmen vel_for, vel_side, vel_yaw (rad/s converted to deg/s), heading';
    m.velocityWindowSec=NaN;
    assert(numel(d.t)==numel(d.dffRaw),'Imaging timestamp and ROI lengths differ.');
    validate_time(tk);valid=d.t>=tk(1) & d.t<=tk(end);
    m.imagingSamplesOutsideBehavior=sum(~valid);
    d.t=d.t(valid);d.dffRaw=d.dffRaw(valid);
    d.fwdRaw=safe_interp(tk,fw,d.t);d.sideRaw=safe_interp(tk,sd,d.t);
    d.yawRaw=safe_interp(tk,yw,d.t);d.heading=safe_interp(tk,hd,d.t);
    d.heat=interp1(tk,double(heat>1e-4),d.t,'nearest',NaN);
    m.behaviorHz=1/median(diff(tk));
else
    assert(isfield(a,'dq') && ~isempty(a.dq),'No imaging-aligned behavioral stream dq(1).');
    q=a.dq(1);d.t=double(q.t(:));d.fwdRaw=double(q.bvf(:));
    d.sideRaw=double(q.bvs(:));d.yawRaw=rad2deg(double(q.bvy(:)));
    d.heading=unwrap(double(q.bh(:)));d.heat=zeros(size(d.t));
    m.behaviorSource=source;m.kinematicSource='imaging-aligned dq(1): bvf, bvs (mm/s), bvy (rad/s), bh (rad)';
    % Odor recordings repurpose bh/bvy for the odor pulse. See the odor
    % section of master_imaging_analysis.m for the documented channel map.
    if contains(lower(folder),'odor')
        assert(isfield(q,'vh') && isfield(q,'vvy'),'Odor trial is missing vh/vvy heading channels.');
        d.heading=unwrap(double(q.vh(:)));d.yawRaw=rad2deg(double(q.vvy(:)));
        m.kinematicSource='odor-specific dq(1): bvf, bvs (mm/s), vvy (rad/s), vh (rad); bh/bvy carry odor, not heading';
    end
    m.dqOptions=q.opt;m.velocityWindowSec=q.opt.dvlensec;
    m.behaviorHz=1/median(diff(d.t));m.imagingSamplesOutsideBehavior=0;
end
validate_time(d.t);
for field=["dffRaw","fwdRaw","sideRaw","yawRaw","heading","heat"]
    assert(numel(d.(field))==numel(d.t),'Mismatched length for %s.',field);
end
m.imagingHz=1/median(diff(d.t));
m.nonfiniteSamples=struct('dff',sum(~isfinite(d.dffRaw)), ...
    'forward',sum(~isfinite(d.fwdRaw)),'side',sum(~isfinite(d.sideRaw)), ...
    'yaw',sum(~isfinite(d.yawRaw)),'heading',sum(~isfinite(d.heading)));
end

function y=column(T,f)
y=T.(f);if iscell(y),assert(numel(y)==1,'Expected one trial row.');y=y{1};end;y=y(:);
end

function validate_time(t)
assert(numel(t)>20 && all(isfinite(t)) && all(diff(t)>0),'Time must be finite, strictly increasing, and contain >20 samples.');
end

function y=safe_interp(t,x,tq)
assert(numel(t)==numel(x),'Behavior timestamps and values have different lengths.');
y=interp1(t,double(x),tq,'linear',NaN);
gap=find(diff(t)>1.5*median(diff(t)));
for k=gap(:)',y(tq>t(k) & tq<t(k+1))=NaN;end
end

function d=prepare_trial(d,p)
d.fs=1/median(diff(d.t));
d.calciumSmoothN=max(3,round(p.calciumSmoothSec*d.fs));
d.kinematicSmoothN=max(3,round(p.kinematicSmoothSec*d.fs));
d.baselineN=max(3,round(p.baselineSec*d.fs));
d.dff=segment_smooth(d.t,d.dffRaw,'gaussian',d.calciumSmoothN,p);
d.baseline=segment_smooth(d.t,d.dff,'movmedian',d.baselineN,p);
d.signal=d.dff-d.baseline;
for f=["fwd","side","yaw"]
    d.(f)=segment_smooth(d.t,d.(f+"Raw"),'movmedian',d.kinematicSmoothN,p);
end
d.accel=nan(size(d.t));blocks=segments(d.t,isfinite(d.fwd),p.maximumGapFactor);
for j=1:size(blocks,1)
    ix=blocks(j,1):blocks(j,2);if numel(ix)>=3,d.accel(ix)=gradient(d.fwd(ix),d.t(ix));end
end
end

function y=segment_smooth(t,x,method,n,p)
y=nan(size(x));blocks=segments(t,isfinite(x),p.maximumGapFactor);
for j=1:size(blocks,1)
    ix=blocks(j,1):blocks(j,2);y(ix)=smoothdata(x(ix),method,n);
end
end

function B=segments(t,valid,gapFactor)
cut=[true;diff(t)>gapFactor*median(diff(t))];
first=find(valid & (cut | [true;~valid(1:end-1)]));
last=find(valid & ([~valid(2:end);true] | [cut(2:end);true]));
B=[first last];
end

function [idx,info]=detect_peaks(d,p)
delta=diff(d.signal);valid=isfinite(delta) & diff(d.t)<=p.maximumGapFactor/d.fs;
sigma=median(abs(delta(valid)-median(delta(valid))))/(sqrt(2)*.6745);
threshold=max(p.minimumAmplitude,p.noiseMultiplier*sigma);
info=struct('noiseSigma',sigma,'threshold',threshold);idx=zeros(0,1);
B=segments(d.t,isfinite(d.signal),p.maximumGapFactor);
for j=1:size(B,1)
    ix=(B(j,1):B(j,2))';sig=d.signal(ix);if numel(ix)<3,continue,end
    [~,loc]=findpeaks(sig,'MinPeakHeight',threshold,'MinPeakProminence',threshold);
    active=sig>=threshold;change=diff([false;active;false]);a=find(change==1);b=find(change==-1)-1;
    if isempty(a),continue,end
    starts=a(1);stops=[];stop=b(1);
    for k=2:numel(a)
        if a(k)-stop-1<=round(p.mergeGapSec*d.fs),stop=b(k);
        else,stops(end+1,1)=stop;starts(end+1,1)=a(k);stop=b(k);end %#ok<AGROW>
    end
    stops(end+1,1)=stop;
    for k=1:numel(starts)
        if stops(k)-starts(k)+1<max(1,round(p.minimumEventSec*d.fs)),continue,end
        c=loc(loc>=starts(k) & loc<=stops(k));if isempty(c),continue,end
        [~,best]=max(sig(c));idx(end+1,1)=ix(c(best)); %#ok<AGROW>
    end
end
idx=sort(idx);
end

function [idx,info]=detect_slowdowns(d,p)
moving=d.accel(isfinite(d.accel) & d.fwd>=p.minimumSpeedBefore);
threshold=NaN;if ~isempty(moving),threshold=prctile(moving,p.decelerationPercentile);end
info=struct('threshold',threshold,'candidates',0,'badWindow',0,'badSpeed',0,'tooClose',0);
idx=zeros(0,1);if ~isfinite(threshold)||threshold>=0,return,end
active=isfinite(d.accel) & d.accel<=threshold & d.fwd>=0;
B=segments(d.t,active,p.maximumGapFactor);info.candidates=size(B,1);
for j=1:size(B,1)
    ix=B(j,1):B(j,2);[~,m]=min(d.accel(ix));i=ix(m);t0=d.t(i);
    if ~window_ok(d,t0,p),info.badWindow=info.badWindow+1;continue,end
    pre=d.fwd(d.t>=t0-p.speedWindowSec & d.t<t0);
    post=d.fwd(d.t>=t0 & d.t<=t0+p.speedWindowSec);
    if median(pre)<p.minimumSpeedBefore || median(pre)-median(post)<p.minimumSpeedDrop
        info.badSpeed=info.badSpeed+1;continue
    end
    idx(end+1,1)=i; %#ok<AGROW>
end
% Keep the strongest deceleration first; rejected speed/window events cannot
% suppress an otherwise eligible slowdown.
[~,order]=sort(d.accel(idx));kept=zeros(0,1);
for j=order(:)'
    if isempty(kept)||all(abs(d.t(idx(j))-d.t(kept))>=p.slowdownSeparationSec)
        kept(end+1,1)=idx(j); %#ok<AGROW>
    else,info.tooClose=info.tooClose+1;end
end
idx=sort(kept);
end

function ok=window_ok(d,t0,p)
start=find(d.t<=t0+p.window(1),1,'last');stop=find(d.t>=t0+p.window(2),1,'first');
ok=~isempty(start) && ~isempty(stop);if ~ok,return,end
ix=start:stop;
ok=all(diff(d.t(ix))<=p.maximumGapFactor/d.fs) && ...
    all(all(isfinite([d.dffRaw(ix) d.dff(ix) d.fwdRaw(ix) d.fwd(ix) ...
    d.sideRaw(ix) d.side(ix) d.yawRaw(ix) d.yaw(ix) d.heading(ix)])));
end

function [E,C,nBad]=extract_events(d,indices,trial,condition,type,p)
C=struct();C.time=(p.window(1):1/p.outputHz:p.window(2))';
fields=["dffRaw","dff","localDff","fwdRaw","fwd","sideRaw","side", ...
    "yawRaw","yaw","absYaw","headingChange","absHeadingChange","x","y","heat"];
for f=fields,C.(f)=zeros(numel(C.time),0);end
rows=cell(0,15);nBad=0;
for i=indices(:)'
    t0=d.t(i);if ~window_ok(d,t0,p),nBad=nBad+1;continue,end
    tq=t0+C.time;one=struct();
    for f=["dffRaw","dff","fwdRaw","fwd","sideRaw","side","yawRaw","yaw"]
        one.(f)=interp1(d.t,d.(f),tq,'linear');
    end
    one.localDff=interp1(d.t,d.signal,tq,'linear');
    heading=interp1(d.t,d.heading,tq,'linear');relative=heading-d.heading(i);
    one.headingChange=rad2deg(relative);one.absHeadingChange=abs(one.headingChange);
    one.absYaw=abs(one.yaw);one.heat=interp1(d.t,d.heat,tq,'nearest');
    % Body-referenced fictive paths: origin and heading aligned at the event.
    one.x=cumtrapz(C.time,one.fwdRaw.*sin(relative)+one.sideRaw.*cos(relative));
    one.y=cumtrapz(C.time,one.fwdRaw.*cos(relative)-one.sideRaw.*sin(relative));
    z=round(-p.window(1)*p.outputHz)+1;one.x=one.x-one.x(z);one.y=one.y-one.y(z);
    for f=fields,C.(f)(:,end+1)=one.(f);end
    h=interp1(d.t,d.heading,t0+p.turnWindowSec,'linear');
    turn=rad2deg(h(2)-h(1));reorient=abs(turn)>=p.turnThresholdDeg;
    pre=median(d.fwd(d.t>=t0-p.speedWindowSec & d.t<t0));
    post=median(d.fwd(d.t>=t0 & d.t<=t0+p.speedWindowSec));
    rows(end+1,:)={trial,condition,type,height_from_rows(rows)+1,t0,i, ...
        reorient,turn,pre,post,pre-post,d.dff(i),d.signal(i),d.heat(i),mean(one.heat)}; %#ok<AGROW>
end
names={'trial','condition','eventType','eventId','eventTimeSec', ...
    'imagingIndex','reorientation','netHeadingChangeDeg','speedBefore','speedAfter', ...
    'speedDrop','smoothedDffAtEvent','localDffAtEvent','heatOnAtEvent','heatFractionInWindow'};
if isempty(rows)
    types=[repmat({'string'},1,3) repmat({'double'},1,3) {'logical'} repmat({'double'},1,8)];
    E=table('Size',[0 numel(names)],'VariableTypes',types,'VariableNames',names);
else,E=cell2table(rows,'VariableNames',names);end
for f=["trial","condition","eventType"],E.(f)=string(E.(f));end
E.reorientation=logical(E.reorientation);
end

function n=height_from_rows(rows),n=size(rows,1);end

function write_trace_table(E,C,path)
n=height(E);nt=numel(C.time);T=table(repelem((1:n)',nt),repmat(C.time,n,1), ...
    'VariableNames',{'eventId','relativeTimeSec'});
fields=fieldnames(C);fields(strcmp(fields,'time'))=[];
for k=1:numel(fields),T.(fields{k})=reshape(C.(fields{k}),[],1);end
writetable(T,path);
end

function plot_trial_peaks(R,out)
fig=new_figure([1480 890]);tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
names=["dff","fwd","absYaw","side","headingChange"];
labels=["VES041 dF/F","Forward (mm/s)","|Yaw| (deg/s)","Side (mm/s)","Heading from t=0 (deg)"];
labels=normalized_labels(R,names,labels);
for j=1:5
    ax=nexttile(tl);trace_plot(ax,R.peaks.time,R.peaks.(names(j)),[.15 .4 .67],labels(j));
    xlabel(ax,'Time from calcium peak (s)');
end
ax=nexttile(tl);path_plot(ax,R.peaks.x,R.peaks.y,[.15 .4 .67]);
title(ax,'Flat paths aligned at peak');
sgtitle(tl,sprintf('%s | %s | %d calcium peaks',R.trial,R.condition,height(R.peakEvents)), ...
    'Interpreter','none','FontSize',16);
if isfield(R.peaks,'xUnmirrored')
    xlabel(ax,'Mirrored lateral displacement (mm)');
    subtitle(tl,'Thin lines: individual events. Thick line: event mean. Each whole path is reflected to x(+2 s) >= 0 before averaging; y is unchanged.');
else
    subtitle(tl,'Thin lines: individual events. Thick line: event mean. Path origin and heading aligned at t=0; left/right turns are not reflected.');
end
if is_normalized(R)
    xlabel(ax,'Mirrored lateral / fly path RMS');ylabel(ax,'Forward / fly path RMS');
    subtitle(tl,'Activity and velocity: within-fly z-scores. Heading and paths: divided by fly RMS. Thin: events; thick: event mean.');
end
save_figure(fig,out,'peaks_kinematics_flat_paths');
end

function plot_trial_slowdowns(R,out)
fig=new_figure([1340 1080]);tl=tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
names=["dff","fwd","absYaw","absHeadingChange"];
labels=["VES041 dF/F","Forward (mm/s)","|Yaw| (deg/s)","|Heading from t=0| (deg)"];
labels=normalized_labels(R,names,labels);
for c=1:2
    choose=R.slowEvents.reorientation==(c==1);color=group_color(R.condition,c==1);
    for j=1:4
        ax=nexttile(tl,(j-1)*2+c);trace_plot(ax,R.slowdowns.time,R.slowdowns.(names(j))(:,choose),color,labels(j));
        if j==1,title(ax,sprintf('%s (n=%d)',turn_label(c==1),sum(choose)));end
        if j==4,xlabel(ax,'Time from strongest deceleration (s)');end
    end
end
sgtitle(tl,R.trial+" | "+R.condition,'Interpreter','none','FontSize',16);
subtitle(tl,'Thin lines: individual slowdowns. Thick line: event mean. Reorientation: |net heading change| >=15 deg from +0.25 to +0.75 s.');
if is_normalized(R)
    subtitle(tl,'Activity and velocity: within-fly z-scores. Heading: divided by fly RMS. Thin: events; thick: event mean. Groups use the original 15 deg criterion.');
    match_signal_limits(fig,labels);
end
save_figure(fig,out,'slowdowns_with_without_reorientation');
end

function plot_combined_peaks(R,out,p)
if nargin<3,p=struct();end
fig=new_figure([1470 1400]);tl=tiledlayout(fig,5,2,'TileSpacing','compact','Padding','compact');
names=["dff","fwd","absYaw","side"];
labels=["VES041 dF/F","Forward (mm/s)","|Yaw| (deg/s)","Side (mm/s)"];
labels=normalized_labels(R,names,labels);
allColors=report_colors(R);summary=table();
for c=1:2
    conditions=["heat","no heat"];ix=find([R.condition]==conditions(c));
    ix=ix(arrayfun(@(i)height(R(i).peakEvents)>0,ix));n=sum(arrayfun(@(i)height(R(i).peakEvents),ix));
    for j=1:4
        ax=nexttile(tl,(j-1)*2+c);hold(ax,'on');M=[];
        for k=ix
            avg=mean(R(k).peaks.(names(j)),2);M(:,end+1)=avg; %#ok<AGROW>
            plot(ax,R(k).peaks.time,avg,'Color',allColors(k,:),'LineWidth',1,'DisplayName',short_name(R(k).trial));
        end
        if ~isempty(M)
            avg=mean(M,2);plot(ax,R(1).peaks.time,avg,'k','LineWidth',2.6,'DisplayName','Equal-trial mean');
            S=table(repmat(conditions(c),numel(avg),1),repmat(names(j),numel(avg),1),R(1).peaks.time,avg, ...
                repmat(numel(ix),numel(avg),1),'VariableNames',{'condition','signal','timeSec','equalTrialMean','nTrials'});
            summary=[summary;S]; %#ok<AGROW>
        else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
        style_axis(ax);xline(ax,0,':','HandleVisibility','off');xlim(ax,[-2 2]);ylabel(ax,labels(j));
        if j==1,title(ax,sprintf('%s | %d trials | %d peaks',conditions(c),numel(ix),n));end
        if j==4,xlabel(ax,'Time from calcium peak (s)');end
    end
    ax=nexttile(tl,8+c);hold(ax,'on');X=[];Y=[];
    for k=ix
        mx=mean(R(k).peaks.x,2);my=mean(R(k).peaks.y,2);X(:,end+1)=mx;Y(:,end+1)=my; %#ok<AGROW>
        plot(ax,mx,my,'Color',allColors(k,:),'LineWidth',1,'DisplayName',short_name(R(k).trial));
    end
    if ~isempty(X)
        plot(ax,mean(X,2),mean(Y,2),'k','LineWidth',2.6,'DisplayName','Equal-trial mean');
        pathMeans=[mean(X,2) mean(Y,2)];pathNames=["x","y"];
        for j=1:2
            S=table(repmat(conditions(c),size(pathMeans,1),1),repmat(pathNames(j),size(pathMeans,1),1), ...
                R(1).peaks.time,pathMeans(:,j),repmat(numel(ix),size(pathMeans,1),1), ...
                'VariableNames',{'condition','signal','timeSec','equalTrialMean','nTrials'});
            summary=[summary;S]; %#ok<AGROW>
        end
    end
    plot(ax,0,0,'ko','MarkerFaceColor','w','HandleVisibility','off');
    style_axis(ax);axis(ax,'equal');xlabel(ax,'Lateral displacement (mm)');ylabel(ax,'Forward displacement (mm)');
    title(ax,'Mean aligned flat paths');
    if isfield(R(1).peaks,'xUnmirrored'),xlabel(ax,'Mirrored lateral displacement (mm)');title(ax,'Mean aligned flat paths (mirrored)');end
    if is_normalized(R),xlabel(ax,'Mirrored lateral / fly path RMS');ylabel(ax,'Forward / fly path RMS');end
    if ~isempty(ix),legend(ax,'Location','eastoutside','FontSize',8,'Interpreter','none');end
end
match_signal_limits(fig,labels);
sgtitle(tl,'VES041 calcium peaks across good-behavior trials','FontSize',18);
subtitle(tl,'Color: one mean per trial. Black: equal weight per contributing trial. Heat labels describe trial type, including heat-off periods within heat trials.');
if isfield(p,'minimumNoHeatPeakAlignedForwardSpeed')
    subtitle(tl,sprintf('No heat: plotted peak-aligned forward mean > %g mm/s. All heat trials retained. Paths mirrored before averaging. Color: trial mean; black: equal-trial mean.',p.minimumNoHeatPeakAlignedForwardSpeed));
end
if is_normalized(R)
    sgtitle(tl,'VES041 calcium peaks | normalized within fly','FontSize',18);
    subtitle(tl,'Activity and velocity: fly z-scores; paths: divided by fly RMS. Original 8 no-heat + 3 heat trials. Color: trial mean; black: equal-trial mean.');
    pathAxes=findall(fig,'Type','axes');
    pathAxes=pathAxes(arrayfun(@(a)string(a.YLabel.String)=="Forward / fly path RMS",pathAxes));
    xl=vertcat(pathAxes.XLim);yl=vertcat(pathAxes.YLim);
    set(pathAxes,'XLim',[min(xl(:,1)) max(xl(:,2))],'YLim',[min(yl(:,1)) max(yl(:,2))]);
end
save_figure(fig,out,'all_trials_peaks_heat_no_heat');
writetable(summary,fullfile(out,'peak_equal_trial_means.csv'));
end

function plot_combined_slowdowns(R,out,p)
if nargin<3,p=struct();end
fig=new_figure([1550 1450]);tl=tiledlayout(fig,8,2,'TileSpacing','compact','Padding','compact');
names=["dff","fwd","absYaw","absHeadingChange"];
labels=["VES041 dF/F","Forward (mm/s)","|Yaw| (deg/s)","|Heading from t=0| (deg)"];
labels=normalized_labels(R,names,labels);
colors=report_colors(R);summary=table();
for group=1:2
    conditions=["heat","no heat"];ix=find([R.condition]==conditions(group));
    for c=1:2
    reorient=c==1;
        for j=1:4
            ax=nexttile(tl,((group-1)*4+j-1)*2+c);hold(ax,'on');M=[];count=0;nt=0;
            for k=ix
                choose=R(k).slowEvents.reorientation==reorient;
                if ~any(choose),continue,end
                nt=nt+1;count=count+sum(choose);
                avg=mean(R(k).slowdowns.(names(j))(:,choose),2);M(:,end+1)=avg; %#ok<AGROW>
                plot(ax,R(k).slowdowns.time,avg,'Color',colors(k,:),'LineWidth',1);
            end
            if ~isempty(M)
                avg=mean(M,2);plot(ax,R(1).slowdowns.time,avg,'k','LineWidth',2.6);
                S=table(repmat(conditions(group),numel(avg),1),repmat(reorient,numel(avg),1), ...
                    repmat(names(j),numel(avg),1),R(1).slowdowns.time,avg,repmat(nt,numel(avg),1), ...
                    'VariableNames',{'condition','reorientation','signal','timeSec','equalTrialMean','nTrials'});
                summary=[summary;S]; %#ok<AGROW>
            else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
            style_axis(ax);xline(ax,0,':');xlim(ax,[-2 2]);ylabel(ax,labels(j));
            if j==1,title(ax,sprintf('%s | %s | %d trials, %d slowdowns',conditions(group),turn_label(reorient),nt,count));end
            if j==4,xlabel(ax,'Time from strongest deceleration (s)');end
        end
    end
end
% Matching limits for each signal make the four groups directly comparable.
match_signal_limits(fig,labels);
sgtitle(tl,'VES041 activity and kinematics at slowdowns','FontSize',18);
subtitle(tl,'Four groups: heat/no heat x with/without reorientation. Color: trial means. Black: equal-trial mean. Heading criterion: >=15 deg from +0.25 to +0.75 s.');
if isfield(p,'minimumNoHeatPeakAlignedForwardSpeed')
    subtitle(tl,sprintf('Same trial selection as peak figure: no-heat peak-aligned forward mean > %g mm/s; all heat trials retained. Color: trial mean; black: equal-trial mean.',p.minimumNoHeatPeakAlignedForwardSpeed));
end
if is_normalized(R)
    sgtitle(tl,'VES041 slowdowns | normalized within fly','FontSize',18);
    subtitle(tl,'Activity and velocity: fly z-scores; heading: divided by fly RMS. Original events and turn groups. Color: trial mean; black: equal-trial mean.');
end
save_figure(fig,out,'all_trials_slowdowns_four_groups');
writetable(summary,fullfile(out,'slowdown_equal_trial_means.csv'));
end

function colors=report_colors(R)
if isfield(R,'plotColor'),colors=vertcat(R.plotColor);else,colors=turbo(numel(R));end
end

function yes=is_normalized(R)
yes=isfield(R,'normalizedByFly') && all([R.normalizedByFly]);
end

function labels=normalized_labels(R,names,labels)
if ~is_normalized(R),return,end
keys=["dff","fwd","absYaw","side","headingChange","absHeadingChange"];
values=["VES041 dF/F (fly z-score)","Forward (fly z-score)","|Yaw| (fly z-score)", ...
    "Side (fly z-score)","Heading / fly RMS","|Heading| / fly RMS"];
for j=1:numel(names),labels(j)=values(keys==names(j));end
end

function plot_detection(d,idx,E,S,info,name,out)
fig=new_figure([1500 850]);tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);plot(ax,d.t,d.dffRaw,'Color',[.75 .75 .75]);hold(ax,'on');
plot(ax,d.t,d.dff,'Color',[.1 .25 .5]);plot(ax,d.t(idx),d.dff(idx),'o','Color',[.7 .5 .1],'MarkerSize',5);
plot(ax,E.eventTimeSec,E.smoothedDffAtEvent,'o','Color',[.7 .05 .15],'MarkerFaceColor',[.7 .05 .15],'MarkerSize',4);
ylabel(ax,'VES041 dF/F');legend(ax,{'Input','Smoothed','All detected','Included'},'Location','eastoutside');style_axis(ax);
if is_normalized(d),ylabel(ax,'VES041 dF/F (fly z-score)');end
ax=nexttile(tl);plot(ax,d.t,d.signal,'Color',[.1 .25 .5]);hold(ax,'on');yline(ax,info.threshold,'r--');
ylabel(ax,'dF/F above local baseline');style_axis(ax);
if is_normalized(d),ylabel(ax,'Local dF/F (fly z-score)');end
ax=nexttile(tl);plot(ax,d.t,d.fwd,'k');hold(ax,'on');
for c=1:2
    choose=S.reorientation==(c==1);v=interp1(d.t,d.fwd,S.eventTimeSec(choose));
    plot(ax,S.eventTimeSec(choose),v,'o','Color',group_color("no heat",c==1),'MarkerFaceColor',group_color("no heat",c==1),'MarkerSize',4);
end
ylabel(ax,'Forward (mm/s)');xlabel(ax,'Trial time (s)');style_axis(ax);
if is_normalized(d),ylabel(ax,'Forward (fly z-score)');end
legend(ax,{'Forward speed','Slowdown with reorientation','Slowdown without reorientation'},'Location','eastoutside');
sgtitle(tl,name+" | event detection overview",'Interpreter','none');
if is_normalized(d),subtitle(tl,'Within-fly normalization of source detections; event times and membership unchanged.');end
save_figure(fig,out,'event_detection_overview');
end

function trace_plot(ax,t,M,color,label)
hold(ax,'on');
if isempty(M),text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');
else
    plot(ax,t,M,'Color',.78+.22*color,'LineWidth',.5);
    plot(ax,t,mean(M,2),'Color',color,'LineWidth',2.5);
end
xline(ax,0,':','Color',[.2 .2 .2]);xlim(ax,[-2 2]);ylabel(ax,label);style_axis(ax);
end

function path_plot(ax,X,Y,color)
hold(ax,'on');
if isempty(X),text(ax,.5,.5,'No eligible paths','Units','normalized','HorizontalAlignment','center');
else
    plot(ax,X,Y,'Color',.78+.22*color,'LineWidth',.5);
    plot(ax,mean(X,2),mean(Y,2),'Color',color,'LineWidth',2.5);
    plot(ax,mean(X(1,:)),mean(Y(1,:)),'s','Color',color,'MarkerFaceColor','w');
    plot(ax,mean(X(end,:)),mean(Y(end,:)),'^','Color',color,'MarkerFaceColor',color);
end
plot(ax,0,0,'ko','MarkerFaceColor','w');axis(ax,'equal');style_axis(ax);
xlabel(ax,'Lateral displacement (mm)');ylabel(ax,'Forward displacement (mm)');
end

function fig=new_figure(sz)
fig=figure('Visible','off','Color','w','Position',[40 40 sz],'DefaultAxesFontSize',11,'DefaultTextInterpreter','none');
end

function style_axis(ax)
ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;
end

function match_signal_limits(fig,labels)
axs=findall(fig,'Type','axes');
for lab=labels
    subset=axs(arrayfun(@(a)string(a.YLabel.String)==lab,axs));
    lim=vertcat(subset.YLim);if ~isempty(lim),set(subset,'YLim',[min(lim(:,1)) max(lim(:,2))]);end
end
end

function c=group_color(condition,reorient)
if condition=="heat",c=[.75 .19 .16];else,c=[.12 .37 .67];end
if ~reorient,c=.25+.65*c;end
end

function s=turn_label(yes)
if yes,s='With reorientation';else,s='Without reorientation';end
end

function s=short_name(s),s=extractBefore(s,'_ves041');end

function save_figure(fig,out,name)
savefig(fig,fullfile(out,string(name)+'.fig'));
exportgraphics(fig,fullfile(out,string(name)+'.png'),'Resolution',160);
close(fig);
end

function write_text(path,s)
fid=fopen(path,'w');assert(fid>=0,'Cannot open output %s.',path);
cleaner=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',s);
end

function write_status(T,out)
fid=fopen(fullfile(out,'included_trials_and_status.txt'),'w');cleaner=onCleanup(@()fclose(fid));
ok=T.status=="included";
fprintf(fid,'VES041 good-behavior trial report\nGenerated %s\n\n',char(datetime('now')));
fprintf(fid,'Finder-tagged trials: %d\nIncluded: %d (%d heat, %d no heat)\nUnavailable: %d\nExcluded by speed: %d\n', ...
    height(T),sum(ok),sum(ok&T.condition=="heat"),sum(ok&T.condition=="no heat"),sum(T.status=="unavailable"),sum(T.status=="excluded_speed"));
fprintf(fid,'Included peaks: %d\nIncluded slowdowns: %d (%d with reorientation, %d without)\n\n', ...
    sum(T.nPeaksIncluded(ok)),sum(T.nSlowdownsIncluded(ok)),sum(T.nSlowdownsWithReorientation(ok)),sum(T.nSlowdownsWithoutReorientation(ok)));
fprintf(fid,'Heat/no heat describes the trial type. No-heat trials include bar jumps, gratings, dark, and odor where tagged.\n');
fprintf(fid,'Reorientation: absolute net heading change >=15 degrees from +0.25 to +0.75 seconds after the slowdown.\n');
fprintf(fid,'Counts below describe events with complete finite -2 to +2 second windows. NaN in the CSV means unavailable, not zero events.\n\n');
for i=1:height(T)
    fprintf(fid,'%s\n  Condition: %s | Status: %s\n',T.trial(i),T.condition(i),T.status(i));
    if ~ok(i),fprintf(fid,'  Reason: %s\n\n',T.detail(i));continue,end
    fprintf(fid,'  Duration: %.2f s | Imaging: %.4f Hz\n',T.durationSec(i),T.imagingHz(i));
    fprintf(fid,'  Peaks: %d detected, %d included, %d excluded for incomplete/nonfinite/gapped windows\n', ...
        T.nPeaksDetected(i),T.nPeaksIncluded(i),T.nPeaksExcludedWindow(i));
    fprintf(fid,'  Slowdowns: %d candidates, %d included; with reorientation %d, without %d\n', ...
        T.nSlowdownCandidates(i),T.nSlowdownsIncluded(i),T.nSlowdownsWithReorientation(i),T.nSlowdownsWithoutReorientation(i));
    fprintf(fid,'  Slowdown exclusions: %d window, %d speed criteria, %d separation\n', ...
        T.nSlowdownsExcludedWindow(i),T.nSlowdownsExcludedSpeed(i),T.nSlowdownsExcludedSeparation(i));
    fprintf(fid,'  Calcium threshold: %.5g dF/F | noise sigma: %.5g | acceleration threshold: %.5g mm/s^2\n\n', ...
        T.calciumThreshold(i),T.noiseSigma(i),T.decelerationThreshold(i));
end
end

function write_methods(T,p,out)
% T and p are also saved verbatim in trial_status.csv / parameters.json.
lines=[
"VES041 peak and slowdown analysis: methods"
""
"Scope and trial inclusion"
"The folder inventory was read from /Users/sophiarenauld/stacks. A trial was selected only if its current macOS Finder tag exactly matched 'good behavior' (after removing the tag color suffix). Tags were read from extended attributes, because Spotlight mdls did not return them reliably. finder_tag_inventory.csv preserves the full selection record. The VES042 reference in the request was confirmed to mean VES041."
"Heat/no heat groups refer to trial names containing 'heat' versus all other selected trial names. A heat trial contributes its full usable recording, including heat-off periods. No sensory epochs were excluded: bar jumps, grating motion, odor, and heat epochs remain. heatOnAtEvent and heatFractionInWindow are saved per event."
"Selected trials without exactly one processed *_s_.mat containing an extracted ROI are listed as unavailable, with NaN counts. They are not treated as recordings with zero events."
""
"Source signals, time bases, and upstream processing"
"Input calcium is the single saved roi.ts trace from each processed *_s_.mat. This is already normalized calcium data; the report does not re-extract ROIs or recompute the original fluorescence normalization. source_metadata.json in each trial records the source file, its size and modification time, the ROI normalization options, and available behavioral options. In the local SCOPA tsnorm implementation, dff005000 first subtracts each trace's minimum, then computes (F-F0)/F0 with F0 set to the fifth percentile over the whole recording. All 16 included files store dff005000, degdtr=0, wavp=[], and channorm=0: no polynomial detrending, wavelet filtering, or channel normalization was requested in their saved ROI normalization options. These are inherited inputs, not additional operations performed by this report."
"For no-heat trials, calcium and behavior use the matching imaging-aligned dq(1) timestamps. dq(1).bvf and bvs are forward and side speed in mm/s; bvy is converted from radians/s to degrees/s; bh is the unwrapped heading in radians. The upstream SCOPA dqmake velocity derivative fits use each recording's saved dvlensec and dvord settings. These already smooth the velocities; per-trial settings are preserved rather than replaced by an assumed common value."
"The two odor trials are a documented exception to this channel map. As specified in the odor section of master_imaging_analysis.m, bh/bvy encode the odor pulse. Heading therefore comes from vh, and yaw speed from vvy (converted from rad/s to deg/s), for trials with 'odor' in their folder name. Forward and side speed remain bvf and bvs. This prevents odor-pulse edges from being interpreted as reorientations or used to reconstruct flat paths."
"For heat trials, calcium timestamps come from daqVolumeTime in the Virmen MAT file. Virmen Seconds, vel_for, vel_side, vel_yaw, and heading supply behavior. vel_yaw is converted from radians/s to degrees/s. Continuous behavior is linearly interpolated onto the actual imaging timestamps; heat state (>0.0001) is nearest-neighbor interpolated. Calcium outside the recorded behavioral time range is removed. No clock shift or rate correction is fitted. The saved Virmen heading and velocity convention, including its recorded yaw gain, is retained; this can differ from physical ball heading in the no-heat files."
"Time vectors must be finite, strictly increasing, and match their signals. Behavior interpolation across gaps greater than 1.5 times its median sampling interval yields missing values. Smoothing is performed separately within finite, contiguous segments and never across a gap. Missing values are not imputed."
""
"Calcium peaks"
"The detector adapts the strict merged-event detector in compare_ves041_peak_amplitude_pivot_all_trials.m. Detection is performed at the native imaging rate for both heat and no-heat recordings. Gaussian smoothing uses max(3, round(0.25*imagingHz)) samples. A centered 30-second moving median (rounded to samples) is subtracted from the smoothed calcium to define local amplitude. Windows shrink at recording/segment edges."
"Noise sigma is median(abs(diff(local_signal)-median(diff(local_signal)))) / (sqrt(2)*0.6745), using only finite adjacent differences without a time gap. The threshold is max(1 dF/F, 4*noise_sigma). A peak must meet both the local height and prominence threshold. Above-threshold runs are merged across interruptions of at most round(0.5*imagingHz) samples. A merged run must span at least max(1, round(0.25*imagingHz)) samples, including any merged interruption. The highest qualifying local peak represents each merged run. There is no additional minimum separation between distinct merged events."
"All detected peaks are counted before window exclusions. Peaks are included in aligned plots only if a complete finite -2 to +2 second calcium and kinematic window exists. No high-amplitude clipping, winsorization, or manual removal is applied. Smoothed absolute dF/F is shown in figures; local-baseline-subtracted dF/F and the input trace are also saved."
""
"Slowdowns"
"The slowdown rule adapts compare_ves041_slowdown_contexts.m. Forward, side, and signed yaw speed each receive a centered moving median of max(3, round(0.25*imagingHz)) samples on the imaging time base. Forward acceleration is MATLAB gradient(smoothed_forward, actual_time). The acceleration threshold is the 10th percentile among finite samples with forward speed >=1 mm/s. A nonnegative threshold yields no candidates."
"Candidate runs have acceleration at or below that negative threshold and forward speed >=0 mm/s. The most negative acceleration in each run is time zero. Median forward speed in [-0.5,0) s must be >=1 mm/s, and the median decrease to [0,+0.5] s must be >=1 mm/s. These rules select forward slowdowns, not every stop and not reverse-walking decelerations."
"Candidates must first have complete finite -2 to +2 second windows and pass the speed criteria. A 0.75-second exclusion distance is then applied in order of strongest deceleration, so an ineligible candidate cannot suppress an eligible one. Candidates removed by the window, speed, or separation rules are counted separately. Slowdowns are detected independently of calcium; they need not contain a calcium peak."
""
"Reorientation"
"The existing 15-degree net-heading criterion and +0.25 to +0.75 second evaluation interval from compare_ves041_peak_amplitude_pivot_all_trials.m are used. Unwrapped heading is interpolated at the two exact interval endpoints. A slowdown is 'with reorientation' when the absolute endpoint difference is >=15 degrees, otherwise 'without reorientation'. No additional <=5 mm/s pivot-speed cutoff is applied, since these categories describe all accepted slowdowns. No calcium criterion enters this classification."
"This is a net-heading definition within the stated interval: a turn that reverses within the interval, or occurs outside it, may fall in the without-reorientation group. Heat results use the recorded Virmen heading, whose yaw gain is recorded in source_metadata.json. These figures are descriptive comparisons, not a claim of causal reorientation."
""
"Aligned traces, flat paths, and averages"
"All included events are interpolated onto -2:0.05:+2 seconds (20 Hz, 81 samples). Interpolation puts events on a common display grid and does not add temporal resolution. Input and smoothed dF/F, forward, side, and signed yaw speed, absolute yaw, heading relative to time zero, local calcium amplitude, and heat state are saved. No event baseline subtraction or amplitude normalization is applied to the displayed calcium."
"Flat paths are reconstructed for each event from the interpolated input forward and side velocities and the unwrapped heading relative to its value at time zero. xdot = forward*sin(heading) + side*cos(heading); ydot = forward*cos(heading) - side*sin(heading). Trapezoidal time integration gives millimeters. Each path is translated to (0,0) at time zero, with its event heading along +y. Left and right turns are not reflected. The curves are fictive paths reconstructed from the saved kinematics, not arena position traces. No separate positional smoothing is added."
"Per-trial figures show every included event as a thin line and their arithmetic mean as a thick line. The path mean averages x and y at matched times; opposite turns can cancel in that mean. The start of the mean path is a square, its end a triangle, and time zero a circle."
"Combined figures show one mean per contributing trial in color, with a black arithmetic mean of those trial means. Each trial therefore has equal weight, regardless of its event count. A trial with zero events in a group does not enter that group's mean. Counts of contributing trials and events are printed on the panels. The plots show no confidence intervals or significance tests. Repeated trials from the same fly are not independent biological replicates; these averages are trial-weighted, not fly-weighted."
""
"Output files and rerunning"
"included_trials_and_status.txt and trial_status.csv give selection status, thresholds, counts, and exclusion reasons. Each trials/<trial>/ folder has PNG and editable MATLAB FIG plots, an event-detection overview, per-event CSV tables, 81-sample-per-event aligned CSV data, extracted_events.mat, and source_metadata.json. The root has the two combined figures, event tables, equal-trial mean CSVs, combined_results.mat, parameters.json, this methods file, and the MATLAB source used."
"Rerun in MATLAB: addpath('/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Imaging_analysis'); build_ves041_event_report('/Users/sophiarenauld/stacks', '/path/to/new/output_folder');"
"The current Finder tags are read on every run, so changes to those tags change the cohort. The saved inventory documents this run. No source recording or existing analysis file is changed."
];
write_text(fullfile(out,'methods.txt'),join(lines,newline));
fid=fopen(fullfile(out,'methods.txt'),'a');cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'\nRecorded processing settings by included trial\n');
fprintf(fid,'Trial | upstream derivative window (s) | upstream polynomial order | added Gaussian/median window (samples; s) | calcium normalization | yaw gain\n');
for i=find(T.status=="included")'
    meta=jsondecode(fileread(fullfile(out,'trials',T.trial(i),'source_metadata.json')));
    n=max(3,round(p.calciumSmoothSec*T.imagingHz(i)));gain=NaN;order=NaN;
    if isfield(meta,'yawGain'),gain=meta.yawGain;end
    if isfield(meta,'dqOptions'),order=meta.dqOptions.dvord;end
    fprintf(fid,'%s | %g | %g | %d; %.6f | %s | %g\n',T.trial(i),T.upstreamVelocityWindowSec(i),order,n,n/T.imagingHz(i),meta.roiOptions.nrm.nrmstr,gain);
end
fprintf(fid,'NaN derivative window/order for heat means the Virmen velocity input has no stored SCOPA derivative-fit setting; the unused processed dq stream is not substituted. NaN yaw gain for no-heat means no Virmen yaw gain applies.\n');
end

function validate_report(R,T,p,out)
ok=T.status=="included";assert(sum(ok)==numel(R),'Mismatch in included trial count.');
assert(all(T.nPeaksDetected(ok)==T.nPeaksIncluded(ok)+T.nPeaksExcludedWindow(ok)),'Peak count mismatch.');
assert(all(T.nSlowdownsIncluded(ok)==T.nSlowdownsWithReorientation(ok)+T.nSlowdownsWithoutReorientation(ok)),'Slowdown category mismatch.');
assert(all(T.nSlowdownCandidates(ok)==T.nSlowdownsIncluded(ok)+T.nSlowdownsExcludedWindow(ok)+T.nSlowdownsExcludedSpeed(ok)+T.nSlowdownsExcludedSeparation(ok)),'Slowdown exclusion count mismatch.');
for i=1:numel(R)
    for type=["peaks","slowdowns"]
        C=R(i).(type);if type=="peaks",E=R(i).peakEvents;else,E=R(i).slowEvents;end
        assert(size(C.dff,2)==height(E),'Trace/event count mismatch.');
        f=fieldnames(C);f(strcmp(f,'time'))=[];
        for j=1:numel(f),assert(all(isfinite(C.(f{j})),'all'),'Nonfinite saved event data.');end
        z=find(abs(C.time)<1e-8);assert(all(abs(C.x(z,:))<1e-10) && all(abs(C.y(z,:))<1e-10),'Path origin mismatch.');
        assert(all(E.reorientation==(abs(E.netHeadingChangeDeg)>=p.turnThresholdDeg)),'Reorientation criterion mismatch.');
        if type=="slowdowns"
            assert(all(diff(E.eventTimeSec)>=p.slowdownSeparationSec-1e-8),'Slowdown separation violation.');
            assert(all(E.speedBefore>=p.minimumSpeedBefore & E.speedDrop>=p.minimumSpeedDrop),'Slowdown speed criterion violation.');
        end
    end
end
write_text(fullfile(out,'validation.txt'),sprintf(['PASS: %d included trials. Event counts, exclusion counts, category partitions, ' ...
    'finite event data, path origins, reorientation labels, and slowdown criteria verified.\n'],numel(R)));
end
