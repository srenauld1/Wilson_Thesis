function report=rerun_ves041_hand_drawn_oct19_oct20()
% Full master-analysis rerun using manual ROI traces, with separate outputs.
% No source MAT or prior analysis folder is modified. Normalization, behavior
% thresholds and peak settings are unchanged; only the ROI source changes.
% Additional heat figures and matching retain rectangular-only inference.
old=get(groot,'defaultFigureVisible');
c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
here=fileparts(mfilename('fullpath'));addpath(here);
addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
sub='hand_drawn_15s_ave_f_20261026';
out=fullfile('/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis',sub);
assert(~isfolder(out),'Output exists; choose a new run name to preserve it.');mkdir(out);
heatout=fullfile(out,'heat_review');mkdir(heatout);
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat'); ...
    dir('/Users/sophiarenauld/stacks/2026-10-20/**/*_or_z0z1_s_.mat')];
assert(numel(files)==16,'Expected 16 recordings; inspect changed input inventory.');
rows=cell(0,7);audit=cell(0,13);hrows=cell(0,10);allP=table;allE=table;
rectIds=["20261019_5_1","20261019_7_1","20261019_8_2", ...
    "20261020_11_2","20261020_12_1","20261020_13_1"];
for k=1:numel(files)
    input=fullfile(files(k).folder,files(k).name);id=string(erase(files(k).name,'_or_z0z1_s_.mat'));
    dest=fullfile(files(k).folder,sub);kind="unknown";acquisition="unknown";tic;
    diary(fullfile(out,id+'.log'));fprintf('MANUAL RUN %d/%d %s\n',k,numel(files),id);
    try
        assert(~isfolder(dest),'Per-trial destination already exists.');
        a=load(input,'dq','roi');before=a.roi.ts;
        assert(strcmp(a.roi.opt.nrm.nrmstr,'dff005015'));
        a=import_ves041_protocol(a,files(k).folder);
        assert(isequaln(a.roi.ts,a.roi.pxfilter.manualts));
        kind=string(a.protocol.kind);acquisition=string(a.protocol.acquisitionStatus);
        identical=isequaln(before,a.roi.ts);
        runMaster(input,dest,here);
        % Save complete trace provenance; do not trust labels alone.
        saved=load(fullfile(dest,'fly_data.mat'),'fly_data');
        assert(isequaln(saved.fly_data.roi.ts,a.roi.pxfilter.manualts), ...
            'Saved master output does not contain the manual trace');clear saved
        a.peakTraceSource=[char(input) ' | full hand-drawn ROI'];
        P=plot_ves041_peak_heat_context(a,dest);
        P=addvars(P,repmat(id,height(P),1),'Before',1,'NewVariableNames','recording');
        allP=[allP;P]; %#ok<AGROW>
        t=a.dq(1).t(:);d=double(a.roi.ts(:));b=double(before(:));
        [nold,sold]=countPeaks(t,b);[nnew,snew]=countPeaks(t,d);
        f=a.dq(1).bvf(:);r=abs(rad2deg(a.dq(1).bvy(:)));
        audit(end+1,:)={id,kind,identical,max(abs(b-d)),corr(b,d,'Rows','complete'), ...
            nold,nnew,sold,snew,corr(b,f,'Rows','complete'),corr(d,f,'Rows','complete'), ...
            corr(b,r,'Rows','complete'),corr(d,r,'Rows','complete')}; %#ok<AGROW>
        if ismember(kind,["heat-barrier","heat-stripes"])
            E=plot_ves041_fly1_heat_exploration(files(k).folder,fullfile(heatout,id),files(k),sub);
            allE=[allE;E]; %#ok<AGROW>
            nbar=NaN;
            if kind=="heat-barrier",p=a.protocol;nbar=numel(unique(p.trial_number(p.barrier_active>0&isfinite(p.trial_number)&p.trial_number>0)));end
            hrows(end+1,:)={id,kind,acquisition,"success",height(P),sum(P.eligibleAlignedWindow),height(E),nbar,string(dest),""}; %#ok<AGROW>
            repair_saved_figure_visibility(fullfile(dest,'heat_exploration'));
        end
        repair_saved_figure_visibility(dest);
        status="success";err="";
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));disp(err);
    end
    elapsed=toc;close all force;diary off;
    rows(end+1,:)={id,string(input),string(dest),kind,status,elapsed,err}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'recording','input','output','context','status','seconds','error'});
    writetable(report,fullfile(out,'run_report.csv'));
    if ~isempty(allP),writetable(allP,fullfile(out,'all_detected_peaks.csv'));end
    if ~isempty(allE),writetable(allE,fullfile(heatout,'all_slowdown_events.csv'));end
    fprintf('FINISHED %s %s %.1fs\n',id,status,elapsed);
end
A=cell2table(audit,'VariableNames',{'recording','context','identicalManualAndPriorTrace','maxDffDifference', ...
    'traceCorrelation','priorPeaks','manualPeaks','priorNoiseSigma','manualNoiseSigma', ...
    'priorForwardR','manualForwardR','priorRotSpeedR','manualRotSpeedR'});
writetable(A,fullfile(out,'roi_comparison.csv'));
H=cell2table(hrows,'VariableNames',{'recording','context','acquisitionStatus','status','nPeaks','nEligiblePeaks','nSlowdowns','nPlacedBarrierEpisodes','output','error'});
writetable(H,fullfile(heatout,'run_report.csv'));
heatP=allP(ismember(allP.recording,H.recording),:);
writetable(heatP,fullfile(heatout,'all_detected_peaks.csv'));
if any(report.status=="error"),error('ves041:ManualRerunErrors','Some trials failed; inspect new run_report.csv.');end
review_ves041_remaining_heat_results(heatout);
% Preserve the exact earlier fly-1 comparison and also extend the same
% behavior-only algorithm to all six rectangular heat recordings.
match_ves041_fly1_heat_slowdowns(fullfile(out,'fly1_matched_heat'),sub);
selected=ismember(string(erase({files.name},'_or_z0z1_s_.mat')),rectIds);
match_ves041_fly1_heat_slowdowns(fullfile(out,'rectangular_matched_heat'),sub,files(selected));
% Rebuild the rectangular summary directly from this run's event catalogs.
metrics=cell(0,9);
for id=rectIds
    for turn=0:1
        for context=["heat","control","no heat"]
            e=allE(allE.recording==id&allE.reorientation==turn&allE.context==context,:);
            p=allP.peakTime(allP.recording==id&logical(allP.eligibleAlignedWindow)&allP.context==context);
            hit=false(height(e),1);post=hit;
            for j=1:height(e),dt=p-e.eventTime(j);hit(j)=any(abs(dt)<=1);post(j)=any(dt>=0&dt<=1);end
            metrics(end+1,:)={id,turn,context,height(e),sum(hit),100*mean(hit),sum(post),100*mean(post),mean(e.deltaDff,'omitnan')}; %#ok<AGROW>
        end
    end
end
S=cell2table(metrics,'VariableNames',{'recording','reorientation','context','nEvents','nPeakWithin1s','percentPeakWithin1s','nPeakFollowing1s','percentPeakFollowing1s','meanDeltaDff'});
writetable(S,fullfile(out,'event_metrics.csv'));plot_ves041_rectangular_heat_summary(out);
disp(A);disp(report);
end

function runMaster(input,dest,here)
fullFilePath=input;analysisOutputFolder=dest;analysisOmitImageStack=true; %#ok<NASGU>
run(fullfile(here,'master_imaging_analysis.m'));
end

function [n,sigma]=countPeaks(t,d)
% Same detector as plot_ves041_peak_heat_context; split invalid samples/gaps.
fs=1/median(diff(t));z=d-movmedian(d,max(3,round(30*fs)),'omitmissing');
sigma=1.4826*mad(z(isfinite(z)),1);good=isfinite(d);
beg=find(good&[true;~good(1:end-1)|diff(t)>2/fs]);
ends=find(good&[~good(2:end)|diff(t)>2/fs;true]);n=0;
for j=1:numel(beg)
    ix=beg(j):ends(j);if numel(ix)<3,continue,end
    [~,p]=findpeaks(d(ix),'MinPeakProminence',max(eps,3*sigma), ...
        'MinPeakDistance',max(1,round(2*fs)),'MinPeakWidth',max(1,round(.25*fs)));
    n=n+numel(p);
end
end
