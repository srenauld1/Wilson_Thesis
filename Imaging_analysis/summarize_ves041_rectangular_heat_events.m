function S = summarize_ves041_rectangular_heat_events()
% Compare identical behavior/peak metrics in rectangular Oct 19/20 trials.
% Rectangular dimensions (512 x 256) were checked in TIFF first-frame tags.
% Behavior catalogs are unchanged. All fluorescence outcomes and peaks are
% recomputed from full hand-drawn ROI traces with dff005015 normalization.
% Peak association requires an eligible peak in the SAME stable context
% within +/-1 s of the slowdown. This is descriptive and not chance-adjusted.
old=get(groot,'defaultFigureVisible');
c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
base='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/15s_ave_f';
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics','rectangular_heat_metrics','hand_drawn');
if ~isfolder(out),mkdir(out);end
E=readtable(fullfile(base,'remaining_heat','all_slowdown_events.csv'),'TextType','string');
F=readtable(fullfile(base,'fly1_matched_heat','candidate_events.csv'),'TextType','string');
ids=["20261019_5_1","20261019_7_1","20261019_8_2", ...
    "20261020_11_2","20261020_12_1","20261020_13_1"];
rows=cell(0,9);
for id=ids
    if startsWith(id,'20261019')
        e=E(E.recording==id,:);
        files=dir(fullfile('/Users/sophiarenauld/stacks/2026-10-19','**',id+'_or_z0z1_s_.mat'));
    else
        files=dir(fullfile('/Users/sophiarenauld/stacks/2026-10-20',id+'*',id+'_or_z0z1_s_.mat'));
        e=F(F.recording==id,:);
    end
        assert(numel(files)==1,'Expected exactly one input for %s',id);
        a=load(fullfile(files.folder,files.name),'dq','roi');
        assert(strcmp(a.roi.opt.nrm.nrmstr,'dff005015'));
        a=import_ves041_protocol(a,files.folder);
        assert(isequaln(a.roi.ts,a.roi.pxfilter.manualts),'Expected full manual ROI');
        % Reuse the same detector and eligibility/context code as Oct 19.
        p=plot_ves041_peak_heat_context(a,fullfile(out,id));
        close all;
        repair_saved_figure_visibility(fullfile(out,id));
        % Behavior selection stays fixed; refresh every fluorescence outcome.
        t=a.dq(1).t(:);d=double(a.roi.ts(:));
        for j=1:height(e)
            et=e.eventTime(j);
            e.baselineDff(j)=mean(d(t>=et-2&t<et-1));
            e.responseDff(j)=mean(d(t>=et&t<=et+1));
            e.deltaDff(j)=e.responseDff(j)-e.baselineDff(j);
        end
        writetable(e,fullfile(out,id,'manual_roi_slowdown_events.csv'));
    for turn=0:1
        for context=["heat","control","no heat"]
            q=e(e.reorientation==turn&e.context==context,:);
            times=p.peakTime(logical(p.eligibleAlignedWindow)&p.context==context);
            hit=false(height(q),1);post=hit;
            for j=1:height(q)
                dt=times-q.eventTime(j);
                hit(j)=any(abs(dt)<=1);post(j)=any(dt>=0&dt<=1);
            end
            rows(end+1,:)={id,turn,context,height(q),sum(hit), ...
                100*mean(hit),sum(post),100*mean(post),mean(q.deltaDff,'omitnan')}; %#ok<AGROW>
        end
    end
end
S=cell2table(rows,'VariableNames',{'recording','reorientation','context', ...
    'nEvents','nPeakWithin1s','percentPeakWithin1s','nPeakFollowing1s', ...
    'percentPeakFollowing1s','meanDeltaDff'});
writetable(S,fullfile(out,'event_metrics.csv'));
disp(S(S.reorientation==1,:));
end
