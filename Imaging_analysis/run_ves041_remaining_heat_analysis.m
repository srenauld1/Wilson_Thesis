function report=run_ves041_remaining_heat_analysis()
% Extend the existing fly-1 heat analyses to 10/19 and 10/20 fly 2.
% Use current dff005015 inputs; preserve all previous smr_analysis results.
% Plot/detection/slowdown thresholds are inherited unchanged from the heat
% functions. Non-heat contexts are logged as skipped, not analyzed as heat.
oldVisible=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@()set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
here=fileparts(mfilename('fullpath'));addpath(here);addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
out=fullfile('/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/15s_ave_f','remaining_heat');
if ~isfolder(out),mkdir(out);end
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat'); ...
    dir('/Users/sophiarenauld/stacks/2026-10-20/*fly2*/*_or_z0z1_s_.mat')];
rows=cell(0,10);allPeaks=table;allSlowdowns=table;
for k=1:numel(files)
    input=fullfile(files(k).folder,files(k).name);id=erase(files(k).name,'_or_z0z1_s_.mat');
    dest=fullfile(files(k).folder,'15s_ave_f');kind="unknown";acquisition="unknown";nPeaks=0;nEligible=0;nSlowdowns=0;nBarrierEpisodes=NaN;
    fprintf('HEAT REVIEW %d/%d %s\n',k,numel(files),id);
    try
        a=load(input,'dq','roi');a=import_ves041_protocol(a,files(k).folder);
        assert(isfield(a,'protocol'),'Missing local protocol metadata');
        kind=string(a.protocol.kind);acquisition=string(a.protocol.acquisitionStatus);
        if ~ismember(kind,["heat-barrier","heat-stripes"])
            rows(end+1,:)={string(id),kind,acquisition,"skipped non-heat",0,0,0,NaN,string(dest),""}; %#ok<AGROW>
            continue
        end
        assert(strcmp(a.roi.opt.nrm.nrmstr,'dff005015'),'Expected updated 15 s dF/F');
        assert(a.protocol.geometryAligned,'Missing measured geometry alignment');
        if kind=="heat-barrier"
            p=a.protocol;nBarrierEpisodes=numel(unique(p.trial_number(p.barrier_active>0 & isfinite(p.trial_number) & p.trial_number>0)));
        end
        if ~isfolder(dest),mkdir(dest);end
        a.peakTraceSource=input;
        plot_ves041_heat_state_traces(a,dest);
        P=plot_ves041_peak_heat_context(a,dest);nPeaks=height(P);nEligible=sum(P.eligibleAlignedWindow);
        P=addvars(P,repmat(string(id),height(P),1),'Before',1,'NewVariableNames','recording');
        allPeaks=[allPeaks;P]; %#ok<AGROW>
        E=plot_ves041_fly1_heat_exploration(files(k).folder,fullfile(out,id),files(k),'15s_ave_f');
        nSlowdowns=height(E);allSlowdowns=[allSlowdowns;E]; %#ok<AGROW>
        repair_saved_figure_visibility(dest);repair_saved_figure_visibility(fullfile(dest,'heat_exploration'));
        status="success";err="";
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));fprintf('%s\n',err);
    end
    close all force;
    rows(end+1,:)={string(id),kind,acquisition,status,nPeaks,nEligible,nSlowdowns,nBarrierEpisodes,string(dest),err}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'recording','context','acquisitionStatus','status','nPeaks','nEligiblePeaks','nSlowdowns','nPlacedBarrierEpisodes','output','error'});
    writetable(report,fullfile(out,'run_report.csv'));
    if ~isempty(allPeaks),writetable(allPeaks,fullfile(out,'all_detected_peaks.csv'));end
    if ~isempty(allSlowdowns),writetable(allSlowdowns,fullfile(out,'all_slowdown_events.csv'));end
end
report=cell2table(rows,'VariableNames',{'recording','context','acquisitionStatus','status','nPeaks','nEligiblePeaks','nSlowdowns','nPlacedBarrierEpisodes','output','error'});
writetable(report,fullfile(out,'run_report.csv'));disp(report);
end
