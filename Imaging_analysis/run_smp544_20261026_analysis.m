function report=run_smp544_20261026_analysis()
% Separate exploratory outputs; fixed source normalization, manual ROIs.
% Do not rerun registration, A2P, or overwrite prior figures/source files.
% Registration QC flags remain attached to every processed output.
old=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');here=fileparts(mfilename('fullpath'));
addpath(here);addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
root='/Users/sophiarenauld/stacks/2026-10-26/a2p_exploratory_20261027';
out='/Users/sophiarenauld/stacks/smp544_20261026_analysis_20261028';
assert(~isfolder(out),'Output already exists; choose a fresh run name.');mkdir(out);
files=dir(fullfile(root,'**','*_or_z0z1_s_.mat'));assert(numel(files)==5,'Expected five exploratory A2P inputs');
rows=cell(0,11);
for k=1:numel(files)
    input=fullfile(files(k).folder,files(k).name);id=erase(files(k).name,'_or_z0z1_s_.mat');
    dest=fullfile(files(k).folder,'smp544_manual_analysis_20261028');kind="unknown";nSlow=NaN;tic;
    diary(fullfile(out,[id '.log']));lastwarn('');
    try
        assert(~isfolder(dest),'Per-trial output already exists');
        registration=jsondecode(fileread(fullfile(files(k).folder,[id '_registration_status.json'])));
        a=runOne(input,dest,here);
        assert(isequaln(a.roi.ts,a.roi.pxfilter.manualts),'Saved analysis must use the full manual trace');
        kind=string(a.protocol.kind);
        if ismember(kind,["heat-barrier","heat-pyramid","heat-stripes"])
            plot_ves041_heat_state_traces(a,dest);
            a.heat_analysis=analyze_ves041_heat_context(a,dest);
            E=plot_ves041_heat_exploration(a,fullfile(dest,'heat_exploration'),id);nSlow=height(E);
            subdirs=dir(fullfile(dest,'heat_exploration','**','*.fig'));
            dirs=unique(string({subdirs.folder}));
            for folder=dirs,repair_saved_figure_visibility(char(folder));end
        elseif kind=="bar-jump" && isfield(a.dq(1),'vh')
            % Shared visual-alignment plotting is cell-agnostic; no VES peak
            % detector is applied to the less-peaky SMP544 fluorescence.
            plot_ves041_visual_event_responses(a,dest,true);
        end
        a.analysisRegistrationQc=registration;
        fly_data=a;save(fullfile(dest,'fly_data.mat'),'fly_data','-v7.3');
        repair_saved_figure_visibility(dest);
        t=a.dq(1).t(:);valid=isfinite(a.dq(1).bvf(:))&logical(a.dq(1).bergii.behaviorVolumeValid(:));
        moving=100*mean(abs(a.dq(1).bvf(valid))>=.5);
        note="Exploratory: input registration QC held; fixed baseline preserved";
        status="success";err="";
        n=numel(t);duration=t(end)-t(1);
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));disp(err);
        moving=NaN;n=NaN;duration=NaN;note="Inspect error log";
    end
    [warningText,~]=lastwarn;close all force;diary off;
    rows(end+1,:)={string(id),kind,status,n,duration,moving,nSlow,string(dest),note,string(warningText),err}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'recording','context','status','nFrames','durationSec','movingPercentValidBehavior','nSlowdowns','output','qcNote','lastWarning','error'});
    writetable(report,fullfile(out,'run_report.csv'));fprintf('SMP544 %s %s %.1fs\n',id,status,toc);
end
disp(report);assert(~any(report.status=="error"),'Some trials failed; inspect run report');
end

function a=runOne(input,dest,here)
smp544_input_file=input;smp544_output_directory=dest;smp544_omit_image_stack=true;smp544_roi_source='manual'; %#ok<NASGU>
run(fullfile(here,'smp544_analysis.m'));
a=a2p_data;
end
