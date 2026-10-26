function audit=verify_ves041_oct19_oct20()
% Refresh heat summaries with current code and audit all local batch outputs.
% This reads processed inputs and does not overwrite those source MAT files.
oldFigureVisible=get(groot,'defaultFigureVisible');
restoreFigures=onCleanup(@()set(groot,'defaultFigureVisible',oldFigureVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat');dir('/Users/sophiarenauld/stacks/2026-10-20/**/*_or_z0z1_s_.mat')];
rows=cell(0,8);
for k=1:numel(files)
    input=fullfile(files(k).folder,files(k).name);out=fullfile(files(k).folder,'smr_analysis');
    kind="unknown";status="unknown";n=0;nj=0;err="";valid=false;
    try
        a=load(input,'dq','roi');assert(size(a.roi.ts,2)==numel(a.dq(1).t),'ROI/DAQ mismatch');
        a=import_ves041_protocol(a,files(k).folder);kind=string(a.protocol.kind);status=string(a.protocol.acquisitionStatus);
        a=process_fictrac_panels_2p(a,.5,.5,strcmp(kind,'bar-jump'));
        if isfield(a,'jump_detected'),nj=sum(a.jump_detected);end
        assert(isfile(fullfile(out,'fly_data.mat')),'Missing saved fly_data.mat');
        assert(isfile(fullfile(out,'fwd_rot_vs_dff_separate.fig')),'Missing primary MATLAB figure');
        if any(strcmp(kind,{'heat-barrier','heat-stripes'}))
            S=analyze_ves041_heat_context(a,out);n=height(S);
        elseif strcmp(kind,'bar-jump')
            plot_ves041_visual_event_responses(a,out,1);
        end
        repair_saved_figure_visibility(out);
        valid=true;
    catch ME
        err=string(getReport(ME,'extended','hyperlinks','off'));fprintf('%s\n',err);
    end
    rows(end+1,:)={string(input),kind,status,valid,n,nj,numel(a.dq(1).t),err}; %#ok<AGROW>
    close all force;
end
audit=cell2table(rows,'VariableNames',{'input','context','acquisitionStatus','outputsVerified','nHeatGeometryEvents','nRecordedBarJumps','nImagingSamples','error'});
writetable(audit,'/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/output_audit.csv');disp(audit);
assert(all(audit.outputsVerified),'Some outputs did not pass verification; see output_audit.csv');
end
