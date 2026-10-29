function audit=verify_ves041_oct24()
% Check completion, full manual ROI use, geometry, gaps and primary outputs.
root='/Users/sophiarenauld/stacks/2026-10-24/registration_xyz_fixed_20261026';
inputs=dir(fullfile(root,'a2p_ready','*','*_s_.mat'));assert(numel(inputs)==7);
report=readtable(fullfile(root,'master_imaging_analysis_20261026','run_report.csv'), ...
    'TextType','string','Delimiter',',','VariableNamingRule','preserve');
assert(height(report)==7&&all(report.status=="success"),'Batch is incomplete or contains failures');rows=cell(0,7);
for k=1:numel(inputs)
    folder=inputs(k).folder;
    receipt=dir(fullfile(folder,'*_behavior_verification.json'));
    metadata=jsondecode(fileread(fullfile(receipt(1).folder,receipt(1).name)));
    out=fullfile(metadata.sourceDirectory,'smr_analysis');
    original=load(fullfile(folder,inputs(k).name),'roi');original=select_ves041_manual_roi(original);
    saved=load(fullfile(out,'fly_data.mat'),'fly_data');x=saved.fly_data;
    assert(isequaln(x.roi.ts,original.roi.ts),'Saved trace differs from full manual ROI');
    assert(~isfield(x,'stack'),'Analysis output unexpectedly duplicates image stack');
    for name=["fwd_rot_vs_dff_separate","ves041_detected_peaks","dff_peaks_vs_kinematics"]
        assert(isfile(fullfile(out,name+".fig"))&&isfile(fullfile(out,name+".png")),'Missing figure %s',name);
    end
    p=x.protocol;world="";geometry="";
    if startsWith(p.kind,'heat-')
        assert(p.geometryAligned);assert(isfile(fullfile(out,'heat_geometry_summary.png')));
        assert(isfile(fullfile(out,'fwd_rot_vs_dff_separate_heat_on_off.fig')));
        bad=~logical(x.dq(1).bergii.behaviorVolumeValid(:));assert(all(isnan(p.pos_x_mm(bad))));
        world=strjoin(unique(p.coordinate_frame(~ismissing(p.coordinate_frame))),';');
        if isfield(p,'heat_geometry'),geometry=strjoin(unique(p.heat_geometry(~ismissing(p.heat_geometry))),';');end
        if isfield(p,'geometryModelMaxErrorV'),assert(p.geometryModelMaxErrorV<1e-6);end
        ev=readtable(fullfile(out,'heat_geometry_events.csv'),'TextType','string','Delimiter',',','VariableNamingRule','preserve');
        if isfield(p,'heat_geometry'),assert(ismember('heat_geometry',ev.Properties.VariableNames));end
        % Give empty event panels an explicit explanation rather than an
        % empty legend. This refresh changes presentation only.
        if isempty(ev),analyze_ves041_heat_context(x,out);end
    end
    figs=dir(fullfile(out,'*.fig'));pngs=dir(fullfile(out,'*.png'));
    for j=1:numel(figs)
        stored=load(fullfile(figs(j).folder,figs(j).name),'hgS_070000','-mat');
        assert(all(arrayfun(@(v)isfield(v.properties,'Visible')&&strcmp(v.properties.Visible,'on'),stored.hgS_070000)),'Figure will open hidden');
    end
    rows(end+1,:)={string(inputs(k).name),string(p.kind),world,geometry,numel(pngs),numel(figs),string(out)}; %#ok<AGROW>
end
audit=cell2table(rows,'VariableNames',{'recording','protocol','coordinateFrame','geometries','pngCount','figCount','output'});
writetable(audit,fullfile(root,'master_imaging_analysis_20261026','output_audit.csv'));disp(audit);
end
