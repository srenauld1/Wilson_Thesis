function report=run_ves041_oct24(recordingPrefixes)
% Run all seven completed October 24 imaging inputs from accepted XYZ bundles.
% Outputs stay in each original recording's smr_analysis folder. No raw data,
% ROI masks, registration, or acquisition hardware is modified.
root='/Users/sophiarenauld/stacks/2026-10-24/registration_xyz_fixed_20261026/a2p_ready';
here=fileparts(mfilename('fullpath'));addpath(here);addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
logRoot=fullfile(fileparts(root),'master_imaging_analysis_20261026');if ~isfolder(logRoot),mkdir(logRoot);end
files=dir(fullfile(root,'*','*_s_.mat'));assert(numel(files)==7,'Expected seven completed imaging recordings');
if nargin>0&&~isempty(recordingPrefixes),files=files(startsWith(string({files.name}),string(recordingPrefixes)));end
% Inspect all inputs and geometry before starting the expensive plot batch.
audit=cell(0,8);
for k=1:numel(files)
    a=load(fullfile(files(k).folder,files(k).name),'dq','roi');a=import_ves041_protocol(a,files(k).folder);
    assert(isfield(a,'protocol'),'Missing acquisition protocol');p=a.protocol;
    assert(numel(a.dq(1).t)==size(a.roi.ts,2),'ROI/DAQ lengths differ');
    world="";geometry="";maxError=NaN;
    if startsWith(p.kind,'heat-')
        assert(p.geometryAligned,'Heat geometry alignment failed');
        world=strjoin(unique(p.coordinate_frame(~ismissing(p.coordinate_frame))),';');
        if isfield(p,'heat_geometry'),geometry=strjoin(unique(p.heat_geometry(~ismissing(p.heat_geometry))),';');end
        if isfield(p,'geometryModelMaxErrorV'),maxError=p.geometryModelMaxErrorV;end
        bad=~logical(a.dq(1).bergii.behaviorVolumeValid(:));
        assert(all(isnan(p.pos_x_mm(bad))),'Geometry bridges invalid behavior volumes');
    end
    audit(end+1,:)={string(files(k).name),string(p.kind),world,geometry,numel(a.dq(1).t),sum(~logical(a.dq(1).bergii.behaviorVolumeValid(:))),maxError,string(p.sourceFolder)}; %#ok<AGROW>
end
writetable(cell2table(audit,'VariableNames',{'recording','kind','coordinateFrame','geometries','volumes','invalidBehaviorVolumes','maxGeometryModelErrorV','acquisitionSource'}),fullfile(logRoot,'input_geometry_audit.csv'));
% Exercise bar jumps, pyramid and changing barriers before the remaining runs.
rank=zeros(numel(files),1);for k=1:numel(files),name=files(k).name;if startsWith(name,'20261024_2_3'),rank(k)=1;elseif startsWith(name,'20261024_5_2'),rank(k)=2;elseif startsWith(name,'20261024_4_1'),rank(k)=3;else,rank(k)=10+k;end,end
[~,order]=sort(rank);rows=cell(0,5);
previousVisibility=get(groot,'defaultFigureVisible');cleanup=onCleanup(@()set(groot,'defaultFigureVisible',previousVisibility)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
for n=1:numel(order)
    k=order(n);input=fullfile(files(k).folder,files(k).name);
    diary(fullfile(logRoot,[files(k).name '_analysis.log']));fprintf('START %s\n',input);timer=tic;
    try
        run_one(input,here);status="success";err="";
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));fprintf('%s\n',err);
    end
    elapsed=toc(timer);fprintf('END %s %s (%.1f seconds)\n',files(k).name,status,elapsed);diary off;close all force;
    receipt=dir(fullfile(files(k).folder,'*_behavior_verification.json'));
    metadata=jsondecode(fileread(fullfile(receipt(1).folder,receipt(1).name)));
    rows(end+1,:)={string(input),status,elapsed,err,string(fullfile(metadata.sourceDirectory,'smr_analysis'))}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'input','status','seconds','error','output'});
    writetable(report,fullfile(logRoot,'run_report.csv'));
end
end
function run_one(input,here)
fullFilePath=input;analysisOmitImageStack=true; %#ok<NASGU>
run(fullfile(here,'master_imaging_analysis.m'));
end
