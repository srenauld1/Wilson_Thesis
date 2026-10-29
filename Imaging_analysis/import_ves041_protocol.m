function a=import_ves041_protocol(a,folder)
% Import integrated stimulus geometry using measured tracking-update times.
% CSV frame identifiers are joined to the imported tracking row map. No
% host-clock offset is guessed. State fields use previous-sample interpolation.
% VES041 analyses use full manual masks, never pixel-selected subregions.
if contains(lower(string(folder)),'ves041') && isfield(a,'roi')
    a=select_ves041_manual_roi(a);
end
% Registered analysis bundles retain a verified link to acquisition logs.
sourceFolder=folder;
f=dir(fullfile(sourceFolder,'*_experiment.json'));
if isempty(f)
    verification=dir(fullfile(folder,'*_behavior_verification.json'));
    if numel(verification)==1
        receipt=jsondecode(fileread(fullfile(verification.folder,verification.name)));
        if isfield(receipt,'sourceDirectory') && isfolder(receipt.sourceDirectory)
            sourceFolder=receipt.sourceDirectory;
            f=dir(fullfile(sourceFolder,'*_experiment.json'));
        end
    end
end
if isempty(f),return,end
assert(numel(f)==1,'ves041:AmbiguousProtocol','Multiple experiment logs in %s',sourceFolder);
j=jsondecode(fileread(fullfile(f(1).folder,f(1).name)));
kind='dark';code=lower(j.socketCode);
if contains(code,'heat_pyramid'),kind='heat-pyramid';elseif contains(code,'heat_barrier'),kind='heat-barrier';elseif contains(code,'heat_stripes'),kind='heat-stripes';elseif contains(code,'bar_jump'),kind='bar-jump';end
a.protocol=struct('kind',kind,'acquisitionStatus',j.status,'geometryAligned',false, ...
    'sourceFolder',sourceFolder,'socketCode',j.socketCode);
if ~any(strcmp(kind,{'heat-barrier','heat-pyramid','heat-stripes'})),return,end
if strcmp(kind,'heat-barrier'),pat='*_heat_barrier.csv';elseif strcmp(kind,'heat-pyramid'),pat='*_heat_pyramid.csv';else,pat='*_heat_stripes.csv';end
f=dir(fullfile(sourceFolder,pat));
if isempty(f),warning('ves041:MissingGeometry','No protocol geometry log in %s',folder);return,end
if ~isfield(a.dq(1),'bergii'),warning('ves041:MissingSync','No measured geometry time mapping in %s',folder);return,end
b=a.dq(1).bergii;required={'csvRowIndicesOneBased','updateTimes','volumeTimes'};
if ~all(isfield(b,required)),warning('ves041:MissingSync','Incomplete measured tracking time map in %s',folder);return,end
g=readtable(fullfile(f(1).folder,f(1).name),'VariableNamingRule','preserve');
assert(numel(f)==1,'ves041:AmbiguousGeometry','Multiple geometry logs');
cf=dir(fullfile(sourceFolder,'*_fictrac_frames.csv'));if isempty(cf),warning('ves041:MissingTrackingCSV','Missing tracking CSV in verified source');return,end
c=readtable(fullfile(cf(1).folder,cf(1).name),'VariableNamingRule','preserve');
rows=double(b.csvRowIndicesOneBased(:));assert(all(rows>=1&rows<=height(c)),'Tracking row map exceeds local CSV');
[tf,ix]=ismember(g.frame,c.frame(rows));
assert(sum(tf)>=2,'ves041:GeometryFrameMismatch','Fewer than two geometry rows match the measured tracking map');
u=double(b.updateTimes(:));assert(numel(u)==numel(rows),'Tracking row/time mapping lengths differ');offset=a.dq(1).t(1)-double(b.volumeTimes(1));
gt=nan(height(g),1);gt(tf)=u(ix(tf))+offset;
fields={'trial_number','trial_type','control_trial','barrier_active','distance_to_barrier_mm', ...
    'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm', ...
    'virtual_heat_voltage','requested_voltage','delivered_voltage','heat_active','pos_x_mm','pos_y_mm', ...
    'stripe_phase_mm','distance_to_boundary_mm','in_heat_stripe','heat_width_mm','cool_width_mm', ...
    'stripe_normal_x','stripe_normal_y','stripe_direction_x','stripe_direction_y', ...
    'period_mm','anchor_x_mm','anchor_y_mm','transverse_position_mm', ...
    'coordinate_frame','coordinate_heading_rad','coordinate_yaw_gain','coordinate_step_valid', ...
    'coordinate_position_assumed','heat_profile','heat_geometry','geometry_mode', ...
    'geometry_cycle_index','geometry_pair_number','ramp_distance_mm','peak_voltage_v', ...
    'forward_zero_edge_beyond_peak_mm','forward_escape_beyond_peak_mm', ...
    'backward_escape_before_ramp_start_mm','trial_end_reason','sequence_complete'};
textFields={'trial_type','coordinate_frame','heat_profile','heat_geometry','geometry_mode','trial_end_reason'};
names=intersect(fields,g.Properties.VariableNames,'stable');t=a.dq(1).t(:);
for k=1:numel(names)
    v=g.(names{k});
    if ismember(names{k},textFields)
        valid=isfinite(gt);[tt,ii]=unique(gt(valid));vv=string(v(valid));
        if numel(tt)<2,continue,end
        row=interp1(tt,(1:numel(tt))',t,'previous',NaN);
        z=strings(size(t));z(:)=missing;ok=isfinite(row);z(ok)=vv(ii(row(ok)));
        if isfield(b,'behaviorVolumeValid'),z(~logical(b.behaviorVolumeValid(:)))=missing;end
        a.protocol.(names{k})=z;continue
    end
    if ~isnumeric(v)&&~islogical(v)
        s=lower(string(v));z=str2double(s);z(s=="true")=1;z(s=="false")=0;
        if all(isnan(z)),continue,end
        v=z;
    end
    % Retain explicit NaNs when geometry is unplaced: do not carry a previous
    % subtrial's line through the next lead-in by deleting missing states.
    valid=isfinite(gt);[tt,ii]=unique(gt(valid));vv=double(v(valid));
    if numel(tt)<2,continue,end
    z=interp1(tt,vv(ii),t,'previous',NaN);
    % Do not carry state through tracking outages or beyond logged support.
    if isfield(b,'behaviorVolumeValid'),z(~logical(b.behaviorVolumeValid(:)))=NaN;end
    a.protocol.(names{k})=z;
end
a.protocol.geometryAligned=true;a.protocol.geometryLog=fullfile(f(1).folder,f(1).name);
a.protocol.t=t;a.protocol.matchedGeometryRows=sum(tf);
% New world-aware profiles explicitly log their shape and active geometry.
% Validate the imported virtual command against that saved profile. Older
% turn-ignoring recordings retain their historical geometry unchanged.
if isfield(a.protocol,'peak_voltage_v') && isfield(a.protocol,'heat_profile')
    p=a.protocol;distance=p.distance_to_barrier_mm;
    distance(abs(distance)<1e-9)=0;
    fraction=max(0,1-distance./p.ramp_distance_mm);
    fraction(distance<0)=0;isPyramid=p.heat_profile=="pyramid";
    fraction(isPyramid)=max(0,1-abs(distance(isPyramid))./p.ramp_distance_mm(isPyramid));
    predicted=fraction.*p.peak_voltage_v;
    valid=p.barrier_active>0 & isfinite(predicted) & isfinite(p.virtual_heat_voltage);
    err=abs(predicted(valid)-p.virtual_heat_voltage(valid));
    assert(~isempty(err)&&max(err)<1e-6,'ves041:HeatProfileMismatch', ...
        'Imported geometry does not reproduce logged virtual heat');
    a.protocol.predicted_virtual_heat_voltage=predicted;
    a.protocol.geometryModelMaxErrorV=max(err);
end
% Persist an explicit stimulus-exclusion mask for downstream spontaneous
% analyses; this never changes the fluorescence or kinematic signals.
heat=a.dq(1).heat(:);mask=isfinite(heat)&heat>0.01;
if isfield(a.protocol,'barrier_active'),mask=mask|a.protocol.barrier_active>0;end
a.protocol.stimulusMask=mask;
a.protocol.spontaneousMask=~(movmax(double(mask),[round(3/median(diff(t))) round(3/median(diff(t)))])>0);
if isfield(b,'behaviorVolumeValid'),a.protocol.spontaneousMask=a.protocol.spontaneousMask&logical(b.behaviorVolumeValid(:));end
end
