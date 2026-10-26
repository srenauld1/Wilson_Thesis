function a=import_ves041_protocol(a,folder)
% Import integrated stimulus geometry using measured tracking-update times.
% CSV frame identifiers are joined to the imported tracking row map. No
% host-clock offset is guessed. State fields use previous-sample interpolation.
% VES041 analyses use full manual masks, never pixel-selected subregions.
if contains(lower(string(folder)),'ves041') && isfield(a,'roi')
    a=select_ves041_manual_roi(a);
end
f=dir(fullfile(folder,'*_experiment.json'));
if isempty(f),return,end
j=jsondecode(fileread(fullfile(f(1).folder,f(1).name)));
kind='dark';code=lower(j.socketCode);
if contains(code,'heat_barrier'),kind='heat-barrier';elseif contains(code,'heat_stripes'),kind='heat-stripes';elseif contains(code,'bar_jump'),kind='bar-jump';end
a.protocol=struct('kind',kind,'acquisitionStatus',j.status,'geometryAligned',false);
if ~any(strcmp(kind,{'heat-barrier','heat-stripes'})),return,end
if strcmp(kind,'heat-barrier'),pat='*_heat_barrier.csv';else,pat='*_heat_stripes.csv';end
f=dir(fullfile(folder,pat));
if isempty(f),warning('ves041:MissingGeometry','No protocol geometry log in %s',folder);return,end
if ~isfield(a.dq(1),'bergii'),warning('ves041:MissingSync','No measured geometry time mapping in %s',folder);return,end
b=a.dq(1).bergii;required={'csvRowIndicesOneBased','updateTimes','volumeTimes'};
if ~all(isfield(b,required)),warning('ves041:MissingSync','Incomplete measured tracking time map in %s',folder);return,end
g=readtable(fullfile(f(1).folder,f(1).name),'VariableNamingRule','preserve');
cf=dir(fullfile(folder,'*_fictrac_frames.csv'));if isempty(cf),warning('ves041:MissingTrackingCSV','Missing local tracking CSV');return,end
c=readtable(fullfile(cf(1).folder,cf(1).name),'VariableNamingRule','preserve');
rows=double(b.csvRowIndicesOneBased(:));assert(all(rows>=1&rows<=height(c)),'Tracking row map exceeds local CSV');
[tf,ix]=ismember(g.frame,c.frame(rows));
assert(sum(tf)>=2,'ves041:GeometryFrameMismatch','Fewer than two geometry rows match the measured tracking map');
u=double(b.updateTimes(:));offset=a.dq(1).t(1)-double(b.volumeTimes(1));
gt=nan(height(g),1);gt(tf)=u(ix(tf))+offset;
fields={'trial_number','trial_type','control_trial','barrier_active','distance_to_barrier_mm', ...
    'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm', ...
    'virtual_heat_voltage','requested_voltage','delivered_voltage','heat_active','pos_x_mm','pos_y_mm', ...
    'stripe_phase_mm','distance_to_boundary_mm','in_heat_stripe','heat_width_mm','cool_width_mm', ...
    'stripe_normal_x','stripe_normal_y','stripe_direction_x','stripe_direction_y'};
names=intersect(fields,g.Properties.VariableNames,'stable');t=a.dq(1).t(:);
for k=1:numel(names)
    v=g.(names{k});
    if ~isnumeric(v)&&~islogical(v)
        s=lower(string(v));z=str2double(s);z(s=="true")=1;z(s=="false")=0;
        if all(isnan(z)),continue,end
        v=z;
    end
    valid=isfinite(gt)&isfinite(double(v));[tt,ii]=unique(gt(valid));vv=double(v(valid));
    if numel(tt)<2,continue,end
    z=interp1(tt,vv(ii),t,'previous',NaN);
    % Do not carry state through tracking outages or beyond logged support.
    if isfield(b,'behaviorVolumeValid'),z(~logical(b.behaviorVolumeValid(:)))=NaN;end
    a.protocol.(names{k})=z;
end
a.protocol.geometryAligned=true;a.protocol.geometryLog=fullfile(f(1).folder,f(1).name);
a.protocol.t=t;a.protocol.matchedGeometryRows=sum(tf);
% Persist an explicit stimulus-exclusion mask for downstream spontaneous
% analyses; this never changes the fluorescence or kinematic signals.
heat=a.dq(1).heat(:);mask=isfinite(heat)&heat>0.01;
if isfield(a.protocol,'barrier_active'),mask=mask|a.protocol.barrier_active>0;end
a.protocol.stimulusMask=mask;
a.protocol.spontaneousMask=~(movmax(double(mask),[round(3/median(diff(t))) round(3/median(diff(t)))])>0);
if isfield(b,'behaviorVolumeValid'),a.protocol.spontaneousMask=a.protocol.spontaneousMask&logical(b.behaviorVolumeValid(:));end
end
