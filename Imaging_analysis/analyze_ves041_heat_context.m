function summary=analyze_ves041_heat_context(a,out)
% Preliminary barrier/stripe events selected from geometry and behavior.
% Heat is measured command voltage, not temperature. dF/F is not used to
% select events. Save all events, including trials with no eligible events.
% Crossings are consecutive valid changes of stripe membership or barrier
% distance sign while the barrier is active. Repeated crossings are retained.
% An approach/retreat is a distance minimum <=5 mm with >=1 mm approach and
% retreat in the adjacent 1 s windows, no crossing within +/-1 s, and >=15
% degrees net heading displacement for the 'turn away' subset. Plots display
% raw dF/F, with no amplitude rejection or per-event baseline subtraction.
p=a.protocol;t=a.dq(1).t(:);d=a.roi.ts(1,:)';f=a.dq(1).bvf(:);r=abs(a.dq(1).bvy_deg(:));h=a.dq(1).bh(:);
p.analysisSettings=struct('approachDistanceMm',5,'minimumApproachRetreatMm',1, ...
    'approachWindowSeconds',1,'turnAwayNetHeadingDeg',15,'eventWindowSeconds',[-2 2], ...
    'mergeCrossings',false,'baselineSubtractDff',false);
assert(numel(t)==numel(d),'Imaging and DAQ lengths differ');
if ~p.geometryAligned,error('ves041:GeometryUnavailable','Cannot analyze heat geometry without measured alignment');end
if strcmp(p.kind,'heat-stripes')
    distance=p.distance_to_boundary_mm;zone=p.in_heat_stripe;
    idx=find(isfinite(zone(1:end-1))&isfinite(zone(2:end))&diff(zone)~=0)+1;
    labels=repmat("heat entry",numel(idx),1);labels(zone(idx)==0)="heat exit";
else
    distance=p.distance_to_barrier_mm;active=p.barrier_active>0;
    distance(~active)=NaN; % Unplaced-barrier sentinel values are not distances.
    adjacent=active(2:end)&active(1:end-1);
    if isfield(p,'trial_number'),adjacent=adjacent&p.trial_number(2:end)==p.trial_number(1:end-1);end
    idx=find(adjacent&isfinite(distance(2:end))&isfinite(distance(1:end-1))&distance(1:end-1).*distance(2:end)<=0)+1;
    labels=repmat("barrier crossing",numel(idx),1);
    if strcmp(p.kind,'heat-pyramid'),labels(:)="pyramid peak crossing";end
end
% Approaches within 5 mm ending in retreat without crossing the boundary.
dist=abs(distance);[~,mins]=findpeaks(-fillmissing(dist,'constant',Inf),'MinPeakDistance',max(1,round(1/median(diff(t)))));
approach=[];turnaway=[];
for k=mins(:)'
    pre=find(t>=t(k)-1&t<t(k));post=find(t>t(k)&t<=t(k)+1);
    if isempty(pre)||isempty(post)||~isfinite(dist(k))||dist(k)>5,continue,end
    if ~strcmp(p.kind,'heat-stripes')
        if ~all(active([pre;k;post])),continue,end
        if isfield(p,'trial_number')&&numel(unique(p.trial_number([pre;k;post])))~=1,continue,end
    end
    if mean(dist(pre),'omitnan')-dist(k)>=1&&mean(dist(post),'omitnan')-dist(k)>=1
        if any(abs(t(idx)-t(k))<=1),continue,end
        approach(end+1,1)=k;dh=atan2(sin(h(post(end))-h(pre(1))),cos(h(post(end))-h(pre(1))));turnaway(end+1,1)=abs(rad2deg(dh))>=15; %#ok<AGROW>
    end
end
alabels=repmat("approach and retreat",numel(approach),1);alabels(logical(turnaway))="turn away";idx=[idx(:);approach(:)];labels=[labels;alabels];
keep=t(idx)>=t(1)+2&t(idx)<=t(end)-2;
for k=1:numel(idx)
    win=t>=t(idx(k))-2&t<=t(idx(k))+2;
    keep(k)=keep(k)&&all(isfinite(f(win))&isfinite(r(win))&isfinite(distance(win)));
    if isfield(a.dq(1).bergii,'behaviorVolumeValid'),keep(k)=keep(k)&&all(a.dq(1).bergii.behaviorVolumeValid(win));end
    if ~strcmp(p.kind,'heat-stripes')&&isfield(p,'trial_number'),keep(k)=keep(k)&&numel(unique(p.trial_number(win)))==1;end
end
idx=idx(keep);labels=labels(keep);[idx,ord]=sort(idx);labels=labels(ord);
control=nan(numel(idx),1);if isfield(p,'control_trial'),control=p.control_trial(idx);end
if ~strcmp(p.kind,'heat-stripes')
    labels(control==1)=labels(control==1)+" (control)";
    labels(control==0)=labels(control==0)+" (heat)";
end
heat=a.dq(1).heat(:);summary=table(t(idx),labels,control,heat(idx),distance(idx),'VariableNames',{'time','event','controlTrial','measuredHeatCommandV','boundaryDistanceMm'});
metadata={'trial_number','coordinate_frame','heat_geometry','ramp_distance_mm','peak_voltage_v'};
for k=1:numel(metadata)
    if isfield(p,metadata{k}),value=p.(metadata{k});summary.(metadata{k})=value(idx);end
end
% Changing geometries are distinct descriptive groups, never pooled as one
% barrier condition. Events remain aligned to geometry, not calcium peaks.
if isfield(p,'heat_geometry')
    labels=labels+" | "+p.heat_geometry(idx);
    summary.event=labels;
end
writetable(summary,fullfile(out,'heat_geometry_events.csv'));
rel=(-2:.05:2)';vars={d,f,r,distance,heat};labs={'VES041 dF/F','Forward speed (mm/s)','Rotational speed (deg/s)','Boundary distance (mm)','Measured heat command (V)'};
if isfield(a,'analysisCellLabel'),labs{1}=[char(a.analysisCellLabel) ' dF/F'];end
if strcmp(p.kind,'heat-pyramid'),labs{4}='Signed distance to pyramid peak (mm)';end
fig=figure('Visible','off','Color','w','Position',[50 50 1500 1250]);tl=tiledlayout(5,2,'TileSpacing','compact');
for v=1:5
    ax=nexttile(tl);plot(ax,t,vars{v},'k');hold(ax,'on');
    if ~isempty(idx)
        yl=ylim(ax);xx=[t(idx)';t(idx)';nan(1,numel(idx))];yy=repmat([yl(:);NaN],1,numel(idx));
        plot(ax,xx(:),yy(:),'Color',[.1 .6 .2],'LineWidth',.5,'HandleVisibility','off');
    end
    ylabel(ax,labs{v});grid(ax,'on');if v==5,xlabel(ax,'Trial time (s)');end
    ax=nexttile(tl);hold(ax,'on');groups=unique(labels);
    for g=1:numel(groups),sel=find(labels==groups(g));X=nan(numel(rel),numel(sel));for k=1:numel(sel),X(:,k)=interp1(t-t(idx(sel(k))),vars{v},rel,'linear',NaN);end;plot(ax,rel,mean(X,2,'omitnan'),'LineWidth',2,'DisplayName',sprintf('%s (n=%d)',groups(g),numel(sel)));end
    xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-2 2]);ylabel(ax,labs{v});grid(ax,'on');
    if v==1
        if isempty(groups),text(ax,.5,.5,'No eligible geometry events','Units','normalized','HorizontalAlignment','center');
        else,legend(ax,'Location','best');end
    end
    if v==5,xlabel(ax,'Time from geometry event (s)');end
end
[parent,~]=fileparts(out);[~,trial]=fileparts(parent);title(tl,strrep([p.kind ' | ' trial],'_','\_'));
if ~feature('ShowFigureWindows'),savefig_visible_on_open(fig,fullfile(out,'heat_geometry_summary.fig'));
else,savefig(fig,fullfile(out,'heat_geometry_summary.fig'));end
exportgraphics(fig,fullfile(out,'heat_geometry_summary.png'),'Resolution',150);close(fig);
save(fullfile(out,'heat_geometry_results.mat'),'summary','p','-v7.3');
end
