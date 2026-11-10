function events=plot_ves041_12_1_boundary_turns()
%PLOT_VES041_12_1_BOUNDARY_TURNS Map the visually identified stripe excursions.
% Read the SAME smr_analysis/fly_data.mat that supplied xy_heat_paths_01.png.
% Restrict to the spatial region the user identified: x=155:175 mm, t=300:340 s.
% Find logged stripe entry/exit pairs. Retain excursions >=0.15 s to separate
% the six clear dips from two one-volume membership toggles near 322 s.
% These are descriptive, region-selected examples, not an unbiased cohort.
% Align each excursion to maximum rotational SPEED after 0.12-s Gaussian
% smoothing between stripe entry and 0.25 s after exit. No dF/F selection.
% Plot the original dF/F and a 0.25-s moving mean; normalization unchanged.
% Show signed yaw alongside |yaw|, forward velocity, unwrapped heading,
% logged Y, and actual heat monitor. Shade LOGGED stripe membership, not an
% inferred temperature. Sampling/geometry and delivered heat are distinct.
% All plotting is headless; saved FIGs open visibly. Source files untouched.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
source='/Users/sophiarenauld/stacks/2026-10-20/20261020_12_1_ves041_syt8m_fly1_stripes/smr_analysis/fly_data.mat';
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics','ves041_20261020_12_1_boundary_turns_20261102');
if ~isfolder(out),mkdir(out);end
% Read numerical HDF5 fields directly: avoid loading the large image stack.
t=read_field(source,'dq/t');dff=read_field(source,'roi/ts');
fwd=read_field(source,'dq/bvf');yaw=read_field(source,'dq/bvy_deg');
heading=read_field(source,'dq/bh');heat=read_field(source,'dq/heat');
x=read_field(source,'protocol/pos_x_mm');y=read_field(source,'protocol/pos_y_mm');
membership=read_field(source,'protocol/in_heat_stripe');
inside=isfinite(membership)&membership~=0;
manual=read_field(source,'roi/pxfilter/manualts');
assert(isequaln(dff,manual),'Source plotted trace is not the saved manual trace');
behaviorValid=read_field(source,'dq/bergii/behaviorVolumeValid');
valid=isfinite(behaviorValid)&behaviorValid~=0&isfinite(membership);
assert(all([numel(dff),numel(fwd),numel(yaw),numel(heading),numel(heat),numel(x),numel(y),numel(inside),numel(valid)]==numel(t)));
fs=1/median(diff(t));sm=movmean(dff,max(1,round(.25*fs)));
signedRot=smoothdata(yaw,'gaussian',max(3,round(.12*fs)));rot=abs(signedRot);
rows=cell(0,11);briefRows=cell(0,3);
change=diff([false;inside;false]);on=find(change==1);off=find(change==-1)-1;
for k=1:numel(on)
    i=on(k);j=off(k);if t(i)<300||t(i)>340||x(i)<155||x(i)>175,continue,end
    duration=t(j)-t(i)+1/fs;
    if duration<.15,briefRows(end+1,:)={t(i),t(min(j+1,numel(t))),duration};continue,end %#ok<AGROW>
    exit=min(j+1,numel(t));window=find(t>=t(i)&t<=t(exit)+.25&valid&isfinite(rot));
    assert(~isempty(window));[r,best]=max(rot(window));p=window(best);
    dh=heading(exit)-heading(i);net=abs(rad2deg(atan2(sin(dh),cos(dh))));
    rows(end+1,:)={size(rows,1)+1,t(i),t(exit),t(p),duration,x(p),y(p),r,yaw(p),net,p}; %#ok<AGROW>
end
events=cell2table(rows,'VariableNames',{'event','stripeEntrySec','stripeExitSec', ...
    'rotationPeakSec','stripeDurationSec','xAtTurnMm','yAtTurnMm','peakRotSpeedDegS', ...
    'signedYawAtPeakDegS','netHeadingEntryToExitDeg','sampleIndex'});
assert(height(events)==6,'Unexpected number of clear spatial excursions; inspect source.');
brief=cell2table(briefRows,'VariableNames',{'entrySec','exitSec','durationSec'});
writetable(events,fullfile(out,'boundary_turn_times.csv'));
writetable(brief,fullfile(out,'brief_membership_toggles.csv'));
% Overview on native time coordinates. Absolute heading is unwrapped within
% this continuous, valid window only, then set to zero at its first sample.
overview=t>=314&t<=325;ix=find(overview);
assert(all(valid(ix))&&all(isfinite(heading(ix))),'Overview contains missing behavior');
hh=rad2deg(unwrap(heading(ix)));hh=hh-hh(1);
fig=figure('Visible','off','Color','w','Position',[10 10 1500 1250]);
tl=tiledlayout(fig,6,1,'TileSpacing','compact','Padding','compact');axs=gobjects(6,1);
values={dff,fwd,yaw,hh,y,heat};
labels={'dF/F','Forward velocity (mm/s)','Yaw (deg/s)','Unwrapped heading (deg)','Logged Y (mm)','Heat monitor (V)'};
for v=1:6
    ax=nexttile(tl);axs(v)=ax;hold(ax,'on');
    if v==4,plot(ax,t(ix),values{v},'k');
    elseif v==1
        plot(ax,t(ix),dff(ix),'Color',[.65 .75 .9]);
        plot(ax,t(ix),sm(ix),'Color',[.05 .25 .75],'LineWidth',1.3);
    elseif v==3
        plot(ax,t(ix),yaw(ix),'Color',[.65 .65 .65]);
        plot(ax,t(ix),abs(yaw(ix)),'Color',[.85 .1 .1],'LineWidth',1.1);
    else,plot(ax,t(ix),values{v}(ix),'k','LineWidth',1.1);end
    ylabel(ax,labels{v});grid(ax,'on');xlim(ax,[314 325]);
    if v==3,yline(ax,0,':k');end
    shade_membership(ax,t(ix),inside(ix));
    for j=1:height(events)
        xline(ax,events.rotationPeakSec(j),'--','Color',[.35 .2 .6]);
        if v==1
            yl=ylim(ax);text(ax,events.rotationPeakSec(j),yl(2),sprintf(' %d',j), ...
                'VerticalAlignment','top','Color',[.35 .2 .6],'FontWeight','bold');
        end
    end
end
linkaxes(axs,'x');xlabel(axs(6),'Time from recording start (s)');
title(tl,'20261020 12_1 | repeated stripe-boundary excursions near X=155-175 mm');
subtitle(tl,'Orange: logged heat-stripe membership | purple dashed: local rotational-speed maximum | dF/F: raw + 0.25 s mean | yaw: signed gray, speed red');
save_plot(fig,out,'boundary_turns_time_traces');

% One aligned column PER excursion, with actual time printed above it.
fig=figure('Visible','off','Color','w','Position',[10 10 2150 1050]);
tl=tiledlayout(fig,4,height(events),'TileSpacing','compact','Padding','compact');
for v=1:4
    rowAxes=gobjects(1,height(events));
    for j=1:height(events)
        ax=nexttile(tl);rowAxes(j)=ax;hold(ax,'on');tr=t-events.rotationPeakSec(j);w=abs(tr)<=2.1;
        switch v
            case 1,plot(ax,tr(w),dff(w),'Color',[.65 .75 .9]);plot(ax,tr(w),sm(w),'Color',[.05 .25 .75],'LineWidth',1.3);ylabel(ax,'dF/F');
            case 2,plot(ax,tr(w),fwd(w),'k');ylabel(ax,'Forward (mm/s)');
            case 3,plot(ax,tr(w),abs(yaw(w)),'Color',[.85 .1 .1]);ylabel(ax,'Rot. speed (deg/s)');
            case 4,plot(ax,tr(w),heat(w),'Color',[.8 .35 .1]);ylabel(ax,'Heat (V)');
        end
        shade_membership(ax,tr(w),inside(w));xline(ax,0,'--','Color',[.35 .2 .6]);xlim(ax,[-2 2]);grid(ax,'on');
        if v==1,title(ax,sprintf('Turn %d | %.3f s',j,events.rotationPeakSec(j)));end
        if v==4,xlabel(ax,'Seconds from turn');end
    end
    linkaxes(rowAxes,'y');
end
title(tl,'20261020 12_1 | each identified excursion aligned to its rotational-speed maximum');
subtitle(tl,'Overlapping +/-2 s windows show neighboring excursions too | same original dF/F as the supplied flat path; no baseline subtraction');
save_plot(fig,out,'boundary_turns_individual_aligned');

% Spatial check linking each numbered time marker to the supplied path.
fig=figure('Visible','off','Color','w','Position',[10 10 1300 500]);ax=axes(fig);hold(ax,'on');
patch(ax,[153 176 176 153],[-9.5 -9.5 -7.5 -7.5],[1 .85 .8],'EdgeColor','none');
z=t>=315&t<=324&x>=153&x<=176;
surface(ax,[x(z) x(z)],[y(z) y(z)],zeros(nnz(z),2),[dff(z) dff(z)],'FaceColor','none','EdgeColor','interp','LineWidth',2);
for j=1:height(events)
    plot(ax,events.xAtTurnMm(j),events.yAtTurnMm(j),'ko','MarkerSize',7);
    text(ax,events.xAtTurnMm(j),events.yAtTurnMm(j)+.12,sprintf('%d: %.2f s',j,events.rotationPeakSec(j)),'FontSize',9);
end
cl=prctile(dff(isfinite(dff)),[1 99]);caxis(ax,cl);colormap(ax,parula);cb=colorbar(ax);cb.Label.String='dF/F';
axis(ax,'equal');xlim(ax,[153 176]);ylim(ax,[-9.5 -5.5]);grid(ax,'on');xlabel(ax,'Logged X (mm)');ylabel(ax,'Logged Y (mm)');
title(ax,'Numbered turn times on the original dF/F-colored path | stripe boundary Y=-7.5 mm');
save_plot(fig,out,'boundary_turns_path_times');
save(fullfile(out,'boundary_turns_data.mat'),'events','brief','source','t','dff','sm','fwd','yaw','heat','x','y','inside','valid');
disp(events);
end

function data=read_field(file,field)
% Saved MATLAB structs may store numeric fields directly or via an HDF5
% object reference. Resolve the latter without reading unrelated image data.
f=H5F.open(file,'H5F_ACC_RDONLY','H5P_DEFAULT');cf=onCleanup(@()H5F.close(f)); %#ok<NASGU>
% dq.bergii is itself a reference to a nested struct, unlike protocol.
% Resolve that intermediate struct before opening its numeric field.
path=['/fly_data/' field];
if startsWith(field,'dq/bergii/')
    parent=H5D.open(f,'/fly_data/dq/bergii');
    cp=onCleanup(@()H5D.close(parent)); %#ok<NASGU>
    ref=H5D.read(parent);g=H5R.dereference(parent,'H5R_OBJECT',ref(:,1));
    cg=onCleanup(@()H5G.close(g)); %#ok<NASGU>
    path=[H5I.get_name(g) '/' extractAfter(field,'dq/bergii/')];
end
d=H5D.open(f,path);cd=onCleanup(@()H5D.close(d)); %#ok<NASGU>
type=H5D.get_type(d);ct=onCleanup(@()H5T.close(type)); %#ok<NASGU>
if H5T.get_class(type)==H5ML.get_constant_value('H5T_REFERENCE')
    ref=H5D.read(d);target=H5R.dereference(d,'H5R_OBJECT',ref(:,1));
    closeTarget=onCleanup(@()H5D.close(target)); %#ok<NASGU>
    data=H5D.read(target);
else,data=H5D.read(d);end
data=double(data(:));
end
function shade_membership(ax,t,in)
yl=ylim(ax);dd=diff([false;in(:);false]);s=find(dd==1);e=find(dd==-1)-1;
for k=1:numel(s)
    right=t(e(k));if e(k)<numel(t),right=t(e(k)+1);end
    p=patch(ax,[t(s(k)) right right t(s(k))],[yl(1) yl(1) yl(2) yl(2)], ...
        [1 .55 .25],'FaceAlpha',.12,'EdgeColor','none');uistack(p,'bottom');
end
ylim(ax,yl);
end
function save_plot(fig,out,stem)
exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',170);
savefig_visible_on_open(fig,fullfile(out,[stem '.fig']));close(fig);
end
