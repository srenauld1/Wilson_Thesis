function E=plot_ves041_heat_exploration(a,out,name)
% Heat exploration using full manual ROI and measured protocol coordinates.
% Same behavior-selected events and outputs as the October 20 exploration.
if nargin<3,[folder,~]=fileparts(out);[~,name]=fileparts(folder);end
assert(isfield(a,'protocol')&&a.protocol.geometryAligned,'Measured geometry required');
assert(ismember(a.protocol.kind,{'heat-barrier','heat-pyramid','heat-stripes'}),'Not a heat recording');
if ~isfolder(out),mkdir(out);end
t=a.dq(1).t(:);d=a.roi.ts(1,:)';p=a.protocol;
assert(numel(t)==numel(d)&&all(diff(t)>0),'Invalid trace alignment');
cl=prctile(d(isfinite(d)),[1 99]);if cl(1)==cl(2),cl=cl+[-.1 .1];end
paths(p,t,d,cl,out,name);
E=slowdowns(a,name);
metadata={'trial_number','coordinate_frame','heat_geometry','ramp_distance_mm','peak_voltage_v'};
[~,ix]=ismember(E.eventTime,t);
for j=1:numel(metadata)
    key=metadata{j};
    if isfield(p,key),values=p.(key);E.(key)=values(ix);
    elseif ismember(key,{'coordinate_frame','heat_geometry'}),E.(key)=strings(height(E),1);
    else,E.(key)=nan(height(E),1);end
end
writetable(E,fullfile(out,'slowdown_events.csv'));
plotName=name;
if isfield(p,'heat_geometry')
    geometries=unique(p.heat_geometry(~ismissing(p.heat_geometry)));
    if numel(geometries)>1,plotName=[name ' | all geometries pooled; see per-geometry folders'];end
end
eventplots(a,E,out,plotName);
if exist('geometries','var')&&numel(geometries)>1
    for j=1:numel(geometries)
        sub=fullfile(out,char(geometries(j)));if ~isfolder(sub),mkdir(sub);end
        subset=E(E.heat_geometry==geometries(j),:);
        writetable(subset,fullfile(sub,'slowdown_events.csv'));
        eventplots(a,subset,sub,[name ' | ' char(geometries(j))]);
    end
end
settings=struct('alignment','Maximum deceleration selected from behavior only', ...
    'smoothWindowSec',.25,'maximumAccelerationMmS2',-2,'minimumSpeedBeforeMmS',1, ...
    'minimumSpeedDropMmS',1,'minimumAcceptedEventSpacingSec',2, ...
    'speedBeforeWindowSec',[-.75 -.25],'speedAfterWindowSec',[.25 .75], ...
    'contextWindowSec',[-.5 .5],'eventWindowSec',[-2 3], ...
    'reorientationAbsoluteYawDeg',30,'heatCommandThresholdV',.05, ...
    'deltaDffBaselineWindowSec',[-2 -1],'deltaDffResponseWindowSec',[0 1], ...
    'source','Supplied imaging ROI (full manual ROI for VES041); no fluorescence amplitude exclusion', ...
    'behaviorMatched',false,'hypothesisTests',false);
if isfield(a,'analysisRoiSource'),settings.source=char(a.analysisRoiSource);end
fid=fopen(fullfile(out,'analysis_settings.json'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true)); %#ok<NASGU>
fprintf('%s: %d behavior-selected slowdowns\n',name,height(E));
end

function paths(p,t,d,cl,out,name)
x=p.pos_x_mm(:);y=p.pos_y_mm(:);ok=isfinite(x)&isfinite(y)&isfinite(d);
if ~strcmp(p.kind,'heat-stripes')
    ids=unique(p.trial_number(p.barrier_active>0&isfinite(p.trial_number)&p.trial_number>0))';
    sets=arrayfun(@(v)find(ok&p.trial_number==v&p.barrier_active>0),ids,'UniformOutput',false);
    labels=arrayfun(@(v)sprintf('%s episode %g',erase(p.kind,'heat-'),v),ids,'UniformOutput',false);
else
    good=find(ok&isfinite(p.stripe_phase_mm));
    edges=t(good(1)):120:t(good(end));edges=[edges t(good(end))+eps(t(good(end)))];
    sets={good};labels={'Full recording'};
    for j=1:numel(edges)-1,sets{end+1}=find(ok&isfinite(p.stripe_phase_mm)&t>=edges(j)&t<edges(j+1));labels{end+1}=sprintf('%.0f-%.0f s',edges(j:j+1));end %#ok<AGROW>
end
for page=1:ceil(numel(sets)/6)
    fig=figure('Visible','off','Color','w','Position',[40 40 1500 950]);tl=tiledlayout(fig,2,3,'TileSpacing','compact');
    for j=(page-1)*6+1:min(page*6,numel(sets))
        ax=nexttile(tl);ix=sets{j};if isempty(ix),continue,end
        hold(ax,'on');bx=[min(x(ix)) max(x(ix))];by=[min(y(ix)) max(y(ix))];
        if ~strcmp(p.kind,'heat-stripes')
            q=ix(1);mx=median(x(ix));my=median(y(ix));
            distance=(p.barrier_x_mm(q)-mx)*p.barrier_normal_x(q)+(p.barrier_y_mm(q)-my)*p.barrier_normal_y(q);
            footx=mx+distance*p.barrier_normal_x(q);footy=my+distance*p.barrier_normal_y(q);
            bx=[min(bx(1),footx) max(bx(2),footx)];by=[min(by(1),footy) max(by(2),footy)];
        end
        pad=max(3,.08*max([diff(bx) diff(by)]));bx=bx+[-pad pad];by=by+[-pad pad];
        [X,Y]=meshgrid(linspace(bx(1),bx(2),250),linspace(by(1),by(2),250));q=ix(1);
        if ~strcmp(p.kind,'heat-stripes')
            dist=(p.barrier_x_mm(q)-X)*p.barrier_normal_x(q)+(p.barrier_y_mm(q)-Y)*p.barrier_normal_y(q);
            Z=max(0,1-abs(dist)/p.barrier_radius_mm(q));
            if isfield(p,'heat_profile')&&p.heat_profile(q)=="barrier"
                Z=max(0,1-dist/p.ramp_distance_mm(q));Z(dist<0)=0;
            end
            if isfield(p,'peak_voltage_v'),Z=Z*p.peak_voltage_v(q)/5;end
            control=p.control_trial(q)==1;
            if control,rgb=repmat(1-.18*Z,1,1,3);suffix=' | control (virtual gradient)';else,rgb=cat(3,ones(size(Z)),1-.24*Z,1-.32*Z);suffix=' | heat';end
            image(ax,bx,by,rgb);contour(ax,X,Y,dist,[0 0],'Color',[.55 .2 .15],'LineWidth',1.2);
            if isfield(p,'heat_geometry')
                suffix=[suffix sprintf(' | %s: %g mm / %g V',p.heat_geometry(q),p.ramp_distance_mm(q),p.peak_voltage_v(q))];
            end
        else
            period=p.heat_width_mm(q)+p.cool_width_mm(q);nx=p.stripe_normal_x(q);ny=p.stripe_normal_y(q);
            phase=mod((X-x(q))*nx+(Y-y(q))*ny+p.stripe_phase_mm(q),period);
            Z=double(phase>=p.cool_width_mm(q));rgb=cat(3,ones(size(Z)),1-.15*Z,1-.22*Z);
            image(ax,bx,by,rgb);contour(ax,X,Y,Z,[.5 .5],'Color',[.75 .55 .5],'LineWidth',.5);suffix='';
            % Check reconstructed zone membership against the aligned log.
            ph=mod((x(ix)-x(q))*nx+(y(ix)-y(q))*ny+p.stripe_phase_mm(q),period);
            assert(mean((ph>=p.cool_width_mm(q))==logical(p.in_heat_stripe(ix)))>.98,'Stripe geometry does not match logged membership');
        end
        % Surface interpolates fluorescence along segments; NaNs break outages.
        xx=x(ix);yy=y(ix);dd=d(ix);gap=[false;diff(ix)>1];xx(gap)=NaN;yy(gap)=NaN;
        surface(ax,[xx xx],[yy yy],zeros(numel(ix),2),[dd dd],'FaceColor','none','EdgeColor','interp','LineWidth',1.8);
        plot(ax,x(ix(1)),y(ix(1)),'go','MarkerFaceColor','g','MarkerSize',6);plot(ax,x(ix(end)),y(ix(end)),'ks','MarkerFaceColor','k','MarkerSize',5);
        axis(ax,'equal');xlim(ax,bx);ylim(ax,by);set(ax,'YDir','normal');colormap(ax,turbo);clim(ax,cl);cb=colorbar(ax);cb.Label.String='dF/F';
        xlabel(ax,'Logged X (mm)');ylabel(ax,'Logged Y (mm)');title(ax,strrep([labels{j} suffix],'_',' '),'FontSize',9);
    end
    title(tl,{strrep(name,'_',' '),'Background: heat command (gray: virtual control), not temperature | green: start; square: end'});
    saveplot(fig,out,sprintf('xy_heat_paths_%02d',page));
end
end

function E=slowdowns(a,name)
t=a.dq(1).t(:);d=a.roi.ts(1,:)';f=a.dq(1).bvf(:);p=a.protocol;h=a.dq(1).heat(:);
if isfield(a.dq(1),'bvy_deg'),r=a.dq(1).bvy_deg(:);else,r=rad2deg(a.dq(1).bvy(:));end
fs=1/median(diff(t));n=max(3,round(.25*fs));f=smoothdata(f,'movmedian',n);r=smoothdata(r,'movmedian',n);acc=gradient(f,t);
valid=isfinite(f)&isfinite(r)&isfinite(h)&isfinite(d)&isfinite(p.pos_x_mm(:));
if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid'),valid=valid&logical(a.dq(1).bergii.behaviorVolumeValid(:));end
candidate=valid&acc<=-2;starts=find(diff([false;candidate])==1);stops=find(diff([candidate;false])==-1);
idx=zeros(size(starts));for j=1:numel(starts),[~,q]=min(acc(starts(j):stops(j)));idx(j)=starts(j)+q-1;end
rows=cell(0,12);last=-Inf;
for q=idx'
    if t(q)-last<2||t(q)<t(1)+2||t(q)>t(end)-3,continue,end
    win=t>=t(q)-2&t<=t(q)+3;pre=t>=t(q)-.75&t<t(q)-.25;post=t>=t(q)+.25&t<=t(q)+.75;core=t>=t(q)-.5&t<=t(q)+.5;
    if ~all(valid(win))||any(diff(t(win))>2/fs),continue,end
    before=median(f(pre));after=median(f(post));drop=before-after;
    if before<1||drop<1,continue,end
    turn=trapz(t(core),abs(r(core)));hasTurn=turn>=30;
    on=h(core)>.05;
    if all(on),context="heat";elseif all(~on),context="no heat";else,context="transition";end
    control=false;
    if ~strcmp(p.kind,'heat-stripes')
        if any(~isfinite(p.control_trial(core)))||numel(unique(p.trial_number(core)))~=1,context="transition";
        elseif all(p.control_trial(core)==1),control=true;context="control";end
    end
    b=t>=t(q)-2&t<t(q)-1;resp=t>=t(q)&t<=t(q)+1;
    rows(end+1,:)={string(name),t(q),context,control,hasTurn,before,after,drop,turn,mean(d(b)),mean(d(resp)),mean(d(resp))-mean(d(b))}; %#ok<AGROW>
    last=t(q);
end
E=cell2table(rows,'VariableNames',{'recording','eventTime','context','control','reorientation','speedBefore','speedAfter','speedDrop','absoluteTurnDeg','baselineDff','responseDff','deltaDff'});
E.recording=string(E.recording);E.context=string(E.context);E.control=logical(E.control);E.reorientation=logical(E.reorientation);
end

function eventplots(a,E,out,name)
t=a.dq(1).t(:);d=a.roi.ts(1,:)';f=a.dq(1).bvf(:);h=a.dq(1).heat(:);
if isfield(a.dq(1),'bvy_deg'),r=abs(a.dq(1).bvy_deg(:));else,r=abs(rad2deg(a.dq(1).bvy(:)));end
rel=(-2:.05:3)';contexts=["no heat" "heat" "transition" "control"];
X=nan(numel(rel),height(E),4);for j=1:height(E),for v=1:4,vars={d,f,r,h};X(:,j,v)=interp1(t,vars{v},E.eventTime(j)+rel,'linear',NaN);end,end
if height(E)>0,lim=prctile(reshape(X(:,:,1),[],1),[2 98]);else,lim=[0 1];end
if any(~isfinite(lim))||diff(lim)<=0,lim=[0 1];end
for mode=1:2
    if mode==1,cs=contexts(1:2);stem='slowdown_heat_comparison';else,cs=contexts(3:4);stem='slowdown_transition_control';end
    fig=figure('Visible','off','Color','w','Position',[20 20 1550 1250]);tl=tiledlayout(fig,5,4,'TileSpacing','compact');
    for c=1:4
        ctx=cs(ceil(c/2));turn=mod(c,2)==0;sel=find(E.context==ctx&E.reorientation==turn);
        for row=1:5
            ax=nexttile(tl,(row-1)*4+c);hold(ax,'on');
            if row==1
                if ~isempty(sel),imagesc(ax,rel,1:numel(sel),X(:,sel,1)');set(ax,'YDir','normal');clim(ax,lim);colormap(ax,parula);colorbar(ax);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
                ylabel(ax,'Event (chronological)');title(ax,sprintf('%s | turn=%d | n=%d',ctx,turn,numel(sel)));xlim(ax,[-2 3]);
            else
                v=row-1;
                if ~isempty(sel),plot(ax,rel,X(:,sel,v),'Color',[.8 .8 .8]);plot(ax,rel,mean(X(:,sel,v),2,'omitnan'),'k','LineWidth',2);end
                labs={'Raw dF/F','Forward (mm/s)','Abs. yaw (deg/s)','Heat command (V)'};ylabel(ax,labs{v});xlim(ax,[-2 3]);
                if v==1,ylim(ax,[min([lim(1);d(isfinite(d))]) max([lim(2);d(isfinite(d))])]);end
            end
            xline(ax,0,':r');if row==5,xlabel(ax,'Time from maximum deceleration (s)');end
        end
    end
    % Identical behavioral y limits across columns within each row.
    for row=3:5,axs=gobjects(1,4);for c=1:4,axs(c)=nexttile(tl,(row-1)*4+c);end;yl=cell2mat(get(axs,'YLim'));set(axs,'YLim',[min(yl(:,1)) max(yl(:,2))]);end
    title(tl,{strrep(name,'_',' '),'Behavior-selected slowdowns | turn: >=30 deg absolute yaw in +/-0.5 s | thin: events; black: mean'});
    saveplot(fig,out,stem);
end
fig=figure('Visible','off','Color','w','Position',[20 20 1200 500]);tl=tiledlayout(fig,1,2);colors=lines(4);
for v=1:2
    ax=nexttile(tl);hold(ax,'on');for c=1:4,sel=E.context==contexts(c);if v==1,xx=E.speedDrop(sel);else,xx=E.absoluteTurnDeg(sel);end;scatter(ax,xx,E.deltaDff(sel),30,colors(c,:),'filled','DisplayName',sprintf('%s (n=%d)',contexts(c),sum(sel)));end
    if v==1,xlabel(ax,'Speed drop (mm/s)');else,xlabel(ax,'Absolute turn in +/-0.5 s (deg)');end
    ylabel(ax,'dF/F change: 0:1 s minus -2:-1 s');yline(ax,0,':k','HandleVisibility','off');legend(ax,'Location','best');grid(ax,'on');
end
title(tl,{strrep(name,'_',' '),'Explore overlapping movement magnitudes; events are not matched'});saveplot(fig,out,'activity_vs_movement');
end
function saveplot(fig,out,stem)
if ~feature('ShowFigureWindows'),savefig_visible_on_open(fig,fullfile(out,[stem '.fig']));else,savefig(fig,fullfile(out,[stem '.fig']));end;exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',140);close(fig);
end
