function report=plot_ves041_barrier_slowdown_heat_bins()
% Behavior-selected slowdowns, pooling turning and nonturning events.
% Control identity comes from the protocol, NOT from low measured voltage.
% Heat bins use median measured command in +/-0.5 s: [0,1), [1,2.5), >=2.5 V.
% Near-zero negative monitor noise is assigned to zero for binning only.
% These are command-voltage bins, not measured temperatures.
% Retain each recording's saved manual ROI normalization; never pool flies
% or changing geometries. No fluorescence-based event selection or clipping.
old=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
out=fullfile(imaging_analysis_artifact_root(),'diagnostics','barrier_slowdown_heat_bins_20261027');
assert(~isfolder(out),'Output exists; choose a fresh output name to preserve prior figures.');mkdir(out);
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat'); ...
    dir('/Users/sophiarenauld/stacks/2026-10-20/**/*_or_z0z1_s_.mat'); ...
    dir('/Users/sophiarenauld/stacks/2026-10-24/registration_xyz_fixed_20261026/a2p_ready/**/*_or_z0z1_s_.mat')];
rows=cell(0,6);allE=table;groupRows=cell(0,9);
for k=1:numel(files)
    id=string(erase(files(k).name,'_or_z0z1_s_.mat'));
    try
        a=load(fullfile(files(k).folder,files(k).name),'dq','roi');
        a=import_ves041_protocol(a,files(k).folder);
        if ~isfield(a,'protocol')||~ismember(a.protocol.kind,{'heat-barrier','heat-pyramid'}),continue,end
        if size(a.roi.mask,1)==size(a.roi.mask,2)
            rows(end+1,:)={id,string(a.protocol.kind),"skipped square",0,"",""};continue %#ok<AGROW>
        end
        assert(a.protocol.geometryAligned,'Geometry not aligned');
        assert(isequaln(a.roi.ts,a.roi.pxfilter.manualts),'Expected manual ROI');
        dest=fullfile(out,id);mkdir(dest);
        [E,X,rel,signals]=selectEvents(a,id);
        writetable(E,fullfile(dest,'slowdown_events.csv'));allE=[allE;E]; %#ok<AGROW>
        geometries=unique(E.geometry(E.included),'stable');
        if isempty(geometries),geometries="no eligible events";end
        for geometry=geometries'
            q=E.included&E.geometry==geometry;
            stem=matlab.lang.makeValidName(char(geometry));
            drawSegments(E(q,:),X(:,q,:),rel,dest,stem,id,a.roi.opt.nrm.nrmstr,geometry);
            for g=1:4
                take=q&E.group==g;n=sum(take);
                groupRows(end+1,:)={id,geometry,g,n,sum(E.reorientation(take)), ...
                    mean(E.medianHeatV(take),'omitnan'),mean(E.deltaDff(take),'omitnan'), ...
                    mean(E.speedBefore(take),'omitnan'),mean(E.speedDrop(take),'omitnan')}; %#ok<AGROW>
            end
        end
        drawTimeline(E,signals,dest,id);
        settings=struct('roiSource',a.analysisRoiSource,'normalization',a.roi.opt.nrm.nrmstr, ...
            'alignment','maximum deceleration within a contiguous candidate run', ...
            'minimumPreSpeed',1,'minimumSpeedDrop',1,'accelerationThreshold',-2, ...
            'smoothingSec',.25,'minimumEventSpacingSec',2,'windowSec',[-2 3], ...
            'heatBinWindowSec',[-.5 .5],'heatBinEdgesV',[0 1 2.5 Inf], ...
            'controlDefinition','control_trial==1 during active placed barrier', ...
            'fullWindowSameTrialAndGeometry',true,'turnsPooled',true);
        save(fullfile(dest,'segmented_events.mat'),'E','X','rel','settings','-v7.3');
        rows(end+1,:)={id,string(a.protocol.kind),"success",sum(E.included),string(dest),""}; %#ok<AGROW>
        fprintf('%s: %d accepted slowdowns; %d geometry groups\n',id,sum(E.included),numel(geometries));
    catch ME
        rows(end+1,:)={id,"unknown","error",0,"",string(getReport(ME,'extended','hyperlinks','off'))}; %#ok<AGROW>
        disp(rows(end,:));
    end
end
report=cell2table(rows,'VariableNames',{'recording','protocol','status','nIncluded','output','error'});
writetable(report,fullfile(out,'run_report.csv'));
if ~isempty(allE),writetable(allE,fullfile(out,'all_slowdown_events.csv'));end
G=cell2table(groupRows,'VariableNames',{'recording','geometry','group','nEvents','nWithTurn','meanMedianHeatV','meanDeltaDff','meanPreSpeed','meanSpeedDrop'});
writetable(G,fullfile(out,'group_summary.csv'));disp(G);
assert(~any(report.status=="error"),'Some recordings failed; inspect run report.');
end

function [E,X,rel,S]=selectEvents(a,id)
t=a.dq(1).t(:);d=double(a.roi.ts(:));rawF=a.dq(1).bvf(:);rawR=rad2deg(a.dq(1).bvy(:));h=a.dq(1).heat(:);p=a.protocol;
side=a.dq(1).bvs(:);fs=1/median(diff(t));smoothN=max(3,round(.25*fs));
f=smoothdata(rawF,'movmedian',smoothN);r=smoothdata(rawR,'movmedian',smoothN);
valid=isfinite(d)&isfinite(rawF)&isfinite(rawR)&isfinite(h)&isfinite(p.distance_to_barrier_mm(:))&logical(a.dq(1).bergii.behaviorVolumeValid(:));
acc=gradient(f,t);cand=valid&acc<=-2;
beg=find(diff([false;cand])==1);ends=find(diff([cand;false])==-1);idx=zeros(size(beg));
for j=1:numel(beg),[~,m]=min(acc(beg(j):ends(j)));idx(j)=beg(j)+m-1;end
geometry=repmat(string(p.kind),size(t));
if isfield(p,'heat_geometry'),geometry=string(p.heat_geometry(:));end
distance=p.distance_to_barrier_mm(:);
virtual=nan(size(t));if isfield(p,'virtual_heat_voltage'),virtual=p.virtual_heat_voltage(:);end
rel=(-2:.05:3)';rows=cell(0,15);X=nan(numel(rel),0,8);last=-Inf;
for q=idx'
    et=t(q);if et-last<2||et<t(1)+2||et>t(end)-3,continue,end
    win=t>=et-2&t<=et+3;pre=t>=et-.75&t<et-.25;post=t>=et+.25&t<=et+.75;core=t>=et-.5&t<=et+.5;
    if ~all(valid(win))||any(diff(t(win))>2/fs),continue,end
    before=median(f(pre));drop=before-median(f(post));if before<1||drop<1,continue,end
    last=et;turn=trapz(t(core),abs(r(core)));voltage=median(h(core));group=0;reason="";
    sameTrial=all(isfinite(p.trial_number(win)))&&numel(unique(p.trial_number(win)))==1;
    sameGeometry=all(~ismissing(geometry(win)))&&numel(unique(geometry(win)))==1;
    placed=all(p.barrier_active(win)>0);control=p.control_trial(win);
    if ~sameTrial||~sameGeometry||~placed||any(~isfinite(control))||numel(unique(control))~=1
        reason="window crosses unplaced barrier, trial, geometry, or control boundary";
    elseif control(1)==1
        group=1;
    elseif control(1)~=0
        reason="unknown control flag";
    elseif voltage<-.05
        reason="unexpected negative measured command";
    elseif max(0,voltage)<1
        group=2;
    elseif voltage<2.5
        group=3;
    else
        group=4;
    end
    b=t>=et-2&t<et-1;response=t>=et&t<=et+1;base=mean(d(b));delta=mean(d(response))-base;
    geo=geometry(q);if ismissing(geo),geo="unplaced";end
    rows(end+1,:)={id,et,geo,p.trial_number(q),group,group>0,reason,voltage, ...
        before,drop,turn,turn>=30,base,mean(d(response)),delta}; %#ok<AGROW>
    vars={d,d-base,rawF,abs(rawR),side,h,virtual,distance};
    event=size(rows,1);
    X(:,event,:)=NaN;
    for v=1:8,X(:,event,v)=interp1(t,vars{v},et+rel,'linear',NaN);end
end
E=cell2table(rows,'VariableNames',{'recording','eventTime','geometry','barrierTrial','group','included','exclusionReason', ...
    'medianHeatV','speedBefore','speedDrop','absoluteTurnDeg','reorientation','baselineDff','responseDff','deltaDff'});
S=struct('t',t,'d',d,'f',rawF,'r',abs(rawR),'h',h);
end

function drawSegments(E,X,rel,out,stem,id,nrm,geometry)
colors=[.25 .45 .65;.75 .60 .12;.95 .45 .10;.78 .12 .12];
names={'Control barrier','Heat trial: 0 to <1 V','Heat trial: 1 to <2.5 V','Heat trial: >=2.5 V'};
labs={'Raw dF/F','dF/F minus -2:-1 s','Forward speed (mm/s)','Rotational speed (deg/s)', ...
    'Signed side speed (mm/s)','Measured heat command (V)','Virtual heat profile (V)','Signed barrier/peak distance (mm)'};
fig=figure('Visible','off','Color','w','Position',[20 20 2100 2100]);tl=tiledlayout(fig,8,4,'TileSpacing','compact','Padding','loose');
for v=1:8
    values=X(:,:,v);finite=values(isfinite(values));lim=[0 1];
    if ~isempty(finite),lim=[min(finite) max(finite)];if diff(lim)==0,lim=lim+[-.1 .1];end;lim=lim+[-1 1]*.05*diff(lim);end
    for g=1:4
        ax=nexttile(tl,(v-1)*4+g);hold(ax,'on');take=E.group==g;
        if any(take)&&any(isfinite(X(:,take,v)),'all')
            plot(ax,rel,X(:,take,v),'Color',.8+[0 0 0]+.2*colors(g,:),'LineWidth',.5);
            plot(ax,rel,mean(X(:,take,v),2,'omitnan'),'Color',colors(g,:),'LineWidth',2.3);
        else,text(ax,.5,.5,'No eligible data','Units','normalized','HorizontalAlignment','center');end
        xline(ax,0,':k');xlim(ax,[-2 3]);ylim(ax,lim);grid(ax,'on');ylabel(ax,labs{v});
        if v==1,title(ax,{names{g},sprintf('n=%d | with turn %d; without %d',sum(take),sum(E.reorientation(take)),sum(~E.reorientation(take)))});end
        if v==8,xlabel(ax,'Time from strongest deceleration (s)');end
    end
end
title(tl,{strrep(id,'_',' ')+" | "+geometry,"Full hand-drawn ROI | saved normalization: "+string(nrm)+" | thin: individual events; thick: mean"},'Interpreter','none');
xlabel(tl,{'Turn and no-turn slowdowns pooled. Full -2:+3 s windows must be valid, within one placed barrier trial and geometry.', ...
    'Bins: median measured command within +/-0.5 s, not temperature. Heat-trial zero voltage stays in low heat. Behavior-selected; no dF/F threshold.'},'FontSize',11);
savefig_visible_on_open(fig,fullfile(out,[stem '_slowdown_segments.fig']));exportgraphics(fig,fullfile(out,[stem '_slowdown_segments.png']),'Resolution',145);close(fig);
end

function drawTimeline(E,S,out,id)
fig=figure('Visible','off','Color','w','Position',[30 30 1700 1000]);tl=tiledlayout(fig,4,1,'TileSpacing','compact');
colors=[.25 .45 .65;.75 .60 .12;.95 .45 .10;.78 .12 .12];vars={S.d,S.f,S.r,S.h};labs={'dF/F','Forward (mm/s)','Rotational speed (deg/s)','Measured heat command (V)'};
axs=gobjects(4,1);
for v=1:4
    ax=nexttile(tl);axs(v)=ax;plot(ax,S.t,vars{v},'k');hold(ax,'on');
    for g=1:4
        for et=E.eventTime(E.included&E.group==g)',xline(ax,et,'Color',colors(g,:),'Alpha',.5);end
    end
    ylabel(ax,labs{v});grid(ax,'on');xlim(ax,[S.t(1) S.t(end)]);
end
linkaxes(axs,'x');xlabel(tl,'Trial time (s)');title(tl,{strrep(id,'_',' '),'Accepted slowdowns: blue control | gold low | orange medium | red high'});
savefig_visible_on_open(fig,fullfile(out,'slowdown_selection_timeline.fig'));exportgraphics(fig,fullfile(out,'slowdown_selection_timeline.png'),'Resolution',145);close(fig);
end
