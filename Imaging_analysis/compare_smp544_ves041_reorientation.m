function results = compare_smp544_ves041_reorientation(outputFolder)
%COMPARE_SMP544_VES041_REORIENTATION Event-aligned comparison of SMP544 and VES041.
%
% SMP544 analysis is restricted to six trials with sustained, varied
% walking. Behavioral events are divided into:
%   1. low-speed reorientation with slowing
%   2. low-speed reorientation without slowing
%   3. slowdown without reorientation
%   4. sustained walking
%
% Reorientation: >=15 deg net heading change around a rotational-speed
% peak and mean forward speed <=5 mm/s from 0 to +0.75 s. Slowing: forward
% speed falls by >=1 mm/s from -1:-0.25 s to 0:+0.75 s. Heat-active epochs
% and events within 3 s of them are omitted. Fluorescence is smoothed over
% 0.25 s, detrended with a 30 s moving median, robustly z-scored within
% each trial, and baseline-subtracted using -3:-2 s before each event.

if nargin<1 || strlength(string(outputFolder))==0
    outputFolder='/Users/sophiarenauld/stacks/smp544_behavior_quantification/reorientation_state_comparison';
end
if ~isfolder(outputFolder),mkdir(outputFolder);end

base='/Users/sophiarenauld/stacks';
smpTrials=["20260720-3_smp544_fly2_dark"; ...
    "20260720-4_smp544_fly2_dark2_volume"; ...
    "20260721-2_smp544_fly1_clbar"; ...
    "20260721-3_smp544_fly1_clbar"; ...
    "20260721-4_smp544_fly1_olheat"; ...
    "20260722-2_smp544_fly1_clbar"];
% All previously screened good-behavior, non-heat VES041 trials. Moving
% optomotor epochs and +/-2 s around detected bar jumps are removed below.
vesTrials=["20250812-20_ves041_fly3_clbarjumps"; ...
    "20250812-21_ves041_fly3_olgrating"; ...
    "20250812-22_ves041_fly3_bbarjumps2"; ...
    "20250813-3_ves041_fly1_olgrating"; ...
    "20250819-10_ves041_fly1_olgrating5"; ...
    "20251001-21_ves041_fly3_clbbar"; ...
    "20260213-2_ves041_clbar"; ...
    "20260213-4_ves041_fly2_clbar"; ...
    "20260217-2_ves041_vt_clbbar1"; ...
    "20260217-3_ves041_vt_clbbar2_actuallydark"; ...
    "20260217-4_ves041_vt_clbbar3"];

categories=["reorientation + slowing","reorientation, no slowing", ...
    "slowdown, no reorientation","sustained walking"];
colors=[.18 .42 .78;.55 .25 .72;.90 .48 .12;.20 .62 .35];
smp=cell(numel(smpTrials),1);rows=cell(0,8);
for i=1:numel(smpTrials)
    smp{i}=load_smp(fullfile(base,smpTrials(i)));
    smp{i}=analyze_trial(smp{i});
    for c=1:4
        rows(end+1,:)={smpTrials(i),categories(c),numel(smp{i}.events{c}), ... %#ok<AGROW>
            mean(smp{i}.eventEffect{c},'omitnan'),median(smp{i}.eventEffect{c},'omitnan'), ...
            mean(smp{i}.eventPreFwd{c},'omitnan'),mean(smp{i}.eventPostFwd{c},'omitnan'), ...
            mean(smp{i}.eventHeading{c},'omitnan')};
    end
    make_trial_figure(smp{i},smpTrials(i),categories,colors,outputFolder);
end
summary=cell2table(rows,'VariableNames',{'trial','category','nEvents','meanDffChangeZ', ...
    'medianDffChangeZ','meanPreFwd','meanPostFwd','meanHeadingChangeDeg'});
summary.trial=string(summary.trial);summary.category=string(summary.category);
writetable(summary,fullfile(outputFolder,'smp544_event_summary.csv'));
make_smp_summary(smp,smpTrials,categories,colors,outputFolder);

ves=cell(numel(vesTrials),1);
for i=1:numel(vesTrials)
    ves{i}=load_ves(fullfile(base,vesTrials(i)));
    ves{i}=analyze_trial(ves{i});
end
comparison=make_cell_comparison(smp,smpTrials,ves,vesTrials,outputFolder);
writetable(comparison,fullfile(outputFolder,'smp544_vs_ves041_reorientation_effects.csv'));
vesCounts=event_count_table(ves,vesTrials);
writetable(vesCounts,fullfile(outputFolder,'ves041_reorientation_event_counts.csv'));
slowdownTrialFolder=fullfile(outputFolder,'smp544_all_slowdowns_by_trial');
if ~isfolder(slowdownTrialFolder),mkdir(slowdownTrialFolder);end
for i=1:numel(smp)
    smp{i}=analyze_slowdowns(smp{i});
    make_slowdown_trial_figure(smp{i},smpTrials(i),slowdownTrialFolder);
end
for i=1:numel(ves),ves{i}=analyze_slowdowns(ves{i});end
slowdownTable=make_slowdown_comparison(smp,smpTrials,ves,vesTrials,outputFolder);
writetable(slowdownTable,fullfile(outputFolder,'smp544_vs_ves041_all_slowdowns.csv'));
results=struct('smpSummary',summary,'comparison',comparison,'outputFolder',string(outputFolder));
results.slowdownComparison=slowdownTable;
disp(summary);disp(comparison);disp(vesCounts);disp(slowdownTable);
end

function D=load_smp(folder)
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));assert(~isempty(f),'No SMP544 input in %s',folder);
q=load(fullfile(f(1).folder,f(1).name),'dq','roi');
t=q.dq(1).t(:);fwd=q.dq(1).bvf(:);rot=rad2deg(q.dq(1).bvy(:));heading=q.dq(1).bh(:);
d=q.roi.ts;if isvector(d),d=d(:)';end
n=min([numel(t),numel(fwd),numel(rot),numel(heading),size(d,2)]);
% Bilateral mean after separately normalizing each ROI prevents the brighter
% side from dominating while retaining shared SMP544 state modulation.
dff=zeros(n,1);for r=1:size(d,1),dff=dff+normalize_dff(d(r,1:n)',t(1:n));end;dff=dff/size(d,1);
stim=false(n,1);
if isfield(q.dq(1),'heat_active') && ~isempty(q.dq(1).heat_active),stim=logical(q.dq(1).heat_active(1:n));
elseif isfield(q.dq(1),'heat') && ~isempty(q.dq(1).heat),h=q.dq(1).heat(1:n);stim=h(:)>median(h,'omitnan')+2*mad(h,1);end
D=pack(t(1:n),dff,fwd(1:n),rot(1:n),heading(1:n),stim);
end

function D=load_ves(folder)
q=load(fullfile(folder,'smr_analysis','fly_data.mat'),'fly_data');x=q.fly_data;
t=x.dq(1).t(:);fwd=x.dq(1).bvf(:);rot=x.dq(1).bvy_deg(:);heading=x.dq(1).bh(:);d=x.roi(1).ts(:);
n=min([numel(t),numel(fwd),numel(rot),numel(heading),numel(d)]);stim=false(n,1);
if contains(lower(folder),'grating') && isfield(x.dq(1),'vvy'),v=x.dq(1).vvy(:);stim(1:min(n,numel(v)))=abs(v(1:min(n,numel(v))))>.02;end
if isfield(x,'jump_detected'),j=find(logical(x.jump_detected(:)));j=j(j<=n);fs=1/median(diff(t));for k=j(:)',stim(max(1,k-round(2*fs)):min(n,k+round(2*fs)))=true;end,end
D=pack(t(1:n),normalize_dff(d(1:n),t(1:n)),fwd(1:n),rot(1:n),heading(1:n),stim);
end

function D=pack(t,dff,fwd,rot,heading,stim)
D=struct('t',t(:),'dff',dff(:),'fwd',fwd(:),'rot',rot(:),'heading',heading(:),'stim',logical(stim(:)));
end

function z=normalize_dff(d,t)
fs=1/median(diff(t));sm=smoothdata(d(:),'gaussian',max(3,round(.25*fs)));
res=sm-movmedian(sm,max(3,round(30*fs)),'omitmissing');s=1.4826*mad(res,1);if ~isfinite(s)||s<=0,s=std(res,'omitnan');end;z=res/s;
end

function D=analyze_trial(D)
fs=1/median(diff(D.t));n=numel(D.t);pre=round(3*fs);post=round(3*fs);rel=(-pre:post)'/fs;
f=smoothdata(D.fwd,'gaussian',max(3,round(.2*fs)));r=smoothdata(D.rot,'gaussian',max(3,round(.12*fs)));
stimBuffered=movmax(double(D.stim),[pre post])>0;edge=false(n,1);edge(1:pre)=true;edge(end-post+1:end)=true;forbidden=stimBuffered|edge;
prom=max(20,2*1.4826*mad(abs(r)-movmedian(abs(r),max(3,round(5*fs))),1));
[~,turns]=findpeaks(abs(r),'MinPeakDistance',max(1,round(1.5*fs)),'MinPeakProminence',prom,'MinPeakHeight',20);
turns=turns(~forbidden(turns));
preA=-round(fs):-round(.25*fs);postA=0:round(.75*fs);headA=-round(.25*fs):round(.75*fs);
[turnSlow,turnNoSlow]=deal([]);
turnHeading=[];turnPre=[];turnPost=[];
for p=turns(:)'
    preF=mean(f(p+preA),'omitnan');postF=mean(f(p+postA),'omitnan');
    hd=heading_change(D.heading,p,headA);
    if postF<=5 && hd>=15
        if preF-postF>=1,turnSlow(end+1)=p;else,turnNoSlow(end+1)=p;end %#ok<AGROW>
        turnHeading(end+1)=hd;turnPre(end+1)=preF;turnPost(end+1)=postF; %#ok<AGROW>
    end
end

acc=[0;diff(f)]*fs;accProm=max(.5,2*1.4826*mad(acc,1));[~,slow]=findpeaks(-acc,'MinPeakDistance',round(1.5*fs),'MinPeakProminence',accProm);
slow=slow(~forbidden(slow));slowOnly=[];
for p=slow(:)'
    preF=mean(f(p+preA),'omitnan');postF=mean(f(p+postA),'omitnan');
    hd=heading_change(D.heading,p,headA);
    if preF-postF>=1 && hd<15,slowOnly(end+1)=p;end %#ok<AGROW>
end
walk=[];step=max(1,round(1.5*fs));
for p=pre+1:step:n-post
    if forbidden(p),continue,end
    preF=mean(f(p+preA),'omitnan');postF=mean(f(p+postA),'omitnan');
    hd=heading_change(D.heading,p,headA);
    if preF>=2 && postF>=2 && abs(preF-postF)<1 && hd<5 && mean(abs(r(p+postA)),'omitnan')<30,walk(end+1)=p;end %#ok<AGROW>
end
D.events={turnSlow(:),turnNoSlow(:),slowOnly(:),walk(:)};D.rel=rel;D.traces=cell(4,3);D.eventEffect=cell(4,1);D.eventPreFwd=cell(4,1);D.eventPostFwd=cell(4,1);D.eventHeading=cell(4,1);
for c=1:4
    e=D.events{c};D.traces{c,1}=extract(D.dff,e,pre,post);D.traces{c,2}=extract(f,e,pre,post);D.traces{c,3}=extract(abs(r),e,pre,post);
    if ~isempty(e)
        baseMask=rel>=-3&rel<=-2;eventMask=rel>=.25&rel<=.75;
        X=D.traces{c,1};X=X-mean(X(baseMask,:),1,'omitnan');D.traces{c,1}=X;D.eventEffect{c}=mean(X(eventMask,:),1,'omitnan')';
        D.eventPreFwd{c}=arrayfun(@(p)mean(f(p+preA),'omitnan'),e);D.eventPostFwd{c}=arrayfun(@(p)mean(f(p+postA),'omitnan'),e);
        D.eventHeading{c}=arrayfun(@(p)heading_change(D.heading,p,headA),e);
    else,D.eventEffect{c}=[];D.eventPreFwd{c}=[];D.eventPostFwd{c}=[];D.eventHeading{c}=[];end
end
end

function D=analyze_slowdowns(D)
% Find all forward-speed drops, independent of whether a turn occurs.
fs=1/median(diff(D.t));n=numel(D.t);pre=round(3*fs);post=round(3*fs);rel=(-pre:post)'/fs;
f=smoothdata(D.fwd,'gaussian',max(3,round(.2*fs)));r=smoothdata(D.rot,'gaussian',max(3,round(.12*fs)));
stimBuffered=movmax(double(D.stim),[pre post])>0;edge=false(n,1);edge(1:pre)=true;edge(end-post+1:end)=true;forbidden=stimBuffered|edge;
acc=[0;diff(f)]*fs;accProm=max(.5,2*1.4826*mad(acc,1));
[~,candidate]=findpeaks(-acc,'MinPeakDistance',round(1.5*fs),'MinPeakProminence',accProm);
candidate=candidate(~forbidden(candidate));preA=-round(fs):-round(.25*fs);postA=0:round(.75*fs);headA=-round(.25*fs):round(.75*fs);
events=[];hasTurn=[];preF=[];postF=[];hd=[];
for p=candidate(:)'
    a=mean(f(p+preA),'omitnan');b=mean(f(p+postA),'omitnan');h=heading_change(D.heading,p,headA);
    if a-b>=1
        events(end+1,1)=p;hasTurn(end+1,1)=h>=15;preF(end+1,1)=a;postF(end+1,1)=b;hd(end+1,1)=h; %#ok<AGROW>
    end
end
D.slowdown=struct('events',events,'hasTurn',logical(hasTurn),'preFwd',preF,'postFwd',postF,'headingChange',hd,'rel',rel);
D.slowdown.traces={extract(D.dff,events,pre,post),extract(f,events,pre,post),extract(abs(r),events,pre,post)};
if ~isempty(events)
    baseMask=rel>=-3&rel<=-2;X=D.slowdown.traces{1};X=X-mean(X(baseMask,:),1,'omitnan');D.slowdown.traces{1}=X;
    eventMask=rel>=.25&rel<=.75;D.slowdown.effect=mean(X(eventMask,:),1,'omitnan')';
else,D.slowdown.effect=[];end
end

function X=extract(x,e,pre,post)
X=nan(pre+post+1,numel(e));for i=1:numel(e),X(:,i)=x(e(i)-pre:e(i)+post);end
end

function hd=heading_change(heading,p,offsets)
delta=heading(p+offsets(end))-heading(p+offsets(1));
hd=abs(rad2deg(atan2(sin(delta),cos(delta))));
end

function make_trial_figure(D,name,cats,colors,out)
fig=figure('Visible','off','Color','w','Position',[50 50 1250 850]);tl=tiledlayout(3,1,'TileSpacing','compact');ylabs={'Baseline-adjusted dF/F (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
for v=1:3,ax=nexttile(tl);hold(ax,'on');for c=1:4,X=D.traces{c,v};if isempty(X),continue,end;plot(ax,D.rel,mean(X,2,'omitnan'),'Color',colors(c,:),'LineWidth',2,'DisplayName',sprintf('%s (n=%d)',cats(c),size(X,2)));end;xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-3 3]);ylabel(ax,ylabs{v});grid(ax,'on');if v==1,legend(ax,'Location','best');end;if v==3,xlabel(ax,'Time from event (s)');end;end
title(tl,strrep(name,'_','\_'));savefig(fig,fullfile(out,name+'_event_classes.fig'));exportgraphics(fig,fullfile(out,name+'_event_classes.png'),'Resolution',220);close(fig);
end

function make_smp_summary(T,names,cats,colors,out)
fig=figure('Visible','off','Color','w','Position',[50 50 1450 900]);tl=tiledlayout(2,2,'TileSpacing','compact');labels={'Baseline-adjusted dF/F (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
for v=1:3,ax=nexttile(tl);hold(ax,'on');for c=1:4,[m,se,t,n]=trial_weighted(T,c,v);if isempty(m),continue,end;shade(ax,t,m,se,colors(c,:));plot(ax,t,m,'Color',colors(c,:),'LineWidth',2.4,'DisplayName',sprintf('%s (%d trials)',cats(c),n));end;xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-3 3]);ylabel(ax,labels{v});grid(ax,'on');if v==1,legend(ax,'Location','best');end;if v>1,xlabel(ax,'Time from event (s)');end;end
ax=nexttile(tl);hold(ax,'on');for c=1:4,vals=[];for i=1:numel(T),if ~isempty(T{i}.eventEffect{c}),vals(end+1)=mean(T{i}.eventEffect{c},'omitnan');end,end;scatter(ax,c*ones(size(vals))+.08*randn(size(vals)),vals,55,colors(c,:),'filled');if ~isempty(vals),plot(ax,[c-.2 c+.2],repmat(mean(vals),1,2),'Color',colors(c,:),'LineWidth',3);end;end;yline(ax,0,':k');xticks(ax,1:4);xticklabels(ax,cats);xtickangle(ax,20);ylabel(ax,'Mean +0.25:+0.75 s dF/F change');grid(ax,'on');title(tl,'SMP544 around behavioral transitions: high-quality trials');savefig(fig,fullfile(out,'smp544_reorientation_state_summary.fig'));exportgraphics(fig,fullfile(out,'smp544_reorientation_state_summary.png'),'Resolution',240);close(fig);
end

function tableOut=make_cell_comparison(S,snames,V,vnames,out)
[sm,seS,t,nS]=trial_weighted_reorientation(S,1);[vm,seV,~,nV]=trial_weighted_reorientation(V,1);
[sf,seSf]=trial_weighted_reorientation(S,2);[vf,seVf]=trial_weighted_reorientation(V,2);
[sr,seSr]=trial_weighted_reorientation(S,3);[vr,seVr]=trial_weighted_reorientation(V,3);
fig=figure('Visible','off','Color','w','Position',[60 60 1350 850]);tl=tiledlayout(2,2,'TileSpacing','compact');C={[.15 .42 .78],[.82 .18 .18]};names={'SMP544','VES041'};Y={sm,vm;sf,vf;sr,vr};E={seS,seV;seSf,seVf;seSr,seVr};labs={'Baseline-adjusted dF/F (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
for p=1:3,ax=nexttile(tl);hold(ax,'on');for g=1:2,shade(ax,t,Y{p,g},E{p,g},C{g});plot(ax,t,Y{p,g},'Color',C{g},'LineWidth',2.5,'DisplayName',names{g});end;xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-3 3]);ylabel(ax,labs{p});xlabel(ax,'Time from reorientation (s)');grid(ax,'on');if p==1,legend(ax,'Location','best');end;end
ax=nexttile(tl);hold(ax,'on');[sv,sUsed,sCount]=trial_effect(S,snames);[vv,vUsed,vCount]=trial_effect(V,vnames);scatter(ax,ones(size(sv))+.06*randn(size(sv)),sv,60,C{1},'filled');scatter(ax,2*ones(size(vv))+.06*randn(size(vv)),vv,60,C{2},'filled');plot(ax,[.8 1.2],repmat(mean(sv),1,2),'Color',C{1},'LineWidth',3);plot(ax,[1.8 2.2],repmat(mean(vv),1,2),'Color',C{2},'LineWidth',3);yline(ax,0,':k');xticks(ax,[1 2]);xticklabels(ax,names);ylabel(ax,'Trial-mean +0.25:+0.75 s dF/F change');grid(ax,'on');title(tl,sprintf('Matched low-speed reorientations: SMP544 (n=%d trials) versus VES041 (n=%d)',nS,nV));savefig(fig,fullfile(out,'smp544_vs_ves041_reorientation.fig'));exportgraphics(fig,fullfile(out,'smp544_vs_ves041_reorientation.png'),'Resolution',240);close(fig);
tableOut=table([repmat("SMP544",numel(sv),1);repmat("VES041",numel(vv),1)],[sUsed;vUsed],[sCount;vCount],[sv(:);vv(:)],'VariableNames',{'cellType','trial','nReorientations','meanDffChangeZ'});
end

function [m,se,t,n]=trial_weighted(T,c,v)
t=(-3:.05:3)';M=[];for i=1:numel(T),X=T{i}.traces{c,v};if isempty(X),continue,end;y=mean(X,2,'omitnan');M(:,end+1)=interp1(T{i}.rel,y,t,'linear',NaN);end %#ok<AGROW>
n=size(M,2);if n==0,m=[];se=[];return,end;m=mean(M,2,'omitnan');se=std(M,0,2,'omitnan')/sqrt(n);
end
function [m,se,t,n]=trial_weighted_reorientation(T,v)
t=(-3:.05:3)';M=[];for i=1:numel(T),X=[T{i}.traces{1,v},T{i}.traces{2,v}];if isempty(X),continue,end;y=mean(X,2,'omitnan');M(:,end+1)=interp1(T{i}.rel,y,t,'linear',NaN);end %#ok<AGROW>
n=size(M,2);m=mean(M,2,'omitnan');se=std(M,0,2,'omitnan')/sqrt(n);
end
function [vals,usedNames,counts]=trial_effect(T,names)
vals=[];usedNames=strings(0,1);counts=[];
for i=1:numel(T)
    x=[T{i}.eventEffect{1};T{i}.eventEffect{2}];
    if ~isempty(x),vals(end+1)=mean(x,'omitnan');usedNames(end+1,1)=names(i);counts(end+1,1)=numel(x);end %#ok<AGROW>
end
end
function C=event_count_table(T,names)
nSlow=zeros(numel(T),1);nNoSlow=nSlow;
for i=1:numel(T),nSlow(i)=numel(T{i}.events{1});nNoSlow(i)=numel(T{i}.events{2});end
C=table(names,nSlow,nNoSlow,nSlow+nNoSlow,'VariableNames',{'trial','reorientationWithSlowing','reorientationWithoutSlowing','totalReorientations'});
end

function tableOut=make_slowdown_comparison(S,sNames,V,vNames,out)
% Pooled traces are means of trial means, giving every recording equal weight.
[sd,seSd,t,nS]=slowdown_weighted(S,1,[]);[vd,seVd,~,nV]=slowdown_weighted(V,1,[]);
[sf,seSf]=slowdown_weighted(S,2,[]);[vf,seVf]=slowdown_weighted(V,2,[]);
[sr,seSr]=slowdown_weighted(S,3,[]);[vr,seVr]=slowdown_weighted(V,3,[]);
[sTurn,seST]=slowdown_weighted(S,1,true);[sNoTurn,seSN]=slowdown_weighted(S,1,false);
[vTurn,seVT]=slowdown_weighted(V,1,true);[vNoTurn,seVN]=slowdown_weighted(V,1,false);
C1=[.15 .42 .78];C2=[.82 .18 .18];
fig=figure('Visible','off','Color','w','Position',[60 60 1450 900]);tl=tiledlayout(2,2,'TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');shade(ax,t,sd,seSd,C1);shade(ax,t,vd,seVd,C2);plot(ax,t,sd,'Color',C1,'LineWidth',2.5,'DisplayName','SMP544');plot(ax,t,vd,'Color',C2,'LineWidth',2.5,'DisplayName','VES041');xline(ax,0,':k','HandleVisibility','off');ylabel(ax,'Baseline-adjusted dF/F (robust z)');title(ax,'All slowdowns');legend(ax,'Location','best');grid(ax,'on');xlim(ax,[-3 3]);
ax=nexttile(tl);hold(ax,'on');shade(ax,t,sf,seSf,C1);shade(ax,t,vf,seVf,C2);plot(ax,t,sf,'Color',C1,'LineWidth',2.5);plot(ax,t,vf,'Color',C2,'LineWidth',2.5);xline(ax,0,':k','HandleVisibility','off');ylabel(ax,'Forward speed (mm/s)');title(ax,'Matched forward-speed transition');grid(ax,'on');xlim(ax,[-3 3]);
ax=nexttile(tl);hold(ax,'on');shade(ax,t,sr,seSr,C1);shade(ax,t,vr,seVr,C2);plot(ax,t,sr,'Color',C1,'LineWidth',2.5);plot(ax,t,vr,'Color',C2,'LineWidth',2.5);xline(ax,0,':k','HandleVisibility','off');ylabel(ax,'Rotational speed (deg/s)');xlabel(ax,'Time from maximum deceleration (s)');title(ax,'Turning was not required');grid(ax,'on');xlim(ax,[-3 3]);
ax=nexttile(tl);hold(ax,'on');plot(ax,t,sTurn,'Color',C1,'LineWidth',2.5,'DisplayName','SMP544: slowdown + turn');plot(ax,t,sNoTurn,'--','Color',C1,'LineWidth',2.2,'DisplayName','SMP544: slowdown only');plot(ax,t,vTurn,'Color',C2,'LineWidth',2.5,'DisplayName','VES041: slowdown + turn');plot(ax,t,vNoTurn,'--','Color',C2,'LineWidth',2.2,'DisplayName','VES041: slowdown only');xline(ax,0,':k','HandleVisibility','off');ylabel(ax,'Baseline-adjusted dF/F (robust z)');xlabel(ax,'Time from maximum deceleration (s)');title(ax,'Secondary split by heading change');legend(ax,'Location','best');grid(ax,'on');xlim(ax,[-3 3]);
title(tl,sprintf('All spontaneous slowdowns: SMP544 (%d trials) versus VES041 (%d trials)',nS,nV));savefig(fig,fullfile(out,'smp544_vs_ves041_all_slowdowns.fig'));exportgraphics(fig,fullfile(out,'smp544_vs_ves041_all_slowdowns.png'),'Resolution',240);close(fig);

rows=cell(0,8);
for g=1:2
    if g==1,T=S;names=sNames;cellName="SMP544";else,T=V;names=vNames;cellName="VES041";end
    for i=1:numel(T)
        n=numel(T{i}.slowdown.events);nt=sum(T{i}.slowdown.hasTurn);effect=mean(T{i}.slowdown.effect,'omitnan');
        rows(end+1,:)={cellName,names(i),n,nt,n-nt,effect,mean(T{i}.slowdown.preFwd,'omitnan'),mean(T{i}.slowdown.postFwd,'omitnan')}; %#ok<AGROW>
    end
end
tableOut=cell2table(rows,'VariableNames',{'cellType','trial','nSlowdowns','nWithTurn','nWithoutTurn','meanDffChangeZ','meanPreFwd','meanPostFwd'});tableOut.cellType=string(tableOut.cellType);tableOut.trial=string(tableOut.trial);
end

function make_slowdown_trial_figure(D,name,out)
% Split slowdowns by heading change; plot every event and each group mean.
labels={'Baseline-adjusted dF/F (robust z)','Forward speed (mm/s)', ...
    'Rotational speed (deg/s)'};
groupNames={'With reorientation (heading change >=15 deg)', ...
    'Without reorientation (heading change <15 deg)'};
groupColors={[.86 .32 .12],[.05 .30 .72]};
groupMasks={D.slowdown.hasTurn,~D.slowdown.hasTurn};
fig=figure('Visible','off','Color','w','Position',[40 40 1550 950]);
tl=tiledlayout(3,2,'TileSpacing','compact','Padding','compact');axs=gobjects(3,2);
for v=1:3
    for g=1:2
        ax=nexttile(tl);axs(v,g)=ax;hold(ax,'on');X=D.slowdown.traces{v};
        if ~isempty(X),X=X(:,groupMasks{g});end
        if ~isempty(X)
            plot(ax,D.slowdown.rel,X,'Color',[.80 .80 .80],'LineWidth',.65, ...
                'HandleVisibility','off');
            plot(ax,D.slowdown.rel,mean(X,2,'omitnan'),'Color',groupColors{g}, ...
                'LineWidth',3,'DisplayName',sprintf('Mean (n=%d)',size(X,2)));
        end
        xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-3 3]);grid(ax,'on');
        if g==1,ylabel(ax,labels{v});end
        if v==1,title(ax,groupNames{g});if ~isempty(X),legend(ax,'Location','best');end,end
        if v==3,xlabel(ax,'Time from maximum deceleration (s)');end
    end
    linkaxes(axs(v,:),'y');
end
title(tl,sprintf('%s | spontaneous slowdowns split by reorientation',strrep(name,'_','\_')));
savefig(fig,fullfile(out,name+'_all_slowdowns_individual_traces.fig'));
exportgraphics(fig,fullfile(out,name+'_all_slowdowns_individual_traces.png'),'Resolution',240);
close(fig);
end

function [m,se,t,n]=slowdown_weighted(T,v,turnState)
t=(-3:.05:3)';M=[];
for i=1:numel(T)
    X=T{i}.slowdown.traces{v};if isempty(X),continue,end
    if ~isempty(turnState),X=X(:,T{i}.slowdown.hasTurn==turnState);end
    if isempty(X),continue,end
    y=mean(X,2,'omitnan');M(:,end+1)=interp1(T{i}.slowdown.rel,y,t,'linear',NaN); %#ok<AGROW>
end
n=size(M,2);if n==0,m=nan(size(t));se=m;return,end;m=mean(M,2,'omitnan');se=std(M,0,2,'omitnan')/sqrt(n);
end
function shade(ax,x,m,se,c)
fill(ax,[x;flipud(x)],[m-se;flipud(m+se)],c,'FaceAlpha',.16,'EdgeColor','none','HandleVisibility','off');
end
