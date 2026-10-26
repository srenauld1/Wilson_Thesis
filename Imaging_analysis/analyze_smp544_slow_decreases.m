function events = analyze_smp544_slow_decreases(outputFolder)
%ANALYZE_SMP544_SLOW_DECREASES Screen gradual SMP544 declines and behavior.
% A primary event is a local minimum in fluorescence smoothed over 0.75 s,
% separated from other minima by >=3 s, for which mean fluorescence falls
% by >=0.75 robust z from -2:-1 s to -0.25:+0.5 s. Stimulus epochs and a
% 3 s buffer are excluded. The event table reports forward-speed change,
% post-event speed, absolute rotational speed, and net heading displacement.

if nargin<1 || strlength(string(outputFolder))==0
    outputFolder='/Users/sophiarenauld/stacks/smp544_behavior_quantification/slow_activity_decreases';
end
if ~isfolder(outputFolder),mkdir(outputFolder);end
base='/Users/sophiarenauld/stacks';
trials=["20260720-3_smp544_fly2_dark";"20260720-4_smp544_fly2_dark2_volume"; ...
    "20260721-2_smp544_fly1_clbar";"20260721-3_smp544_fly1_clbar"; ...
    "20260721-4_smp544_fly1_olheat";"20260722-2_smp544_fly1_clbar"];
rows=cell(0,12);allTraces={};allRel={};
for i=1:numel(trials)
    D=load_trial(fullfile(base,trials(i)));E=find_declines(D,.75);
    for k=1:numel(E.index)
        rows(end+1,:)={trials(i),E.time(k),E.drop(k),E.preFwd(k),E.postFwd(k), ... %#ok<AGROW>
            E.preFwd(k)-E.postFwd(k),E.meanRot(k),E.maxRot(k),E.headingChange(k), ...
            E.postFwd(k)<=5,E.preFwd(k)-E.postFwd(k)>=1,E.headingChange(k)>=15};
    end
    make_trial_plot(E,trials(i),outputFolder);allTraces{end+1}=E.traces;allRel{end+1}=E.rel; %#ok<AGROW>
end
events=cell2table(rows,'VariableNames',{'trial','time','dffDropZ','preFwd','postFwd', ...
    'fwdDrop','meanAbsRot','maxAbsRot','headingChangeDeg','lowPostSpeed','slowdown','reorientation'});
events.trial=string(events.trial);writetable(events,fullfile(outputFolder,'smp544_slow_decrease_events.csv'));
make_summary(events,allTraces,allRel,outputFolder);
thresholds=[.5 .75 1];sens=cell(numel(thresholds)*numel(trials),4);q=0;
for a=1:numel(thresholds),for i=1:numel(trials),q=q+1;D=load_trial(fullfile(base,trials(i)));E=find_declines(D,thresholds(a));sens(q,:)={thresholds(a),trials(i),numel(E.index),mean(E.drop,'omitnan')};end,end
sensitivity=cell2table(sens,'VariableNames',{'thresholdZ','trial','nEvents','meanDropZ'});sensitivity.trial=string(sensitivity.trial);
writetable(sensitivity,fullfile(outputFolder,'smp544_slow_decrease_threshold_sensitivity.csv'));
disp(events);disp(groupcounts(events,'trial'));
end

function D=load_trial(folder)
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name),'dq','roi');
t=q.dq(1).t(:);fwd=q.dq(1).bvf(:);rot=rad2deg(q.dq(1).bvy(:));heading=q.dq(1).bh(:);d=q.roi.ts;if isvector(d),d=d(:)';end
n=min([numel(t),numel(fwd),numel(rot),numel(heading),size(d,2)]);z=zeros(n,1);
for r=1:size(d,1),z=z+norm_trace(d(r,1:n)',t(1:n));end;z=z/size(d,1);stim=false(n,1);
if isfield(q.dq(1),'heat_active')&&~isempty(q.dq(1).heat_active),stim=logical(q.dq(1).heat_active(1:n));
elseif isfield(q.dq(1),'heat')&&~isempty(q.dq(1).heat),h=q.dq(1).heat(1:n);stim=h(:)>median(h,'omitnan')+2*mad(h,1);end
D=struct('t',t(1:n),'dff',z,'fwd',fwd(1:n),'rot',rot(1:n),'heading',heading(1:n),'stim',stim);
end

function z=norm_trace(d,t)
fs=1/median(diff(t));sm=smoothdata(d(:),'gaussian',max(3,round(.25*fs)));res=sm-movmedian(sm,max(3,round(30*fs)),'omitmissing');s=1.4826*mad(res,1);if ~isfinite(s)||s<=0,s=std(res,'omitnan');end;z=res/s;
end

function E=find_declines(D,threshold)
fs=1/median(diff(D.t));pre=round(3*fs);post=round(3*fs);rel=(-pre:post)'/fs;
z=smoothdata(D.dff,'gaussian',max(3,round(.75*fs)));f=smoothdata(D.fwd,'gaussian',max(3,round(.2*fs)));r=smoothdata(abs(D.rot),'gaussian',max(3,round(.12*fs)));
forbidden=movmax(double(D.stim),[pre post])>0;forbidden(1:pre)=true;forbidden(end-post+1:end)=true;
[~,p]=findpeaks(-z,'MinPeakDistance',round(3*fs),'MinPeakProminence',threshold/2);
p=p(~forbidden(p));a=-round(2*fs):-round(fs);b=-round(.25*fs):round(.5*fs);h=-round(fs):round(fs);
keep=false(size(p));drop=nan(size(p));preF=drop;postF=drop;meanR=drop;maxR=drop;hd=drop;
for k=1:numel(p)
    drop(k)=mean(z(p(k)+a),'omitnan')-mean(z(p(k)+b),'omitnan');keep(k)=drop(k)>=threshold;
    preF(k)=mean(f(p(k)+a),'omitnan');postF(k)=mean(f(p(k)+b),'omitnan');meanR(k)=mean(r(p(k)+h),'omitnan');maxR(k)=max(r(p(k)+h));
    delta=D.heading(p(k)+h(end))-D.heading(p(k)+h(1));hd(k)=abs(rad2deg(atan2(sin(delta),cos(delta))));
end
p=p(keep);E=struct('index',p,'time',D.t(p),'drop',drop(keep),'preFwd',preF(keep),'postFwd',postF(keep),'meanRot',meanR(keep),'maxRot',maxR(keep),'headingChange',hd(keep),'rel',rel);
E.traces={extract(z,p,pre,post),extract(f,p,pre,post),extract(r,p,pre,post)};
end

function X=extract(x,e,pre,post)
X=nan(pre+post+1,numel(e));for i=1:numel(e),X(:,i)=x(e(i)-pre:e(i)+post);end
end

function make_trial_plot(E,name,out)
labs={'SMP544 activity (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};fig=figure('Visible','off','Color','w','Position',[50 50 1250 900]);tl=tiledlayout(3,1,'TileSpacing','compact');
for v=1:3,ax=nexttile(tl);hold(ax,'on');X=E.traces{v};if ~isempty(X),plot(ax,E.rel,X,'Color',[.82 .82 .82],'LineWidth',.6);plot(ax,E.rel,mean(X,2,'omitnan'),'Color',[.16 .38 .72],'LineWidth',3);end;xline(ax,0,':k');xlim(ax,[-3 3]);ylabel(ax,labs{v});grid(ax,'on');if v==3,xlabel(ax,'Time from fluorescence minimum (s)');end;end
title(tl,sprintf('%s | gradual SMP544 decreases (n=%d)',strrep(name,'_','\_'),numel(E.index)));savefig(fig,fullfile(out,name+'_slow_decreases.fig'));exportgraphics(fig,fullfile(out,name+'_slow_decreases.png'),'Resolution',220);close(fig);
end

function make_summary(T,traces,rels,out)
t=(-3:.05:3)';M=cell(1,3);for v=1:3,for i=1:numel(traces),X=traces{i}{v};if ~isempty(X),M{v}(:,end+1)=interp1(rels{i},mean(X,2,'omitnan'),t,'linear',NaN);end,end,end %#ok<AGROW>
fig=figure('Visible','off','Color','w','Position',[50 50 1450 850]);tl=tiledlayout(2,2,'TileSpacing','compact');labs={'SMP544 activity (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
for v=1:3,ax=nexttile(tl);hold(ax,'on');m=mean(M{v},2,'omitnan');se=std(M{v},0,2,'omitnan')/sqrt(size(M{v},2));fill(ax,[t;flipud(t)],[m-se;flipud(m+se)],[.16 .38 .72],'FaceAlpha',.18,'EdgeColor','none');plot(ax,t,m,'Color',[.16 .38 .72],'LineWidth',3);xline(ax,0,':k');xlim(ax,[-3 3]);ylabel(ax,labs{v});xlabel(ax,'Time from fluorescence minimum (s)');grid(ax,'on');end
ax=nexttile(tl);hold(ax,'on');scatter(ax,T.fwdDrop,T.headingChangeDeg,45,T.meanAbsRot,'filled');xline(ax,1,':k');yline(ax,15,':k');xlabel(ax,'Forward-speed decrease (mm/s)');ylabel(ax,'Heading change (deg)');cb=colorbar(ax);cb.Label.String='Mean rotational speed (deg/s)';grid(ax,'on');title(tl,sprintf('Gradual SMP544 decreases across good-behavior trials (n=%d)',height(T)));savefig(fig,fullfile(out,'smp544_slow_decreases_summary.fig'));exportgraphics(fig,fullfile(out,'smp544_slow_decreases_summary.png'),'Resolution',240);close(fig);
end
