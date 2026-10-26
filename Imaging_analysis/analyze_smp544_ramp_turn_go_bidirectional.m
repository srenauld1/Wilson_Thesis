function results=analyze_smp544_ramp_turn_go_bidirectional(outputFolder)
% Independent physiology/behavior selection for SMP544 turn-and-go motifs.
% ROI 1 is right and ROI 2 is left, matching smp544_analysis.m.
if nargin<1||strlength(string(outputFolder))==0,outputFolder='/Users/sophiarenauld/stacks/smp544_behavior_quantification/ramp_turn_go_bidirectional';end
if ~isfolder(outputFolder),mkdir(outputFolder);end
base='/Users/sophiarenauld/stacks';trials=["20260720-3_smp544_fly2_dark";"20260720-4_smp544_fly2_dark2_volume";"20260721-2_smp544_fly1_clbar";"20260721-3_smp544_fly1_clbar";"20260721-4_smp544_fly1_olheat";"20260722-2_smp544_fly1_clbar"];
rng(544);R=cell(numel(trials),1);rows=cell(0,14);
for i=1:numel(trials)
 D=load_trial(fullfile(base,trials(i)));R{i}=analyze_trial(D);make_trial_figure(R{i},trials(i),outputFolder);make_diagnostic(D,R{i},trials(i),outputFolder);
 arms={R{i}.physRight,R{i}.physLeft,R{i}.beh};names=["right-physiology","left-physiology","behavior"];
 for a=1:3,E=arms{a};for k=1:numel(E.index),rows(end+1,:)={trials(i),names(a),E.time(k),E.rightDrop(k),E.leftDrop(k),E.preFwd(k),E.postFwd(k),E.postFwd(k)-E.preFwd(k),E.meanRot(k),E.maxRot(k),E.headingChange(k),E.rampDuration(k),E.hasTurnGo(k),E.hasPrecedingTrough(k)};end,end %#ok<AGROW>
end
events=cell2table(rows,'VariableNames',{'trial','selectionArm','time','rightDffDropZ','leftDffDropZ','preFwd','postFwd','fwdIncrease','meanAbsRot','maxAbsRot','headingChangeDeg','rampDurationS','hasTurnGo','hasPrecedingTrough'});events.trial=string(events.trial);events.selectionArm=string(events.selectionArm);writetable(events,fullfile(outputFolder,'smp544_ramp_turn_go_events.csv'));make_summary(R,outputFolder);results=struct('events',events,'trials',{R},'outputFolder',string(outputFolder));disp(groupsummary(events,{'selectionArm','trial'}));
end

function D=load_trial(folder)
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name),'dq','roi');t=q.dq(1).t(:);fwd=q.dq(1).bvf(:);rot=rad2deg(q.dq(1).bvy(:));heading=q.dq(1).bh(:);d=q.roi.ts;if isvector(d),d=d(:)';end
n=min([numel(t),numel(fwd),numel(rot),numel(heading),size(d,2)]);right=norm_trace(d(1,1:n)',t(1:n));if size(d,1)>=2,left=norm_trace(d(2,1:n)',t(1:n));else,left=nan(n,1);end;stim=false(n,1);
if isfield(q.dq(1),'heat_active')&&~isempty(q.dq(1).heat_active),stim=logical(q.dq(1).heat_active(1:n));elseif isfield(q.dq(1),'heat')&&~isempty(q.dq(1).heat),h=q.dq(1).heat(1:n);stim=h(:)>median(h,'omitnan')+2*mad(h,1);end
D=struct('t',t(1:n),'right',right,'left',left,'fwd',fwd(1:n),'rot',rot(1:n),'heading',heading(1:n),'stim',stim);
end
function z=norm_trace(d,t),fs=1/median(diff(t));sm=smoothdata(d(:),'gaussian',max(3,round(.25*fs)));res=sm-movmedian(sm,max(3,round(30*fs)),'omitmissing');s=1.4826*mad(res,1);if ~isfinite(s)||s<=0,s=std(res,'omitnan');end;z=res/s;end

function A=analyze_trial(D)
fs=1/median(diff(D.t));pre=round(5*fs);post=round(5*fs);rel=(-pre:post)'/fs;n=numel(D.t);zr=smoothdata(D.right,'gaussian',max(3,round(1.5*fs)));zl=smoothdata(D.left,'gaussian',max(3,round(1.5*fs)));f=smoothdata(D.fwd,'gaussian',max(3,round(.2*fs)));r=smoothdata(abs(D.rot),'gaussian',max(3,round(.12*fs)));
forbidden=movmax(double(D.stim),[round(3*fs) round(3*fs)])>0;forbidden(1:pre)=true;forbidden(end-post+1:end)=true;
[pr,rr,rejr]=find_ramps(zr,forbidden,fs);[pl,rl,rejl]=find_ramps(zl,forbidden,fs);pb=find_turn_go(f,r,D.heading,forbidden,fs);
A=struct('rel',rel);A.physRight=describe(pr,rr,D,zr,zl,f,r,fs,pre,post);A.physLeft=describe(pl,rl,D,zr,zl,f,r,fs,pre,post);A.beh=describe(pb,nan(size(pb)),D,zr,zl,f,r,fs,pre,post);A.rejectedRight=rejr;A.rejectedLeft=rejl;
allEvents=[pr;pl;pb];A.physRight=add_controls(A.physRight,forbidden,allEvents,pre,post,fs,zr,zl,f,r);A.physLeft=add_controls(A.physLeft,forbidden,allEvents,pre,post,fs,zr,zl,f,r);A.beh=add_controls(A.beh,forbidden,allEvents,pre,post,fs,zr,zl,f,r);
end

function [accepted,duration,rejected]=find_ramps(z,forbidden,fs)
[~,cand]=findpeaks(-z,'MinPeakDistance',round(5*fs),'MinPeakProminence',.4);cand=cand(~forbidden(cand));accepted=[];duration=[];rejected=[];
for p=cand(:)'
 preRange=p-round(5*fs):p-round(fs);postRange=p+round(.5*fs):p+round(3*fs);[preHigh,ii]=max(z(preRange));postHigh=max(z(postRange));drop=preHigh-z(p);recovery=postHigh-z(p);dur=(numel(preRange)-ii+1+round(fs))/fs;segment=z(preRange(ii):p);negativeFraction=mean(diff(segment)<=0);
 if drop>=1&&recovery>=.75&&dur>=1&&dur<=5&&negativeFraction>=.6,accepted(end+1,1)=p;duration(end+1,1)=dur;else,rejected(end+1,1)=p;end %#ok<AGROW>
end
end
function idx=find_turn_go(f,r,h,forbidden,fs)
[~,cand]=findpeaks(r,'MinPeakHeight',30,'MinPeakProminence',20,'MinPeakDistance',round(2*fs));cand=cand(~forbidden(cand));idx=[];
for p=cand(:)',preF=mean(f(p-round(2*fs):p-round(.75*fs)),'omitnan');postF=mean(f(p+round(.25*fs):p+round(1.5*fs)),'omitnan');hd=head_change(h,p,-round(.5*fs):round(.75*fs));if preF<=5&&postF-preF>=1&&hd>=15,idx(end+1,1)=p;end,end %#ok<AGROW>
end

function E=describe(idx,dur,D,zr,zl,f,r,fs,pre,post)
idx=idx(:);m=numel(idx);E=struct('index',idx,'time',D.t(idx),'rampDuration',dur(:),'traces',{{extract(zr,idx,pre,post),extract(zl,idx,pre,post),extract(f,idx,pre,post),extract(r,idx,pre,post)}});E.rightDrop=nan(m,1);E.leftDrop=nan(m,1);E.preFwd=nan(m,1);E.postFwd=nan(m,1);E.meanRot=nan(m,1);E.maxRot=nan(m,1);E.headingChange=nan(m,1);E.hasTurnGo=false(m,1);E.hasPrecedingTrough=false(m,1);
for k=1:m,p=idx(k);E.rightDrop(k)=mean(zr(p-round(3*fs):p-round(2*fs)))-min(zr(p-round(2*fs):p+round(.25*fs)));E.leftDrop(k)=mean(zl(p-round(3*fs):p-round(2*fs)))-min(zl(p-round(2*fs):p+round(.25*fs)));E.preFwd(k)=mean(f(p-round(2*fs):p-round(.75*fs)));E.postFwd(k)=mean(f(p+round(.25*fs):p+round(1.5*fs)));E.meanRot(k)=mean(r(p-round(.5*fs):p+round(.75*fs)));E.maxRot(k)=max(r(p-round(.5*fs):p+round(.75*fs)));E.headingChange(k)=head_change(D.heading,p,-round(.5*fs):round(.75*fs));E.hasTurnGo(k)=E.preFwd(k)<=5&&E.postFwd(k)-E.preFwd(k)>=1&&E.headingChange(k)>=15;E.hasPrecedingTrough(k)=max(E.rightDrop(k),E.leftDrop(k))>=1;end
end
function E=add_controls(E,forbidden,events,pre,post,fs,zr,zl,f,r),eligible=find(~forbidden);for p=events(:)',eligible(abs(eligible-p)<round(5*fs))=[];end;k=numel(E.index);if k&&~isempty(eligible),c=eligible(randperm(numel(eligible),min(k,numel(eligible))));else,c=[];end;E.control=c(:);E.controlTraces={extract(zr,c,pre,post),extract(zl,c,pre,post),extract(f,c,pre,post),extract(r,c,pre,post)};end
function hd=head_change(h,p,o),d=h(p+o(end))-h(p+o(1));hd=abs(rad2deg(atan2(sin(d),cos(d))));end
function X=extract(x,e,pre,post),X=nan(pre+post+1,numel(e));for i=1:numel(e),X(:,i)=x(e(i)-pre:e(i)+post);end,end

function make_trial_figure(A,name,out)
arms={A.physRight,A.physLeft,A.beh};titles={'Right dF/F troughs','Left dF/F troughs','Behavioral turn-and-go'};labs={'Right SMP544 (robust z)','Left SMP544 (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};cols={[.05 .32 .75],[.10 .65 .72],[.82 .24 .12]};fig=figure('Visible','off','Color','w','Position',[20 20 1800 1100]);tl=tiledlayout(4,3,'TileSpacing','compact');
for v=1:4,for a=1:3,ax=nexttile(tl);hold(ax,'on');X=arms{a}.traces{v};if ~isempty(X),plot(ax,A.rel,X,'Color',[.84 .84 .84],'LineWidth',.45,'HandleVisibility','off');plot(ax,A.rel,mean(X,2,'omitnan'),'Color',cols{a},'LineWidth',2.6,'DisplayName',sprintf('Events (n=%d)',size(X,2)));end;C=arms{a}.controlTraces{v};if ~isempty(C),plot(ax,A.rel,mean(C,2,'omitnan'),'--','Color',[.25 .25 .25],'LineWidth',1.5,'DisplayName','Random times');end;xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-5 5]);grid(ax,'on');if a==1,ylabel(ax,labs{v});end;if v==1,title(ax,titles{a});legend(ax,'Location','best');end;if v==4,xlabel(ax,'Time from selected trough or turn (s)');end;end,end
title(tl,sprintf('%s | independent selection, 1.5 s dF/F filter',strrep(name,'_','\_')));savefig(fig,fullfile(out,name+'_bidirectional.fig'));exportgraphics(fig,fullfile(out,name+'_bidirectional.png'),'Resolution',210);close(fig);
end

function make_diagnostic(D,A,name,out)
fs=1/median(diff(D.t));zr=smoothdata(D.right,'gaussian',max(3,round(1.5*fs)));zl=smoothdata(D.left,'gaussian',max(3,round(1.5*fs)));fig=figure('Visible','off','Color','w','Position',[40 40 1550 900]);tl=tiledlayout(4,1,'TileSpacing','compact');Y={zr,zl,D.fwd,abs(D.rot)};labs={'Right dF/F','Left dF/F','Forward speed','Rotational speed'};
for v=1:4,ax=nexttile(tl);plot(ax,D.t,Y{v},'k','LineWidth',.8);hold(ax,'on');for p=A.rejectedRight(:)',xline(ax,D.t(p),'Color',[.82 .82 .82],'HandleVisibility','off');end;for p=A.physRight.index(:)',xline(ax,D.t(p),'Color',[.05 .32 .75],'LineWidth',1.3,'HandleVisibility','off');end;for p=A.physLeft.index(:)',xline(ax,D.t(p),'Color',[.10 .65 .72],'LineWidth',1.3,'HandleVisibility','off');end;for p=A.beh.index(:)',xline(ax,D.t(p),'Color',[.82 .24 .12],'LineWidth',1.5,'HandleVisibility','off');end;ylabel(ax,labs{v});grid(ax,'on');if v==4,xlabel(ax,'Time (s)');end;end
title(tl,sprintf('%s | gray rejected, blue right trough, cyan left trough, red turn-and-go',strrep(name,'_','\_')));savefig(fig,fullfile(out,name+'_candidate_diagnostic.fig'));exportgraphics(fig,fullfile(out,name+'_candidate_diagnostic.png'),'Resolution',210);close(fig);
end

function make_summary(R,out)
t=(-5:.05:5)';arms={'physRight','physLeft','beh'};titles={'Right dF/F-selected troughs','Left dF/F-selected troughs','Behavior-selected turn-and-go'};labs={'Right SMP544 (robust z)','Left SMP544 (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};cols={[.05 .32 .75],[.10 .65 .72],[.82 .24 .12]};fig=figure('Visible','off','Color','w','Position',[20 20 1800 1100]);tl=tiledlayout(4,3,'TileSpacing','compact');
for v=1:4,for a=1:3,M=[];C=[];for i=1:numel(R),E=R{i}.(arms{a});if ~isempty(E.index),M(:,end+1)=interp1(R{i}.rel,mean(E.traces{v},2,'omitnan'),t,'linear',NaN);end;if ~isempty(E.control),C(:,end+1)=interp1(R{i}.rel,mean(E.controlTraces{v},2,'omitnan'),t,'linear',NaN);end,end %#ok<AGROW>
ax=nexttile(tl);hold(ax,'on');if ~isempty(M),m=mean(M,2,'omitnan');se=std(M,0,2,'omitnan')/sqrt(size(M,2));fill(ax,[t;flipud(t)],[m-se;flipud(m+se)],cols{a},'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');plot(ax,t,m,'Color',cols{a},'LineWidth',2.7,'DisplayName',sprintf('Events (%d trials)',size(M,2)));end;if ~isempty(C),plot(ax,t,mean(C,2,'omitnan'),'--','Color',[.2 .2 .2],'LineWidth',1.7,'DisplayName','Random times');end;xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-5 5]);grid(ax,'on');if a==1,ylabel(ax,labs{v});end;if v==1,title(ax,titles{a});legend(ax,'Location','best');end;if v==4,xlabel(ax,'Time from selected trough or turn (s)');end;end,end
title(tl,'SMP544 trough and turn-and-go analysis: independent selection');savefig(fig,fullfile(out,'smp544_ramp_turn_go_bidirectional_summary.fig'));exportgraphics(fig,fullfile(out,'smp544_ramp_turn_go_bidirectional_summary.png'),'Resolution',230);close(fig);
end
