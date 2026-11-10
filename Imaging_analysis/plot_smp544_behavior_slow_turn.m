function plot_smp544_behavior_slow_turn()
out=fullfile(fileparts(fileparts(mfilename('fullpath'))),'outputs','smp544_behavior_slow_turn_20261102');
if ~isfolder(out),mkdir(out);end
base='/Users/sophiarenauld/stacks';
trials=["20260720-3_smp544_fly2_dark";"20260720-4_smp544_fly2_dark2_volume";"20260721-2_smp544_fly1_clbar";"20260721-3_smp544_fly1_clbar";"20260721-4_smp544_fly1_olheat";"20260722-2_smp544_fly1_clbar"];
flyIds=[1 1 2 2 2 3];t=(-5:.05:5)';X=cell(3,4);rows=cell(0,4);rng(544);
for i=1:numel(trials)
 D=load_trial(fullfile(base,trials(i)));A=analyze_trial(D);E=A.beh;
 for v=1:4
  if ~isempty(E.index),X{flyIds(i),v}=[X{flyIds(i),v},interp1(A.rel,E.traces{v},t,'linear',NaN)];end
 end
 for k=1:numel(E.index),rows(end+1,:)={trials(i),flyIds(i),E.time(k),E.headingChange(k)};end
end
means=nan(numel(t),4,3);counts=zeros(3,1);
for j=1:3
 counts(j)=size(X{j,1},2);
 for v=1:4,if ~isempty(X{j,v}),means(:,v,j)=mean(X{j,v},2,'omitnan');end,end
end
colors=[.05 .32 .75;.05 .32 .75;0 0 0;.85 .12 .10];
ylabels={'Right SMP544 (robust z)','Left SMP544 (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
fig=figure('Visible','off','Color','w','Position',[100 100 730 1150]);tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
for v=1:4
 a=nexttile(tl);hold(a,'on');c=colors(v,:);Y=squeeze(means(:,v,:));m=mean(Y,2,'omitnan');n=sum(isfinite(Y),2);se=std(Y,0,2,'omitnan')./sqrt(n);se(n<2)=NaN;
 fill(a,[t;flipud(t)],[m-se;flipud(m+se)],c,'FaceAlpha',.18,'EdgeColor','none','DisplayName','SEM across flies');
 for j=1:3,plot(a,t,Y(:,j),'-','Color',.65+.35*c,'LineWidth',.8,'HandleVisibility','off');end
 plot(a,t,m,'-','Color',c,'LineWidth',2.8,'DisplayName',sprintf('Mean (%d flies; %d events)',sum(counts>0),sum(counts)));
 xlim(a,[-5 5]);ylabel(a,ylabels{v});set(a,'FontSize',14,'LineWidth',1,'Box','off');grid(a,'off');
 if v==1,title(a,'Behavior-selected slowdown + turn','FontSize',18);legend(a,'Location','northoutside','FontSize',12,'Box','off');end
 if v==4,xlabel(a,'Time from rotational-speed peak (s)');end
end
stem=fullfile(out,'smp544_behavior_slow_turn_mean_sem');savefig(fig,[stem '.fig']);exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');close(fig);
T=cell2table(rows,'VariableNames',{'Trial','Fly','TimeS','HeadingChangeDeg'});writetable(T,[stem '_events.csv']);save([stem '.mat'],'t','means','counts','trials','flyIds');disp(counts);disp(height(T));
end
function idx=find_slow_turn(f,r,h,forbidden,fs)
[~,cand]=findpeaks(r,'MinPeakHeight',30,'MinPeakProminence',20,'MinPeakDistance',round(2*fs));cand=cand(~forbidden(cand));idx=[];
for p=cand(:)'
 pre=mean(f(p-round(2*fs):p-round(fs)),'omitnan');near=mean(f(p-round(.5*fs):p+round(.25*fs)),'omitnan');hd=head_change(h,p,-round(.5*fs):round(.75*fs));
 if pre-near>=1 && hd>=15,idx(end+1,1)=p;end
end
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
[pr,rr,rejr]=find_ramps(zr,forbidden,fs);[pl,rl,rejl]=find_ramps(zl,forbidden,fs);pb=find_slow_turn(f,r,D.heading,forbidden,fs);
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

