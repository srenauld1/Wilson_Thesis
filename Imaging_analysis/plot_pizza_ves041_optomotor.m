function result=plot_pizza_ves041_optomotor(root)
% Explore motion-onset and same-direction turn responses, retaining all bouts.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'gratings','optomotor_exploration');if ~isfolder(out),mkdir(out);end
M=readtable(fullfile(root,'manifest.csv'),'TextType','string');M=M(M.status=="complete"&M.category=="gratings",:);
rel=(-2:.05:4)';R=struct([]);allE=table();
settings=struct('motionThresholdRadS',.02,'minimumMotionSec',1.5,'turnPeakThresholdDegS',50, ...
 'turnSmoothingSec',.2,'windowSec',[-2 4],'baselineSec',[-1 -.25], ...
 'direction','Same signed fly yaw and visual-pattern velocity, as in existing optomotor analysis; no absolute CW/CCW label', ...
 'activity','Original saved ROI dF/F; baseline change shown in additional summary row', ...
 'turnSelection','One strongest direction-matched, 0.2 s median-smoothed yaw peak within each stimulus-motion bout', ...
 'averaging','Gray individual bouts; black within-trial mean; combined black mean gives equal weight to contributing trials (one fly per date)', ...
 'control','No compensatory turn means matched peak <50 deg/s; may include smaller matched turns or opposite-direction turns', ...
 'validity','Entire -2:+4 s event window finite and valid, with no gap >2 median sample intervals; invalid windows excluded from alignment, retained in bout table', ...
 'baseline','Each event subtracts its own mean dF/F at -1:-0.25 s relative to its alignment time; raw dF/F also shown');
for k=1:height(M)
 name=M.trial(k);fprintf('OPTOMOTOR %s\n',name);
 S=load(M.source(k),'fly_data');q=S.fly_data.dq(1);
 D=load(fullfile(root,'analysis_data',name+"_events.mat"),'d','groups');d=D.d;t=d.t;
 assert(isequal(t,double(q.t(:))));stim=rad2deg(double(q.vvy(:)));yaw=movmedian(d.yaw,max(1,round(.2/median(diff(t)))));
 valid=d.valid(:)&isfinite(stim)&isfinite(yaw)&isfinite(d.fwd)&isfinite(d.dff);
 starts=[];stops=[];directions=[];
 for direction=[1 -1]
  [s,e]=runs(valid&direction*stim>rad2deg(.02));
  keep=t(e)-t(s)+median(diff(t))>=1.5;starts=[starts;s(keep)];stops=[stops;e(keep)];directions=[directions;repmat(direction,nnz(keep),1)]; %#ok<AGROW>
 end
 [starts,ix]=sort(starts);stops=stops(ix);directions=directions(ix);nb=numel(starts);
 onset=nan(numel(rel),nb,5);turn=onset;rows=cell(0,14);
 for b=1:nb
  direction=directions(b);s=starts(b);e=stops(b);
  [peak,j]=max(direction*yaw(s:e));p=s+j-1;compensatory=peak>=50;
  vars=[direction*stim,direction*d.yaw,d.fwd,d.dff];
  onsetValid=complete(t,valid,t(s),rel);turnValid=compensatory&&complete(t,valid,t(p),rel);
  if onsetValid,onset(:,b,:)=align(t,vars,t(s),rel);end
  if turnValid,turn(:,b,:)=align(t,vars,t(p),rel);end
  peakLag=NaN;if ~isempty(D.groups{1}),lags=t(D.groups{1})-t(p);[~,ip]=min(abs(lags));peakLag=lags(ip);end
  pre=t>=t(s)-1&t<t(s)&valid;during=t>=t(s)&t<=t(e)&valid;
  delta=mean(d.dff(during))-mean(d.dff(pre));
  rows(end+1,:)={name,b,direction,t(s),t(e),t(p),peak,compensatory,onsetValid,turnValid,delta,peakLag,mean(d.fwd(during)),mean(direction*d.yaw(during))}; %#ok<AGROW>
 end
 E=cell2table(rows,'VariableNames',{'trial','bout','stimulusSign','onsetSec','offsetSec','matchedTurnPeakSec','matchedTurnPeakDegS','compensatory','onsetWindowValid','turnWindowValid','duringMinusPreDff','nearestVesPeakLagFromTurnSec','meanForwardDuring','meanMatchedYawDuring'});
 r=struct('trial',name,'source',M.source(k),'t',t,'stimDegS',stim,'d',d,'valid',valid,'starts',starts,'stops',stops,'directions',directions,'E',E,'onset',onset,'turn',turn);
 if isempty(R),R=r;allE=E;else,R(end+1)=r;allE=[allE;E];end %#ok<AGROW>
 folder=fullfile(out,name);if ~isfolder(folder),mkdir(folder);end
 writetable(E,fullfile(folder,'motion_bouts.csv'));save(fullfile(folder,'optomotor_events.mat'),'r','rel','settings');
 % Whole recording: stimulus, signed fly yaw, forward speed, and neural signal.
 fig=figure('Visible','off','Color','w','Position',[20 20 1550 950]);tl=tiledlayout(fig,4,1,'TileSpacing','compact');
 vals={stim,d.yaw,d.fwd,d.dff};labs={'Stimulus velocity (deg/s)','Fly yaw (deg/s)','Forward (mm/s)','VES041 \DeltaF/F'};
 for v=1:4
  ax=nexttile(tl);plot(ax,t,vals{v},'k');hold(ax,'on');shade(ax,t(starts),t(stops),directions);ylabel(ax,labs{v});xlim(ax,[t(1),t(end)]);grid(ax,'on');
 end
 xlabel(ax,'Recording time (s)');title(tl,strrep(name,'_','\_'));subtitle(tl,'Blue/orange shading: positive/negative stimulus motion; original neural and kinematic traces');
 saveplots(fig,fullfile(folder,'whole_recording'));
 % Every stimulus bout is shown, six columns per page, with no outcome selection.
 for page=1:ceil(nb/6)
  ids=(page-1)*6+1:min(page*6,nb);fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);tl=tiledlayout(fig,4,numel(ids),'TileSpacing','compact');
  for v=1:4
   for b=ids
    ax=nexttile(tl);m=t>=t(starts(b))-2&t<=t(stops(b))+2;
    plot(ax,t(m)-t(starts(b)),vals{v}(m),'k');hold(ax,'on');xline(ax,0,'b:');xline(ax,t(stops(b))-t(starts(b)),'b:');
    xline(ax,E.matchedTurnPeakSec(b)-E.onsetSec(b),'--','Color',[.85 .2 .2]);grid(ax,'on');
    if b==ids(1),ylabel(ax,labs{v});end
    if v==1,title(ax,sprintf('Bout %d | sign %+d\nmatched peak %.0f deg/s | turn %d',b,directions(b),E.matchedTurnPeakDegS(b),E.compensatory(b)));end
    if v==4,xlabel(ax,'Time from stimulus onset (s)');end
   end
  end
  title(tl,strrep(name,'_','\_'));subtitle(tl,'Every motion bout | Blue dotted: motion start/end; red dashed: strongest matched turn (including subthreshold)');
  saveplots(fig,fullfile(folder,sprintf('all_motion_bouts_page_%02d',page)));
 end
 summary_figure(r,rel,fullfile(folder,'motion_and_compensatory_turn_alignment'),false);
 fprintf('  %d motion bouts; %d compensatory; %d below cutoff; %d valid turn windows\n',nb,nnz(E.compensatory),nnz(~E.compensatory),nnz(E.turnWindowValid));
end
writetable(allE,fullfile(out,'all_motion_bouts.csv'));
summary_figure(R,rel,fullfile(out,'ves041_optomotor_equal_fly_summary'),true);
summary_figure(R,rel,fullfile(out,'ves041_optomotor_by_stimulus_direction'),true,true);
save(fullfile(out,'optomotor_summary.mat'),'R','allE','rel','settings');
fid=fopen(fullfile(out,'methods.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'trials',numel(R),'bouts',height(allE),'compensatoryBouts',nnz(allE.compensatory));disp(result);
end
function summary_figure(R,rel,stem,pooled,byDirection)
if nargin<5,byDirection=false;end
names={'All motion onsets','Motion onsets: compensatory turn','Motion onsets: below turn cutoff','Compensatory turn peaks'};
if byDirection,names={'Positive stimulus onsets','Negative stimulus onsets','Positive-stimulus turn peaks','Negative-stimulus turn peaks'};end
labs={'Direction-matched stimulus (deg/s)','Direction-matched fly yaw (deg/s)','Forward velocity (mm/s)','Original VES041 \DeltaF/F','Baseline-subtracted \DeltaF/F'};
fig=figure('Visible','off','Color','w','Position',[20 20 1650 1250]);tl=tiledlayout(fig,5,4,'TileSpacing','compact');
for v=1:5
 for c=1:4
  ax=nexttile(tl);hold(ax,'on');means=nan(numel(rel),numel(R));count=0;
  for k=1:numel(R)
   r=R(k);select=r.E.onsetWindowValid;turnAligned=c==4;
   if byDirection
    turnAligned=c>=3;if turnAligned,select=r.E.turnWindowValid;end
    direction=1;if mod(c,2)==0,direction=-1;end
    select=select&r.E.stimulusSign==direction;
   else
    if c==2,select=select&r.E.compensatory;elseif c==3,select=select&~r.E.compensatory;elseif c==4,select=r.E.turnWindowValid;end
   end
   if turnAligned,X=r.turn(:,select,v);else,X=r.onset(:,select,v);end
   count=count+nnz(select);
   if ~isempty(X)
    if ~pooled,plot(ax,rel,X,'Color',[.8 .8 .8],'LineWidth',.5);end
    means(:,k)=mean(X,2,'omitnan');
   end
  end
  if pooled,plot(ax,rel,means,'LineWidth',1);end
  plot(ax,rel,mean(means,2,'omitnan'),'k','LineWidth',2);xline(ax,0,'r:');grid(ax,'on');xlim(ax,[-2 4]);
  if c==1,ylabel(ax,labs{v});end
  if v==1,title(ax,sprintf('%s\n%d bouts | %d contributing flies',names{c},count,nnz(any(isfinite(means),1))));end
  if v==5,xlabel(ax,'Time from onset / turn peak (s)');yline(ax,0,':');end
 end
end
if pooled,title(tl,'VES041 optomotor exploration | 3 trials / 3 flies | equal-fly mean');subtitle(tl,'Thin colors: fly means; black: equal-fly mean | Positive matched yaw follows stimulus direction | Raw dF/F and event-baseline change shown');
else,title(tl,strrep(R.trial,'_','\_'));subtitle(tl,'Gray: every eligible bout; black: mean | Matched turn cutoff 50 deg/s | Motion threshold 0.02 rad/s for >=1.5 s');end
saveplots(fig,stem);
end
function X=align(t,vars,event,rel)
X=nan(numel(rel),1,5);
for v=1:4,X(:,1,v)=interp1(t,vars(:,v),event+rel);end
X(:,1,5)=X(:,1,4)-mean(X(rel>=-1&rel<=-.25,1,4));
end
function ok=complete(t,valid,event,rel)
m=t>=event+rel(1)&t<=event+rel(end);idx=find(m);
ok=event+rel(1)>=t(1)&&event+rel(end)<=t(end)&&numel(idx)>2&&all(valid(max(1,idx(1)-1):min(numel(t),idx(end)+1)))&&all(diff(t(idx))<=2*median(diff(t)));
end
function [s,e]=runs(mask),z=diff([false;mask(:);false]);s=find(z==1);e=find(z==-1)-1;end
function shade(ax,s,e,d)
limits=ylim(ax);cc=[0 .447 .741;.85 .325 .098];
for j=1:numel(s),color=cc(1+(d(j)<0),:);h=patch(ax,[s(j),e(j),e(j),s(j)],limits([1 1 2 2]),color,'FaceAlpha',.1,'EdgeColor','none');uistack(h,'bottom');end
ylim(ax,limits);
end
function saveplots(fig,stem)
exportgraphics(fig,string(stem)+".png",'Resolution',150);savefig_visible_on_open(fig,string(stem)+".fig");close(fig);
end
