function result=plot_pizza_optomotor_publication(root,useZscore,subtractBaseline,semBandsBoth)
% Two-column presentation view; source arrays and prior figures stay intact.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
if nargin<2,useZscore=true;end
if nargin<3,subtractBaseline=false;end
if nargin<4,semBandsBoth=false;end
out=fullfile(root,'gratings','optomotor_exploration');
S=load(fullfile(out,'ves041_optomotor_equal_fly_summary_zscored.mat'));
Q=S.summary;t=S.rel;columns=[1 4];variables=[2 3 4];
normalization=S.methods.normalization;activityLabel='VES041 activity (z-score)';suffix='zscored';
baselineTable=table();
if ~useZscore
 % Undo each fly's own transform before recomputing the equal-fly mean/SEM.
 for k=1:height(S.normalization)
  Q.flyMeans(:,k,4,:)=Q.flyMeans(:,k,4,:).*S.normalization.fullRecordingSdDff(k)+S.normalization.fullRecordingMeanDff(k);
 end
 n=sum(isfinite(Q.flyMeans),2);Q.mean=mean(Q.flyMeans,2,'omitnan');
 Q.sem=std(Q.flyMeans,0,2,'omitnan')./sqrt(n);Q.sem(n<2)=NaN;
 raw=load(fullfile(out,'optomotor_summary.mat'),'R');
 for k=1:numel(raw.R)
  expected=mean(raw.R(k).onset(:,raw.R(k).E.onsetWindowValid,4),2);
  assert(max(abs(expected-Q.flyMeans(:,k,4,1)))<1e-10,'Raw dF/F reconstruction mismatch');
  expected=mean(raw.R(k).turn(:,raw.R(k).E.turnWindowValid,4),2);
  assert(max(abs(expected-Q.flyMeans(:,k,4,4)))<1e-10,'Raw turn-aligned dF/F mismatch');
 end
 normalization='Original source ROI dF/F; no z scoring or per-event baseline subtraction';
 activityLabel='VES041 activity (\DeltaF/F)';suffix='raw_dff';
end
if subtractBaseline
 assert(~useZscore,'Baseline-subtracted display uses original dF/F units');
 for k=1:numel(raw.R)
  r=raw.R(k);nb=height(r.E);baseline=nan(1,nb);
  for b=1:nb
   before=r.valid&r.t>=r.E.onsetSec(b)-1&r.t<=r.E.onsetSec(b)-.25;
   assert(any(before),'No pre-stimulus baseline samples');
   baseline(b)=mean(r.d.dff(before));
  end
  for c=1:4
   select=r.E.onsetWindowValid;
   if c==2,select=select&r.E.compensatory;elseif c==3,select=select&~r.E.compensatory;elseif c==4,select=r.E.turnWindowValid;end
   if c==4,X=r.turn(:,select,4);else,X=r.onset(:,select,4);end
   corrected=X-baseline(select);
   Q.flyMeans(:,k,4,c)=mean(corrected,2,'omitnan');
   % Identical bout-specific PRE-STIMULUS offsets in both alignments.
   expected=mean(X,2,'omitnan')-mean(baseline(select));
   assert(max(abs(expected-Q.flyMeans(:,k,4,c)))<1e-10);
  end
  baselineTable=[baselineTable;table(repmat(string(r.trial),nb,1),r.E.bout,r.E.onsetSec,baseline(:), ...
   'VariableNames',{'trial','bout','stimulusOnsetSec','preStimulusMeanDff'})]; %#ok<AGROW>
 end
 n=sum(isfinite(Q.flyMeans),2);Q.mean=mean(Q.flyMeans,2,'omitnan');
 Q.sem=std(Q.flyMeans,0,2,'omitnan')./sqrt(n);Q.sem(n<2)=NaN;
 normalization='Original dF/F minus each bout mean over native valid samples -1:-0.25 s before stimulus onset; same baseline used for turn alignment; no z scoring';
 activityLabel='VES041 activity\newline(baseline-subtracted \DeltaF/F)';suffix='baseline_subtracted';
end
colors=[.15 .43 .70;.86 .34 .17;.48 .31 .65];motionColor=[.75 .87 .88];
fig=figure('Visible','off','Color','w','Position',[20 20 1340 1320]);
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');
labels={'Fly yaw (deg/s)','Forward speed (mm/s)',activityLabel};
axs=gobjects(3,2);yLimits=zeros(3,2,2);legendHandles=[];
for row=1:3
 v=variables(row);
 for col=1:2
  c=columns(col);ax=nexttile(tl);axs(row,col)=ax;hold(ax,'on');
  m=Q.mean(:,1,v,c);e=Q.sem(:,1,v,c);flies=Q.flyMeans(:,:,v,c);
  z=[flies(:);m-e;m+e];z=z(isfinite(z));limits=[min(z),max(z)];
  if diff(limits)==0,limits=limits+[-.5 .5];end
  % Individual panel limits expand the forward/neural traces relative to
  % the old shared limits dominated by the weak/no-turn comparison.
  padding=.10;if col==2&&row>=2,padding=.16;end
  limits=limits+[-1 1]*padding*diff(limits);
  if subtractBaseline&&row==3,limits=[-1 1];end
  yLimits(row,col,:)=reshape(limits,1,1,2);
  if col==1
   motion=patch(ax,[0 2 2 0],limits([1 1 2 2]),motionColor,'EdgeColor','none','FaceAlpha',.48);
  end
  if col==1||semBandsBoth
   ok=isfinite(m)&isfinite(e);
   fill(ax,[t(ok);flipud(t(ok))],[m(ok)-e(ok);flipud(m(ok)+e(ok))],[.2 .2 .2],'FaceAlpha',.14,'EdgeColor','none','HandleVisibility','off');
  end
  h=gobjects(1,4);
  for k=1:3,h(k)=plot(ax,t,flies(:,k),'Color',colors(k,:),'LineWidth',1.7);end
  h(4)=plot(ax,t,m,'Color',[.04 .04 .04],'LineWidth',3.2);
  if col==2&&~semBandsBoth
   % Sparse error bars avoid any shaded objects in the turn-aligned column.
   ix=1:10:numel(t);errorbar(ax,t(ix),m(ix),e(ix),'LineStyle','none','Color',[.12 .12 .12],'LineWidth',1.15,'CapSize',5,'HandleVisibility','off');
  end
  if row==1&&col==1,legendHandles=[h motion];end
  xline(ax,0,':','Color',[.4 .4 .4],'LineWidth',1.2,'HandleVisibility','off');
  if row==3,yline(ax,0,':','Color',[.65 .65 .65],'LineWidth',1,'HandleVisibility','off');end
  xlim(ax,[-2 4]);ylim(ax,limits);xticks(ax,-2:1:4);
  if subtractBaseline&&row==3,yticks(ax,-1:.5:1);end
  set(ax,'FontName','Arial','FontSize',20,'LineWidth',1.1,'Box','off','TickDir','out','TickLength',[.012 .012]);
  ylabel(ax,labels{row},'FontSize',23);
  if row<3,ax.XTickLabel=[];
  elseif col==1,xlabel(ax,'Time from visual motion onset (s)','FontSize',22);
  else,xlabel(ax,'Time from matched turn peak (s)','FontSize',22);end
  if row==1
   headings={'All visual motion','Matched optomotor turns'};
   title(ax,{headings{col},sprintf('%d bouts | 3 flies',sum(Q.boutCounts(:,c)))},'FontSize',25,'FontWeight','normal');
  end
 end
end
lg=legend(legendHandles,{'Aug 12, fly 3','Aug 13, fly 1','Aug 19, fly 1','Equal-fly mean +/- SEM','Visual motion on'}, ...
 'Orientation','horizontal','NumColumns',3,'FontSize',18,'Box','off');lg.Layout.Tile='south';
title(tl,'VES041 during optomotor stimulation','FontSize',30,'FontWeight','bold');
if semBandsBoth
 subtitle(tl,{'Mean +/- SEM across flies', 'Baseline: -1 to -0.25 s before visual motion onset, for both columns'},'FontSize',18);
else
 subtitle(tl,'Mean +/- SEM across flies: band at left, error bars at right','FontSize',18);
end
assert(numel(findall(fig,'Type','axes'))==6);
for row=1:3
 if semBandsBoth
  assert(numel(findall(axs(row,2),'Type','patch'))==1,'Right column should contain only the SEM band');
  assert(isempty(findall(axs(row,2),'Type','errorbar')));
 else
  assert(isempty(findall(axs(row,2),'Type','patch')),'Right column must have no shading');
  assert(numel(findall(axs(row,2),'Type','errorbar'))==1);
 end
end
stem=fullfile(out,['ves041_optomotor_publication_' suffix]);
exportgraphics(fig,[stem '.png'],'Resolution',300);savefig_visible_on_open(fig,[stem '.fig']);close(fig);
methods=struct('source',fullfile(out,'ves041_optomotor_equal_fly_summary_zscored.mat'), ...
 'columns','Original columns 1 and 4 only; all-motion onset and matched-turn peak alignment', ...
 'rows','Original rows 2,3,4 only: yaw, forward speed, neural dF/F (see normalization)', ...
 'shading','Visual-motion interval 0-2 s at left only; schematic interval, actual bout boundaries vary; left SEM band; no shading at right', ...
 'error','Between-fly SEM reused unchanged; right error bars every 0.5 s, with unchanged continuous mean', ...
 'axes','Panel-specific limits include all fly means and mean +/- SEM; 10 percent padding, 16 percent in right forward/neural panels', ...
 'normalization',normalization,'mean',S.methods.mean);
if semBandsBoth
 methods.error='Shaded +/-SEM across fly means in both columns; no error bars';
 methods.shading='Visual-motion interval 0-2 s at left only; SEM bands in both columns';
end
save([stem '.mat'],'methods','yLimits','columns','variables','Q','t','baselineTable');
if subtractBaseline,writetable(baselineTable,[stem '_baselines.csv']);end
result=struct('stem',stem,'axes',6);disp(result);
end
