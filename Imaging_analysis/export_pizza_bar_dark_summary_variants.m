function export_pizza_bar_dark_summary_variants(includeGratings,zeroBaselineOnly,barOnly)
% SVG-only single-column, equal-fly summaries. Reuse saved event selections;
% extend their windows to -3:+3 only if all bracketing samples are valid.
% Pool events within each fly before equal weighting across flies. Baseline
% correction is per-event mean(-3:-2 s) subtraction, offset to baseline 1.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
if nargin<1,includeGratings=false;end
if nargin<2,zeroBaselineOnly=false;end
if nargin<3,barOnly=false;end
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
bar=load(fullfile(root,'bar','ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3.mat'));
dark=load(fullfile(root,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns.mat'));
rel=bar.rel(:);flies=dark.flies;trials=dark.F.trial;trialFlies=dark.F.flyId;
if includeGratings
 trials=[trials;"20250813-3_ves041_fly1_olgrating";"20250819-10_ves041_fly1_olgrating5"];
 extra=["20250813_fly1";"20250819_fly1"];
 assert(~any(ismember(extra,flies)));trialFlies=[trialFlies;extra];flies=[flies;extra];
end
pooled=cell(numel(flies),1);counts=zeros(numel(flies),1);
for k=1:numel(trials)
 trial=trials(k);f=find(flies==trialFlies(k));assert(isscalar(f));
 A=load(fullfile(root,'analysis_data',trial+"_events.mat"),'d');d=A.d;
 E=load(fullfile(root,'analysis_data',trial+"_events_with_rotation_peaks_and_speedups.mat"),'groups');
 ids=E.groups{4};valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
 vars={d.dff,d.fwd,abs(d.yaw)};dt=median(diff(d.t));
 for j=1:numel(ids)
  t0=d.t(ids(j));lo=find(d.t<=t0-3,1,'last');hi=find(d.t>=t0+3,1,'first');
  if isempty(lo)||isempty(hi)||any(~valid(lo:hi))||any(diff(d.t(lo:hi))>2*dt),continue,end
  X=nan(numel(rel),1,3);
  for v=1:3,X(:,1,v)=interp1(d.t,vars{v},t0+rel,'linear',NaN);end
  assert(all(isfinite(X(:))));pooled{f}=cat(2,pooled{f},X);counts(f)=counts(f)+1;
 end
end
combined=nan(numel(rel),numel(flies),3);
for f=1:numel(flies),if counts(f)>0,combined(:,f,:)=mean(pooled{f},2);end,end
% Preserve the existing bar panel's physical-unit kinematic axes.
src=fullfile(root,'bar','ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3_dff_event_baseline_slowdowns_rotation.fig');
ref=openfig(src,'invisible');axesRef=findall(ref,'Type','axes');yl=zeros(3,2);
for ax=axesRef',yl(ax.Layout.Tile,:)=ax.YLim;end
close(ref);yl(1,:)=[0 4];
data={bar.flyMeans(:,:,:,2),combined};cohorts={'bar','bar_dark'};
selected=1:2;window=[-3 3];
if includeGratings,selected=2;cohorts{2}='bar_dark_gratings';window=[-2 2];end
if barOnly,selected=1;window=[-2 2];end
for cohort=selected
 X=data{cohort};n=nnz(any(isfinite(X(:,:,1)),1));
 modes=[true false];if zeroBaselineOnly,modes=true;end
 for corrected=modes
  Y=X;
  if corrected
   baseline=mean(Y(rel>=-3&rel<=-2,:,1),1,'omitnan');
   baselineTarget=1;if zeroBaselineOnly,baselineTarget=0;end
   Y(:,:,1)=Y(:,:,1)-baseline+baselineTarget;
   assert(all(abs(mean(Y(rel>=-3&rel<=-2,:,1),1,'omitnan')-baselineTarget)<1e-10));
   suffix='baseline_dff';
  else,suffix='original_dff';end
  if includeGratings||barOnly,suffix=[suffix '_minus2_plus2'];end
  name=sprintf('ves041_%s_across_%d_flies_slowdowns_rotation_%s.svg',cohorts{cohort},n,suffix);
  limits=yl;if zeroBaselineOnly,limits(1,:)=[-1 1];end
  render(Y,rel,limits,corrected,fullfile(out,name),window,zeroBaselineOnly);
  fprintf('%s: %d flies\n',name,n);
 end
end
fprintf('Combined cohort: %d recordings, %d estimated flies, %d retained events.\n',numel(trials),numel(flies),sum(counts));
disp(table(flies,counts,'VariableNames',{'fly','events'}));
end

function render(X,t,limits,corrected,path,window,zeroBaseline)
fig=figure('Visible','off','Color','w','Position',[20 20 600 1000]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','loose');
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};
labels={'dF/F',{'Forward velocity','(mm/s)'},{'Rotational speed','(deg/s)'}};
if corrected
 if zeroBaseline,labels{1}={'Baseline-subtracted','dF/F'};else,labels{1}={'dF/F','(baseline = 1)'};end
end
for v=1:3
 ax=nexttile(tl);hold(ax,'on');V=X(:,:,v);mu=mean(V,2,'omitnan');n=sum(isfinite(V),2);
 se=std(V,0,2,'omitnan')./sqrt(n);
 fill(ax,[t;flipud(t)],[mu-se;flipud(mu+se)],.8+.2*colors{v},'EdgeColor','none');
 plot(ax,t,V,'Color',.45+.55*colors{v},'LineWidth',.65);
 plot(ax,t,mu,'Color',colors{v},'LineWidth',3);
 xline(ax,0,':','Color',[.45 .45 .45]);ylim(ax,limits(v,:));xlim(ax,window);
 ax.FontSize=19;ax.Box='off';grid(ax,'on');ax.GridAlpha=.1;ax.XTick=window(1):window(2);
 ylabel(ax,labels{v},'FontSize',22);
 if v==1,title(ax,'Slowdowns + rotation','FontSize',23);end
 if v<3,ax.XTickLabel=[];else,xlabel(ax,'Time (s)','FontSize',22);end
end
print(fig,path,'-dsvg','-painters');close(fig);
end
