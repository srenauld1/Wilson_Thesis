function export_pizza_0812_29_event_summaries(onlyTwoColumn,matchOneColumnToMinimal)
% SVG-only, one-recording summaries using original pizza-talk event rules.
% Thin lines: events, thick: arithmetic event mean, band: event SEM.
% Kinematics are the original saved imaging-aligned values; selection alone
% uses 0.25-s median behavior. No dF/F cap or amplitude normalization.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
if nargin<1,onlyTwoColumn=false;end
if nargin<2,matchOneColumnToMinimal=false;end
source='/Users/sophiarenauld/stacks/20250812-29_ves041_fly3_bbarjumps2 copy/20250812_29_1_or_z0z1_s_.mat';
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
a=load(source,'dq','roi');t=double(a.dq(1).t(:));
d=struct('t',t,'dff',double(a.roi.ts(1,:))','fwd',double(a.dq(1).bvf(:)), ...
 'yaw',rad2deg(double(a.dq(1).bvy(:))));
d.valid=isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);fs=1/median(diff(t));
base=movmedian(d.dff,max(3,round(30*fs)),'omitmissing');res=d.dff-base;
sigma=1.4826*mad(res(isfinite(res)),1);prominence=3*sigma;
[~,pk]=findpeaks(d.dff,'MinPeakProminence',prominence,'MinPeakDistance',max(1,round(2*fs)), ...
 'MinPeakWidth',max(1,round(.25*fs)));
pk=pk(arrayfun(@(i)complete(d,i,-2,3),pk));
f=smoothdata(d.fwd,'movmedian',max(3,round(.25*fs)));acc=gradient(f,t);
candidate=d.valid&acc<=-2;change=diff([false;candidate;false]);s=find(change==1);e=find(change==-1)-1;
slow=[];last=-Inf;
for j=1:numel(s)
 [~,q]=min(acc(s(j):e(j)));i=s(j)+q-1;t0=t(i);
 if t0-last<2||~complete(d,i,-2,3),continue,end
 pre=t>=t0-.75&t<t0-.25;post=t>=t0+.25&t<=t0+.75;
 before=median(f(pre));drop=before-median(f(post));if before<1||drop<1,continue,end
 slow(end+1,1)=i;last=t0; %#ok<AGROW>
end
angle=pizza_ves041_turn_angle(d,slow,[0 .5]);slow=slow(angle>=15);
rel=(-2:.05:3)';X={extract(d,pk,rel),extract(d,slow,rel)};
if ~matchOneColumnToMinimal
 draw(X,rel,{'VES041 peaks','Slowdowns + rotation'},false, ...
  fullfile(out,'20250812-29_summary_peaks_and_slowdowns_rotation.svg'));
end
if onlyTwoColumn
 fprintf('Restored two-column SVG only: %d peaks; %d slowdown+rotation events.\n',numel(pk),numel(slow));
 return
end
% A -3:-2 s per-event baseline requires the extra second to be valid.
baselineIdx=slow(arrayfun(@(i)complete(d,i,-3,3),slow));
rel=(-3:.05:3)';Y=extract(d,baselineIdx,rel);
b=mean(Y(rel>=-3&rel<=-2,:,1),1);Y(:,:,1)=Y(:,:,1)-b;
draw({Y},rel,{'Slowdowns + rotation'},true, ...
 fullfile(out,'20250812-29_summary_slowdowns_rotation_dff_baseline_subtracted.svg'),matchOneColumnToMinimal);
fprintf('8/12-29: %d peaks; %d slowdown+rotation events; %d baseline-eligible. Peak prominence %.4f dF/F.\n', ...
 numel(pk),numel(slow),numel(baselineIdx),prominence);
fprintf('Baseline: per-event mean -3:-2 s. Bands: SEM across events in this one recording, not across flies.\n');
end

function yes=complete(d,i,lo,hi)
t=d.t;left=find(t<=t(i)+lo,1,'last');right=find(t>=t(i)+hi,1,'first');
yes=~isempty(left)&&~isempty(right);
if yes,yes=all(d.valid(left:right))&&all(diff(t(left:right))<=2*median(diff(t)));end
end

function X=extract(d,idx,rel)
X=nan(numel(rel),numel(idx),3);vals={d.dff,d.fwd,abs(d.yaw)};
for j=1:numel(idx),for v=1:3,X(:,j,v)=interp1(d.t,vals{v},d.t(idx(j))+rel,'linear');end,end
assert(all(isfinite(X(:))));
end

function draw(X,rel,titles,baseline,filename,matchMinimal)
if nargin<6,matchMinimal=false;end
ncol=numel(X);fig=figure('Visible','off','Color','w','Position',[20 20 490*ncol 1050]);
tl=tiledlayout(fig,3,ncol,'TileSpacing','compact','Padding','loose');
colors={[.9 .02 .02],[0 0 0],[.05 .15 .95]};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
if baseline,labels{1}={'Baseline-subtracted','dF/F'};end
for v=1:3
 lim=[];
 for c=1:ncol
  ax=nexttile(tl,(v-1)*ncol+c);hold(ax,'on');V=X{c}(:,:,v);assert(size(V,2)>1);
  mu=mean(V,2);se=std(V,0,2)/sqrt(size(V,2));
  bandColor=.82+.18*colors{v};thinColor=.75+.25*colors{v};
  fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],bandColor,'EdgeColor','none');
  plot(ax,rel,V,'Color',thinColor,'LineWidth',.45);
  plot(ax,rel,mu,'Color',colors{v},'LineWidth',3);
  xline(ax,0,':','Color',[.45 .45 .45]);xlim(ax,rel([1 end]));
  ax.FontSize=18;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=ceil(rel(1)):floor(rel(end));
  if c==1,ylabel(ax,labels{v},'FontSize',21);end
  if v==1,title(ax,titles{c},'FontSize',22);end
  if v==3,xlabel(ax,'Time (s)','FontSize',21);else,ax.XTickLabel=[];end
  lim=[lim;V(:);mu-se;mu+se]; %#ok<AGROW>
 end
 bounds=[min(lim),max(lim)];pad=.05*diff(bounds);
 for c=1:ncol
  ax=nexttile(tl,(v-1)*ncol+c);
  if v==1
   if baseline&&~matchMinimal,ylim(ax,[-1 1]);else,ylim(ax,[0 4]);end
  elseif v==2&&(~baseline||matchMinimal),ylim(ax,[0 10]);
  elseif v==3
   if baseline&&~matchMinimal,ylim(ax,[0 100]);else,ylim(ax,[0 80]);end
  else,ylim(ax,bounds+[-pad pad]);end
 end
end
print(fig,filename,'-dsvg','-painters');close(fig);
end
