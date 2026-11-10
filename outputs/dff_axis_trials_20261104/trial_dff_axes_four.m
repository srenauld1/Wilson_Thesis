function trial_dff_axes_four
addpath('/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Imaging_analysis');
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


B=load('/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/bar/ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3.mat');
keep=B.rel>=-2 & B.rel<=3;
assert(max(abs(B.rel(keep)-rel))<1e-8);
X{3}=B.flyMeans(keep,:,:,1);X{4}=B.flyMeans(keep,:,:,2);
nPeak=sum(any(isfinite(X{3}(:,:,1)),1));nSlow=sum(any(isfinite(X{4}(:,:,1)),1));
fprintf('Cohort peaks: %d flies, %d events; slowdowns: %d flies, %d events.\n',nPeak,sum(B.counts(:,1)),nSlow,sum(B.counts(:,2)));
out='/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/outputs/dff_axis_trials_20261104';
for zoom=0
 fig=figure('Visible','off','Color','w','Position',[20 20 2000 1300]);
 tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
 colors={[.9 .02 .02],[0 0 0],[.05 .15 .95]};
 for v=1:3
  for c=1:4
   if v==1,ax=nexttile(tl,c,[2 1]);else,ax=nexttile(tl,(v)*4+c);end
   hold(ax,'on');V=X{c}(:,:,v);mu=mean(V,2,'omitnan');se=std(V,0,2,'omitnan')./sqrt(sum(isfinite(V),2));
   fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],.82+.18*colors{v},'EdgeColor','none');
   % Cohort columns show equal-fly means and SEM only.
   plot(ax,rel,mu,'Color',colors{v},'LineWidth',2.5);
   xline(ax,0,':','Color',[.5 .5 .5]);xlim(ax,[-2 3]);
   set(ax,'FontName','Helvetica Neue','FontSize',17,'Box','off');
   if v==1
    if zoom
     lims={[0 3.6],[.6 2],[0 2.8]};ylim(ax,lims{c});
    else,ylim(ax,[0 3.5]);yticks(ax,0:.5:3.5);end
    titles={'Single recording: fluorescence peaks','Single recording: slowdowns + rotation',sprintf('Across %d flies: fluorescence peaks',nPeak),sprintf('Across %d flies: slowdowns + rotation',nSlow)};title(ax,titles{c},'FontSize',15);
    if c==1,ylabel(ax,'\DeltaF/F');end
   elseif v==2,ylim(ax,[0 10]);if c==1,ylabel(ax,{'Forward velocity','(mm/s)'});end
   else,ylim(ax,[0 80]);xlabel(ax,'Time from event (s)');if c==1,ylabel(ax,{'Rotational speed','(°/s)'});end
   end
   if v<3,ax.XTickLabel=[];end
  end
 end
 if zoom,name='dff_taller_individual_scales';else,name='dff_four_columns_shared_0_3p5';end
 print(fig,fullfile(out,[name '.svg']),'-dsvg','-painters');
 print(fig,fullfile(out,[name '.png']),'-dpng','-r100');close(fig);
end
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

