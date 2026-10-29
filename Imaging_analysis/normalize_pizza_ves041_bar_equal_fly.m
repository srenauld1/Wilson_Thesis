function result=normalize_pizza_ves041_bar_equal_fly(out,includeDark)
% Scale each metric by full-recording within-fly mean and sample SD.
% Linear normalization of fly means is identical to normalizing each event first.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
if nargin<2,includeDark=false;end
rawStem=fullfile(out,'bar','ves041_bar_pooled_equal_fly_peaks_slowdowns');cohortLabel='bar';
if includeDark,rawStem=fullfile(out,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns');cohortLabel='dark, closed-loop bar, and bar-jump';end
S=load([rawStem '.mat']);
flies=S.flies;nf=numel(flies);reference=cell(nf,3);
for k=1:height(S.F)
 f=find(flies==S.F.flyId(k));A=load(fullfile(out,'analysis_data',S.F.trial(k)+"_events.mat"),'d');d=A.d;
 valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);vars={d.dff,d.fwd,abs(d.yaw)};
 for v=1:3,reference{f,v}=[reference{f,v};vars{v}(valid)];end
end
center=zeros(nf,3);scale=zeros(nf,3);rows=cell(0,5);metric=["dff","forwardVelocity","absoluteYawSpeed"];
for f=1:nf
 for v=1:3
  x=reference{f,v};center(f,v)=mean(x);scale(f,v)=std(x,0);assert(numel(x)>1&&isfinite(scale(f,v))&&scale(f,v)>0,'Degenerate normalization reference');
  z=(x-center(f,v))/scale(f,v);assert(abs(mean(z))<1e-10&&abs(std(z,0)-1)<1e-10);
  rows(end+1,:)={flies(f),metric(v),center(f,v),scale(f,v),numel(x)}; %#ok<AGROW>
 end
end
flyMeans=S.flyMeans;
for f=1:nf,for v=1:3,flyMeans(:,f,v,:)=(S.flyMeans(:,f,v,:)-center(f,v))/scale(f,v);end,end
counts=S.counts;contributingFlies=S.contributingFlies;rel=S.rel;
groupMean=nan(size(S.groupMean));
for c=1:4,groupMean(:,:,c)=squeeze(mean(flyMeans(:,:,:,c),2,'omitnan'));assert(nnz(counts(:,c)>0)==contributingFlies(c));end
assert(isequal(counts,S.counts)&&all(counts(:,2)==counts(:,3)+counts(:,4)));
% Check a separately normalized event mean for every fly/metric/group.
for f=1:nf
 allX=cell(1,4);trialIds=find(S.F.flyId==flies(f));
 for k=trialIds'
  A=load(fullfile(out,'analysis_data',S.F.trial(k)+"_events.mat"),'X');
  for c=1:4,allX{c}=cat(2,allX{c},A.X{c});end
 end
 for c=1:4,for v=1:3
  expected=mean((allX{c}(:,:,v)-center(f,v))/scale(f,v),2);actual=flyMeans(:,f,v,c);assert(isequal(isnan(expected),isnan(actual)));q=abs(expected-actual);assert(all(q(isfinite(q))<1e-10));
 end,end
end
normalization=cell2table(rows,'VariableNames',{'flyId','metric','fullRecordingMean','fullRecordingSD','validSamples'});
z=flyMeans(isfinite(flyMeans));bound=max(1,ceil(2*max(abs(z)))/2);sharedLimits=[-bound bound];
fig=figure('Visible','off','Color','w','Position',[10 10 1800 1300]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
labels={'VES041 peaks','All slowdowns','Slowdowns without turns','Slowdowns with turns'};
flyLabel=replace(flies,'_fly',' F');flyLabel=replace(flyLabel,'_unlabeled_inferred',' *');flyLabel=replace(flyLabel,'_inferred','*');
% Blue-white-red map: zero is white and common limits apply to every metric.
cmap=[linspace(.15,1,128)' linspace(.3,1,128)' ones(128,1);ones(128,1) linspace(1,.2,128)' linspace(1,.15,128)'];
for row=1:4
 for c=1:4
  ax=nexttile(tl);hold(ax,'on');
  if row==1
   im=imagesc(ax,rel,1:nf,flyMeans(:,:,1,c)');im.AlphaData=isfinite(im.CData);set(ax,'YDir','normal','Color',[.9 .9 .9]);ylim(ax,[.5 nf+.5]);
   yticks(ax,1:nf);yticklabels(ax,flyLabel);ax.YAxis.FontSize=8;caxis(ax,sharedLimits);colormap(ax,cmap);cb=colorbar(ax);cb.Label.String='Fly mean \DeltaF/F (z-score)';
   title(ax,{labels{c},sprintf('%d flies | %d events',contributingFlies(c),sum(counts(:,c)))});ylabel(ax,'Fly (* inferred identity)');
  else
   v=row-1;plot(ax,rel,flyMeans(:,:,v,c),'Color',[.7 .7 .7],'LineWidth',.9);plot(ax,rel,groupMean(:,v,c),'k','LineWidth',2.3);
   names={'\DeltaF/F (within-fly z-score)','Forward velocity (within-fly z-score)','|Yaw| speed (within-fly z-score)'};ylabel(ax,names{v});ylim(ax,sharedLimits);grid(ax,'on');yline(ax,0,':','Color',[.45 .45 .45]);
  end
  xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
 end
end
title(tl,['VES041 ' cohortLabel ' trials | normalized metrics | equal-fly pooled peaks and slowdowns'],'FontWeight','bold');
subtitle(tl,{'Each metric: (value - fly mean) / fly SD, using all valid samples from that fly''s full included recordings | All panels share the same z-score scale', ...
 sprintf('%d recordings | %d estimated flies (%d inferred) | Gray: fly means; black: equal-fly mean | Slowdown turn: >=15 deg over 0:+0.5 s',height(S.F),nf,nnz(contains(lower(flies),'inferred')))},'FontSize',10);
stem=[rawStem '_normalized'];exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
axs=findall(fig,'Type','axes');assert(numel(axs)==16);for ax=axs',if ax.Layout.Tile>4,assert(isequal(ax.YLim,sharedLimits));else,assert(isequal(ax.CLim,sharedLimits));end,end
methods=struct('reference',['All valid native imaging-aligned samples from full included ' cohortLabel ' recordings pooled within each fly; each sample counted once'], ...
 'transform','Each metric separately: (value - full-recording fly mean) / full-recording fly sample SD (N-1)', ...
 'absoluteYaw','Compute abs(yaw) first; use its own mean and SD. Negative z-scores mean below-usual rotational speed, not opposite turn direction.', ...
 'averaging','Unchanged source events: average within fly, then weight each contributing fly equally', ...
 'scale','All trace y-limits and heatmap color limits are identical; all finite fly-mean values retained', ...
 'interpretation','z=1 means one within-fly SD above usual activity or speed; equal standardized changes need not imply equal physical magnitudes');
writetable(normalization,[stem '_parameters.csv']);
save([stem '.mat'],'flies','flyMeans','groupMean','counts','contributingFlies','rel','normalization','center','scale','sharedLimits','methods','rawStem','-v7.3');
fid=fopen([stem '_methods.txt'],'w');assert(fid>=0);fields=fieldnames(methods);for j=1:numel(fields),fprintf(fid,'%s: %s\n',fields{j},methods.(fields{j}));end;fprintf(fid,'\nVerified full-recording references yield mean zero and sample SD one. Every normalized fly mean agrees with independently normalized saved events. All 16 data panels and shared scales verified. Original unnormalized pooled outputs retained.\n');fclose(fid);
result=struct('stem',stem,'sharedLimits',sharedLimits,'flies',nf,'events',sum(counts,1));disp(result);
end
