function result=plot_pizza_bar_scatter_summary_mean_sd(root)
% Add 20-degree mean +/- sample SD to the summary only; reuse saved samples.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'bar_tuning','bar_velocity_scatter');
S=load(fullfile(out,'bar_scatter_summary.mat'),'R');R=S.R;
edges=-180:20:180;centers=edges(1:end-1)+10;
colors=[0 .447 .741;.85 .325 .098];signs=[1 -1];rows=cell(0,8);
fig=figure('Visible','off','Color','w','Position',[20 20 1650 1100]);cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,ceil(numel(R)/3),3,'TileSpacing','compact');
for k=1:numel(R)
 r=R(k);ax=nexttile(tl);hold(ax,'on');valid=r.valid;
 cat=r.category(valid);c=repmat([.65 .65 .65],nnz(valid),1);
 for j=1:2,c(cat==signs(j),:)=repmat(colors(j,:),nnz(cat==signs(j)),1);end
 h=scatter(ax,r.barPositionDeg(valid),r.dff(valid),8,c,'filled','MarkerFaceAlpha',.18,'MarkerEdgeAlpha',.18,'HandleVisibility','off');
 assert(numel(h.XData)==nnz(valid));
 bin=discretize(r.barPositionDeg,edges);
 for j=1:2
  mu=nan(18,1);sd=mu;
  for b=1:18
   y=r.dff(valid&r.category==signs(j)&bin==b);n=numel(y);
   if n>0,mu(b)=mean(y);end
   if n>=2,sd(b)=std(y,0);end
   rows(end+1,:)={r.trial,signs(j),edges(b),edges(b+1),centers(b),n,mu(b),sd(b)}; %#ok<AGROW>
  end
  errorbar(ax,centers,mu,sd,'-o','Color',colors(j,:),'LineWidth',1.5,'MarkerSize',3,'MarkerFaceColor',colors(j,:),'CapSize',3);
 end
 xlim(ax,[-180 180]);xticks(ax,-180:90:180);xlabel(ax,'Bar position (deg)');ylabel(ax,'VES041 \DeltaF/F');grid(ax,'on');
 title(ax,{strrep(r.trial,'_','\_'),char(r.reference)},'FontSize',10);
end
title(tl,'VES041 versus bar position: all points with mean +/- SD in 20-degree bins');
subtitle(tl,'Blue: bar velocity >+1 deg/s | Orange: <-1 deg/s | Gray: near-stationary points | Bars: sample SD, not SEM; directions calculated separately');
stem=fullfile(out,'ves041_all_bar_trials_all_points');
% Keep a points-only copy once; individual trial files are never modified.
for ext=[".png",".fig"]
 backup=stem+"_points_only"+ext;
 if ~isfile(backup),copyfile(stem+ext,backup);end
end
exportgraphics(fig,stem+".png",'Resolution',150);savefig_visible_on_open(fig,stem+".fig");
T=cell2table(rows,'VariableNames',{'trial','barVelocitySign','binStartDeg','binEndDeg','binCenterDeg','samples','meanDff','sdDff'});
assert(all(T.binEndDeg-T.binStartDeg==20));
for k=1:numel(R)
 for j=1:2
  sel=T.trial==R(k).trial&T.barVelocitySign==signs(j);
  assert(sum(T.samples(sel))==nnz(R(k).valid&R(k).category==signs(j)),'Binned sample count mismatch');
 end
end
writetable(T,stem+"_20deg_mean_sd.csv");
methods=struct('bins','[-180,-160), ... [160,180]; 20 degrees, plotted at bin centers', ...
 'mean','Arithmetic mean of valid calcium samples per trial, position bin, and bar-velocity direction', ...
 'sd','Sample standard deviation (N-1); not SEM or confidence interval; undefined for N<2', ...
 'minimumSamples','Mean shown for >=1 sample; SD shown for >=2; no occupancy threshold', ...
 'scope','Summary PNG/FIG only; individual trial figures unchanged; all valid scatter points retained', ...
 'source','Existing bar_scatter_summary.mat; same angle references, exclusions, velocity categories, and original dF/F');
save(stem+"_20deg_mean_sd.mat",'T','edges','centers','methods');
result=struct('summary',stem,'trials',numel(R));disp(result);
end
