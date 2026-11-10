function plot_pizza_flatpath_palette_options(trial)
% Identical samples/limits/geometry: only the fluorescence colormap changes.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
if nargin<1,trial='20261020_12_1';end
S=load(fullfile(out,[trial '_flat_path_dff_native_p05_p95_scatter.mat']));
u=linspace(0,1,256)';
maps={parula(256), ...
 interp1([0 .35 .65 .85 1],[.04 .06 .22;.08 .30 .58;.12 .65 .66;1 .45 .02;1 1 .08],u), ...
 interp1([0 .3 .6 .82 1],[.035 .015 .12;.28 .05 .42;.7 .14 .45;.97 .48 .26;1 .98 .64],u), ...
 interp1([0 .55 .75 .9 1],[.20 .27 .36;.48 .57 .64;.76 .79 .8;.95 .22 .08;1 .02 .0],u)};
names={'1. Blue to yellow (current)','2. Blue / teal to bright yellow', ...
 '3. Purple / magenta to pale yellow','4. Cool gray to vivid red'};
fig=figure('Visible','off','Color','w','Position',[20 20 1550 1000]);
isStripe=strcmp(trial,'20261020_12_1');
if isStripe
 tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
else
 fig.Position=[20 20 1500 1150];
 tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
end
for k=1:4
 ax=nexttile(tl);hold(ax,'on');
 pad=max(.4,.04*max([range(S.x),range(S.y)]));
 xl=[min(S.x)-pad max(S.x)+pad];yl=[min(S.y)-pad max(S.y)+pad];
 if isStripe
  patch(ax,[xl fliplr(xl)],[yl(1) yl(1) -7.5 -7.5],[1 .928 .88],'EdgeColor','none');
 end
 markerSize=22;if ~isStripe,markerSize=7;end
 scatter(ax,S.x,S.y,markerSize,S.dn,'filled','MarkerEdgeColor','none');
 axis(ax,'equal');xlim(ax,xl);ylim(ax,yl);colormap(ax,maps{k});caxis(ax,S.limits);
 cb=colorbar(ax);cb.Label.String='dF/F';cb.FontSize=14;
 title(ax,names{k},'FontSize',19);ylabel(ax,'Y (mm)');ax.FontSize=14;ax.Box='off';
 if k==4||~isStripe,xlabel(ax,'X position (mm)');end
end
title(tl,[strrep(trial,'_',' ') ' | identical data and 5th–95th percentile limits'],'FontSize',21);
stem=fullfile(out,[trial '_flat_path_palette_options']);
exportgraphics(fig,[stem '.png'],'Resolution',150);print(fig,[stem '.svg'],'-dsvg','-painters');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '.mat'],'maps','names');
end
