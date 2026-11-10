function plot_pizza_0812_29_short_median()
% Display-only 5-sample centered median applied identically to both native
% velocity traces. dF/F remains the original linear interpolation, unchanged.
% NaNs remain missing and the full filtered range is displayed (no clipping).
assert(~feature('ShowFigureWindows'),'Run headlessly.');
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
source=fullfile(out,'20250812-29_time_traces_176_226_native_frame_differences.mat');
S=load(source);t=S.nativeTime(:);windowSamples=5;
windowSeconds=windowSamples*median(diff(t));
forward=movmedian(S.forward(:),windowSamples,'includenan');
rotation=movmedian(S.rotation(:),windowSamples,'includenan');
dff=S.dff(:);nativeTime=t;
fig=figure('Visible','off','Color','w','Position',[40 40 1550 850]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
values={dff,forward,rotation};colors={'r','k','b'};
labels={'dF/F','Forward velocity (mm/s)','Rotational velocity (deg/s)'};
for k=1:3
 ax=nexttile(tl);plot(ax,t,values{k},colors{k},'LineWidth',1.4);
 grid(ax,'on');ax.FontSize=17;xlim(ax,[176 226.5]);ylabel(ax,labels{k});
 if k==3,yline(ax,0,'Color',[.3 .3 .3]);xlabel(ax,'Time (s)');end
end
title(tl,'20250812-29 | dF/F, forward and rotational velocity','FontSize',24);
subtitle(tl,sprintf('Behavior: identical 5-sample median (%.0f ms) | dF/F unchanged',windowSeconds*1000),'FontSize',16);
stem=fullfile(out,'20250812-29_time_traces_176_226_native_median_5samples');
print(fig,[stem '.svg'],'-dsvg','-painters');
exportgraphics(fig,[stem '.png'],'Resolution',150);
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '.mat'],'nativeTime','dff','forward','rotation','windowSamples','windowSeconds','source');
fprintf('Window %.1f ms; forward maximum %.2f -> %.2f mm/s; absolute rotation maximum %.2f -> %.2f deg/s\n', ...
 windowSeconds*1000,max(S.forward),max(forward),max(abs(S.rotation)),max(abs(rotation)));
end
