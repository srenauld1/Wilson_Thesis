function plot_pizza_0812_29_imaging_scatter()
% Acquired imaging samples only, using their matching saved dq(1) path.
% No upsampling or smoothing. Color saturation does not remove any samples.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
source='/Users/sophiarenauld/stacks/20250812-29_ves041_fly3_bbarjumps2 copy/20250812_29_1_or_z0z1_s_.mat';
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
a=load(source,'dq','roi');
t=double(a.dq(1).t(:));x=double(a.dq(1).pxb(:));y=double(a.dq(1).pyb(:));dff=double(a.roi.ts(1,:))';
assert(isequal(size(t),size(x),size(y),size(dff)));
valid=isfinite(t)&isfinite(x)&isfinite(y)&isfinite(dff);
limits=prctile(dff(isfinite(dff)),[5 95]);markerSize=15;
fig=figure('Visible','off','Color','w','Position',[40 40 1200 850]);
ax=axes(fig);hold(ax,'on');scatter(ax,x(valid),y(valid),markerSize,dff(valid),'filled','MarkerEdgeColor','none');
axis(ax,'equal');pad=.04*max([range(x(valid)),range(y(valid))]);
xlim(ax,[min(x(valid))-pad,max(x(valid))+pad]);ylim(ax,[min(y(valid))-pad,max(y(valid))+pad]);
colormap(ax,parula(256));caxis(ax,limits);cb=colorbar(ax);cb.FontSize=18;cb.Label.String='dF/F';
ax.FontSize=18;ax.Box='off';xlabel(ax,'X position (mm)','FontSize',22);ylabel(ax,'Y position (mm)','FontSize',22);
title(ax,'20250812-29 | imaging-rate samples','FontSize',24);
stem=fullfile(out,'20250812-29_flat_path_dff_imaging_rate_scatter15_p05_p95');
exportgraphics(fig,[stem '.png'],'Resolution',150);print(fig,[stem '.svg'],'-dsvg','-painters');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '.mat'],'t','x','y','dff','valid','limits','markerSize','source');
fprintf('Imaging %.3f Hz; %d/%d samples plotted; marker area %g points squared; limits %.4f %.4f\n', ...
 1/median(diff(t)),nnz(valid),numel(t),markerSize,limits);
end
