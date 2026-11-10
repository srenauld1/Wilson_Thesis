function plot_ves041_boundary_path_314_325()
% Exact 314-325-s path; original dF/F, interval-specific 5th/95th color limits.
% Saturate colors only: do not discard samples or change fluorescence data.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics', ...
    'ves041_20261020_12_1_boundary_turns_20261102');
S=load(fullfile(out,'boundary_turns_data.mat'),'t','x','y','dff','inside');
keep=S.t>=314&S.t<=325&isfinite(S.x)&isfinite(S.y)&isfinite(S.dff);
ix=find(keep);x=S.x(ix);y=S.y(ix);d=S.dff(ix);
limits=prctile(d,[5 95]);assert(diff(limits)>0);
fig=figure('Visible','off','Color','w','Position',[20 20 1450 480]);
ax=axes(fig);hold(ax,'on');
patch(ax,[min(x)-.5 max(x)+.5 max(x)+.5 min(x)-.5], ...
    [-9.5 -9.5 -7.5 -7.5],[1 .85 .8],'EdgeColor','none');
xx=x;yy=y;gap=[false;diff(ix)>1];xx(gap)=NaN;yy(gap)=NaN;
surface(ax,[xx xx],[yy yy],zeros(numel(ix),2),[d d], ...
    'FaceColor','none','EdgeColor','interp','LineWidth',2.8);
plot(ax,x(1),y(1),'go','MarkerFaceColor','g','MarkerSize',7);
plot(ax,x(end),y(end),'ks','MarkerFaceColor','k','MarkerSize',6);
axis(ax,'equal');xlim(ax,[min(x)-.5 max(x)+.5]);ylim(ax,[min(y)-.4 max(y)+.4]);
colormap(ax,parula(256));caxis(ax,limits);cb=colorbar(ax);cb.Label.String='dF/F (5th-95th percentile color limits)';
xlabel(ax,'Logged X (mm)');ylabel(ax,'Logged Y (mm)');grid(ax,'on');
title(ax,'20261020 12_1 | 314-325 s | original dF/F');
subtitle(ax,sprintf('Color limits %.3f to %.3f, calculated within 314-325 s | green circle: start; black square: end',limits));
stem=fullfile(out,'flat_path_dff_314_325_p05_p95');
exportgraphics(fig,[stem '.png'],'Resolution',180);
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '.mat'],'ix','limits');disp(limits);
end
