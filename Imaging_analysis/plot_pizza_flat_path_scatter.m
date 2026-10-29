function info=plot_pizza_flat_path_scatter(d,name,stem,pathSource)
% Keep all valid points; compress color limits and draw brightest points last.
good=d.valid&isfinite(d.x)&isfinite(d.y)&isfinite(d.dff);ix=find(good);assert(~isempty(ix),'No valid path samples');
cl=prctile(d.dff(ix),[5 95]);if cl(1)==cl(2),cl=cl+[-.5 .5];end
[~,order]=sort(d.dff(ix),'ascend');ix=ix(order);markerSize=14;
fig=figure('Visible','off','Color','w','Position',[20 20 1100 950]);clean=onCleanup(@()close(fig));ax=axes(fig);
scatter(ax,d.x(ix),d.y(ix),markerSize,d.dff(ix),'filled','MarkerEdgeColor','none');axis(ax,'equal');grid(ax,'on');
colormap(ax,parula(256));caxis(ax,cl);cb=colorbar(ax);cb.Label.String='Original \DeltaF/F (5th-95th percentile colors)';
xlabel(ax,'X (mm)');ylabel(ax,'Y (mm)');title(ax,{strrep(name,'_',' '),'Flat path colored by original \DeltaF/F'});
notes={'Scatter: all valid samples; higher \DeltaF/F drawn last; colors saturated outside 5th-95th percentiles'};
if contains(pathSource,'.fig'),notes{end+1}='Coordinates recovered from original editable flat-path figure';end
subtitle(ax,notes,'FontSize',10);
exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
info=struct('low',cl(1),'high',cl(2),'markerSize',markerSize,'sampleCount',numel(ix),'drawOrder','ascending dF/F','percentiles',[5 95]);
end
