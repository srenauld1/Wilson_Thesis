function export_pizza_barrier_episodes_13_14()
assert(~feature('ShowFigureWindows'),'Use -nodisplay');
base='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
source=fullfile(base,'smr_analysis/heat_exploration/xy_heat_paths_03.fig');
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples/barrier_episodes_13_14_20261104';if ~isfolder(out),mkdir(out);end
old=openfig(source,'invisible');oa=findall(old,'Type','axes');
chosen=gobjects(1,2);for k=1:numel(oa),for j=1:2,if startsWith(string(oa(k).Title.String),sprintf('Barrier episode %d |',12+j)),chosen(j)=oa(k);end,end,end
assert(all(isgraphics(chosen)));lims=zeros(2,4);for j=1:2,lims(j,:)=[chosen(j).XLim chosen(j).YLim];end
spanX=max(lims(:,2)-lims(:,1));spanY=max(lims(:,4)-lims(:,3));
A=load(fullfile(base,'smr_analysis/fly_data.mat'),'fly_data');p=A.fly_data.protocol;
cmap=interp1([0 .5 1],[.05 .20 .8;.70 .70 .70;.9 .04 .03],linspace(0,1,256));
f=figure('Visible','off','Color','w','Position',[20 20 1700 510]);axesNew=gobjects(1,2);
for j=1:2
 ax=copyobj(chosen(j),f);axesNew(j)=ax;ax.Units='normalized';ax.Position=[.065+(j-1)*.48 .20 .385 .65];
 cx=mean(lims(j,1:2));cy=mean(lims(j,3:4));bx=cx+[-.5 .5]*spanX;by=cy+[-.5 .5]*spanY;
 % Redraw only the background for expanded equal-size viewports. Saved path
 % fluorescence, sample positions and start/end symbols are kept unchanged.
 delete(findall(ax,'Type','image'));delete(findall(ax,'Type','contour'));
 ix=find(p.trial_number==12+j&p.barrier_active>0,1);assert(~isempty(ix));
 [X,Y]=meshgrid(linspace(bx(1),bx(2),350),linspace(by(1),by(2),180));
 dist=(p.barrier_x_mm(ix)-X)*p.barrier_normal_x(ix)+(p.barrier_y_mm(ix)-Y)*p.barrier_normal_y(ix);
 Z=max(0,1-abs(dist)/p.barrier_radius_mm(ix));
 if j==2,rgb=repmat(1-.18*Z,1,1,3);else,rgb=cat(3,ones(size(Z)),1-.24*Z,1-.32*Z);end
 hold(ax,'on');im=image(ax,bx,by,rgb);uistack(im,'bottom');contour(ax,X,Y,dist,[0 0],'Color',[.55 .2 .15],'LineWidth',1.2);
 colormap(ax,cmap);axis(ax,'equal');xlim(ax,bx);ylim(ax,by);set(ax,'YDir','normal');ax.FontSize=14;ax.Box='off';
 ax.XLabel.String='Logged X (mm)';ax.YLabel.String='Logged Y (mm)';ax.Title.FontSize=18;
 sx=bx(1)+.06*spanX;sy=by(1)+.13*spanY;plot(ax,[sx sx+10],[sy sy],'k-','LineWidth',2.5);text(ax,sx+5,sy-.04*spanY,'10 mm','HorizontalAlignment','center','VerticalAlignment','top','FontSize',13);
 cb=colorbar(ax);cb.FontSize=12;cb.Label.String='ΔF/F';cb.Label.FontSize=16;
 % Reset after colorbar so both axes have the same physical dimensions.
 ax.Position=[.065+(j-1)*.48 .20 .385 .65];cb.Position=[.458+(j-1)*.48 .265 .012 .48];
 assert(max(abs(ax.CLim-chosen(j).CLim))<1e-12);
end
stem=fullfile(out,'20261020_11_2_barrier_episodes_13_heat_14_control_equal_xy_blue_red');
exportgraphics(f,[stem '.png'],'Resolution',180);exportgraphics(f,[stem '.pdf'],'ContentType','vector');print(f,[stem '.svg'],'-dsvg','-painters');savefig(f,[stem '.fig']);close(f);close(old);
methods=struct('sourceFigure',source,'episodes',[13 14],'sharedViewportMm',[spanX spanY],'xyAspect','One mm on X equals one mm on Y; equal physical panel dimensions and identical spans across panels','palette','Study blue-gray-red; low blue, high red','fluorescence','Original saved path CData and shared numeric color limits retained','background','Original programmed gradient, control virtual gradient and barrier centerlines reconstructed on expanded viewports','scaleBarMm',10);
fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);copyfile([mfilename('fullpath') '.m'],out);disp(methods);disp(stem);
end
