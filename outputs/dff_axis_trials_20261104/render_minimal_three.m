function render_minimal_three
out='/Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/outputs/dff_axis_trials_20261104';
S=load(fullfile(out,'bar_dark_optomotor_aligned_data.mat'),'X','rel');X=S.X([1 2 4]);t=S.rel;
fig=figure('Visible','off','Color','w','Position',[20 20 1300 1000]);
cols={[.9 .02 .02],[0 0 0],[.05 .15 .95]};
left=[.15 .435 .72];w=.235;bottom=[.57 .305 .055];h=[.39 .225 .225];
for v=1:3
 for c=1:3
  ax=axes(fig,'Position',[left(c) bottom(v) w h(v)]);hold(ax,'on');V=X{c}(:,:,v);
  n=sum(isfinite(V),2);mu=mean(V,2,'omitnan');se=std(V,0,2,'omitnan')./sqrt(n);
  fill(ax,[t;flipud(t)],[mu-se;flipud(mu+se)],.82+.18*cols{v},'EdgeColor','none');
  plot(ax,t,mu,'Color',cols{v},'LineWidth',2.6);xlim(ax,[-2 3]);
  xline(ax,0,':','Color',[.25 .25 .25],'LineWidth',1.3);
  set(ax,'FontName','Helvetica Neue','FontSize',17,'Box','off','XTick',[],'YTick',[],'Color','none','Clipping','off');
  if v==1
   if c==1,yl=[0 3.5];else,yl=[.1 .1+3.5/1.5];end
   ylim(ax,yl);ax.XColor='none';ax.YColor='none';
   barbottom=yl(1)+.13*diff(yl);xx=-2.13;
   plot(ax,[xx xx],[barbottom barbottom+.5],'Color',[.15 .15 .15],'LineWidth',2.5);
   text(ax,xx-.10,barbottom+.25,'0.5','FontName','Helvetica Neue','FontSize',16,'HorizontalAlignment','right','VerticalAlignment','middle');
   if c==1,text(ax,-3.35,mean(yl),'\DeltaF/F','FontName','Helvetica Neue','FontSize',23,'HorizontalAlignment','center');end
  else
   if v==2,ylim(ax,[0 10]);ticks=0:2:10;lab={'Forward','velocity','(mm/s)'};else,ylim(ax,[0 80]);ticks=0:20:80;lab={'Rotational','speed (°/s)'};end
   ax.XColor='none';
   if c==1,ax.YTick=ticks;ylabel(ax,lab,'FontName','Helvetica Neue','FontSize',21,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');else,ax.YColor='none';end
   if v==3
    y=-3;plot(ax,[-2 -1],[y y],'Color',[.15 .15 .15],'LineWidth',2);
    text(ax,-1.5,y-2,'1 sec','FontName','Helvetica Neue','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','top');
   end
  end
 end
end
print(fig,fullfile(out,'dff_three_columns_minimal_double_scale.svg'),'-dsvg','-painters');
print(fig,fullfile(out,'dff_three_columns_minimal_double_scale.png'),'-dpng','-r120');close(fig);
end
