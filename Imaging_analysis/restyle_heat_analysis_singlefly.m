function restyle_heat_analysis_singlefly()
assert(~feature('ShowFigureWindows'),'Use -nodisplay');
base='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/pooled_four_original_pizza_rules_20261102';
S=load(fullfile(base,'crunch_data.mat'));k=find(S.ids=="20261020_11_2");assert(isscalar(k)&&isequal(S.counts(k,:),[20 19]));
out='/Users/sophiarenauld/Downloads/heat_analysis_singlefly_reference_style';
f=figure('Visible','off','Color','w','Position',[20 20 1100 850]);
colors={[1 .16 .12],[.08 .08 .08],[.05 .18 .9]};bands={[1 .79 .79],[.80 .80 .80],[.81 .84 1]};
ylimits={[0 2],[-1 10],[0 150]};ticks={0:1:2,0:5:10,0:50:150};mu=cell(1,2);se=mu;
for c=1:2,mu{c}=squeeze(mean(S.X{k,c},2));se{c}=squeeze(std(S.X{k,c},0,2))/sqrt(S.counts(k,c));end
for r=1:3
 for c=1:2
  x=.27+(c-1)*.40;y=.72-(r-1)*.295;
  ax=axes(f,'Position',[x y .305 .22]);hold(ax,'on');t=S.rel(:);m=mu{c}(:,r);e=se{c}(:,r);
  fill(ax,[t;flipud(t)],[m-e;flipud(m+e)],bands{r},'EdgeColor','none','FaceAlpha',1);
  plot(ax,t,m,'Color',colors{r},'LineWidth',2.1);xline(ax,0,':','Color',[.55 .55 .55],'LineWidth',1);
  xlim(ax,[-2 3]);ylim(ax,ylimits{r});ax.YTick=ticks{r};ax.XTick=[];ax.Box='off';ax.FontSize=18;ax.LineWidth=1;ax.XColor=[.12 .12 .12];ax.YColor=[.12 .12 .12];
  ax.XAxis.Visible='off';if c==2,ax.YTickLabel=[];end
  if r==3
   if c==1,ylim(ax,[0 80]);ax.YTick=0:20:80;else,ylim(ax,[0 150]);ax.YTick=0:50:150;ax.YTickLabel={'0','50','100','150'};end
  end
  if r==3
   line(ax,[-2 -1],[-.23 -.23]*ax.YLim(2),'Color',[.12 .12 .12],'LineWidth',2.5,'Clipping','off');text(ax,-2,-.30*ax.YLim(2),'1 s','FontSize',16,'Color',[.12 .12 .12],'VerticalAlignment','top','Clipping','off');
  end
 end
end
annotation(f,'textbox',[.02 .79 .20 .06],'String','ΔF/F','FontSize',25,'FontWeight','bold','HorizontalAlignment','center','EdgeColor','none');
annotation(f,'textbox',[.02 .445 .20 .16],'String',{'Forward','Velocity','(mm/s)'},'FontSize',22,'FontWeight','bold','HorizontalAlignment','center','EdgeColor','none');
annotation(f,'textbox',[.02 .15 .20 .16],'String',{'Rotational','Speed','(°/s)'},'FontSize',22,'FontWeight','bold','HorizontalAlignment','center','EdgeColor','none');
exportgraphics(f,[out '.png'],'Resolution',180);exportgraphics(f,[out '.pdf'],'ContentType','vector','BackgroundColor','none');print(f,[out '.svg'],'-dsvg','-painters');savefig(f,[out '.fig']);close(f);
methods=struct('inputSVG','/Users/sophiarenauld/Downloads/heat_analysis_singlefly.svg','source',fullfile(base,'crunch_data.mat'),'recording','20261020_11_2','counts',[20 19],'columns',{{'Heat-off','Heat-on'}},'display','Mean and event SEM only, original arrays unchanged; individual events and headings hidden to match reference','axisLimits',{{[0 2],[-1 10],[0 150]}},'timeWindow',[-2 3],'scaleBarSeconds',1);
save([out '_data.mat'],'mu','se','methods');fid=fopen([out '_methods.json'],'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);disp(out);
end
