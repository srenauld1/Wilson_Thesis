function pizza_ves041_render_eight_columns(X,rel,titleText,subtitleText,stem,isNormalized,flies,eventCounts)
% X contains original event traces or fly-mean traces, one cell per column.
if nargin<7,flies=strings(0,1);end
if nargin<8,eventCounts=cellfun(@(x)size(x,2),X);end
pooled=~isempty(flies);assert(numel(X)==8);
fig=figure('Visible','off','Color','w','Position',[10 10 3300 1400]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,8,'TileSpacing','compact','Padding','loose');
labels={{'VES041 peaks'},{'All slowdowns'},{'Slowdowns without turns'},{'Slowdowns with turns'}, ...
 {'Rotation peaks','Forward drop <1 mm/s'},{'Rotation peaks','Forward drop >=1 mm/s'}, ...
 {'Speed-ups','Without rotation (<15 deg)'},{'Speed-ups','With rotation (>=15 deg)'}};
lims=cell(1,3);allZ=[];
for v=1:3
 z=[];for c=1:8,q=X{c}(:,:,v);z=[z;q(:)];end;z=z(isfinite(z)); %#ok<AGROW>
 if isNormalized,allZ=[allZ;z];elseif isempty(z),lims{v}=[0 1];else,lims{v}=[min(z) max(z)];if diff(lims{v})==0,lims{v}=lims{v}+[-.5 .5];else,lims{v}=lims{v}+[-1 1]*.025*diff(lims{v});end,end %#ok<AGROW>
end
if isNormalized,bound=max(1,ceil(2*max(abs(allZ)))/2);for v=1:3,lims{v}=[-bound bound];end,end
flyLabel=replace(flies,'_fly',' F');flyLabel=replace(flyLabel,'_unlabeled_inferred',' *');flyLabel=replace(flyLabel,'_inferred','*');
for row=1:4,for c=1:8
 ax=nexttile(tl);hold(ax,'on');n=size(X{c},2);
 if row==1
  if n>0,im=imagesc(ax,rel,1:n,X{c}(:,:,1)');im.AlphaData=isfinite(im.CData);set(ax,'YDir','normal','Color',[.9 .9 .9]);ylim(ax,[.5 n+.5]);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');ylim(ax,[.5 1.5]);end
  caxis(ax,lims{1});cb=colorbar(ax);if isNormalized,cb.Label.String='Fly mean \DeltaF/F (z-score)';else,cb.Label.String='Original \DeltaF/F';end
  if pooled,yticks(ax,1:n);yticklabels(ax,flyLabel);ax.YAxis.FontSize=8;ylabel(ax,'Fly (* inferred identity)');count=sprintf('%d flies | %d events',nnz(any(isfinite(X{c}(:,:,1)),1)),eventCounts(c));else,ylabel(ax,'Event (chronological)');count=sprintf('n=%d',eventCounts(c));end
  title(ax,[labels{c} {count}]);
 else
  v=row-1;if n>0,plot(ax,rel,X{c}(:,:,v),'Color',[.75 .75 .75],'LineWidth',.5);plot(ax,rel,mean(X{c}(:,:,v),2,'omitnan'),'k','LineWidth',2.2);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
  if isNormalized,names={'\DeltaF/F (z-score)','Forward velocity (z-score)','|Yaw| speed (z-score)'};yline(ax,0,':','Color',[.45 .45 .45]);else,names={'Original \DeltaF/F','Forward velocity (mm/s)','|Yaw| speed (deg/s)'};end
  ylabel(ax,names{v});ylim(ax,lims{v});grid(ax,'on');
 end
 xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
end,end
if isNormalized,cmap=[linspace(.15,1,128)' linspace(.3,1,128)' ones(128,1);ones(128,1) linspace(1,.2,128)' linspace(1,.15,128)'];colormap(fig,cmap);else,colormap(fig,parula(256));end
title(tl,titleText,'FontWeight','bold');subtitle(tl,subtitleText,'FontSize',10);assert(numel(findall(fig,'Type','axes'))==32);
exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
end
