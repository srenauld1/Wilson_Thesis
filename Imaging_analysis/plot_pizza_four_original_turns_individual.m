function plot_pizza_four_original_turns_individual()
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/pooled_four_requested_trials_oldcriteria_20261102';S=load(fullfile(root,'crunch_data.mat'));out=fullfile(root,'individual_trials');if ~isfolder(out),mkdir(out);end
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};light={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};labs={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};names={'Control turns','Heat turns'};
for k=1:numel(S.ids)
 fig=figure('Visible','off','Color','w','Position',[20 20 1120 1050]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');
 for row=1:3
  z=[];for s=1:2,q=S.X{k,s}(:,:,row);z=[z;q(:)];end
  if isempty(z),lim=[0 1];else,pad=max(.01,.06*range(z));lim=[min(z)-pad max(z)+pad];end
  for s=1:2
   ax=nexttile(tl);hold(ax,'on');n=S.counts(k,s);
   if n>0
    X=S.X{k,s}(:,:,row);mu=mean(X,2);se=std(X,0,2)/sqrt(n);
    plot(ax,S.rel,X,'Color',light{row},'LineWidth',.65,'Tag','eventTrace');
    if n>1,fill(ax,[S.rel;flipud(S.rel)],[mu-se;flipud(mu+se)],colors{row},'FaceAlpha',.2,'EdgeColor','none');end
    plot(ax,S.rel,mu,'Color',colors{row},'LineWidth',3);
   else,text(ax,.5,.5,'No qualifying turns','Units','normalized','HorizontalAlignment','center','FontSize',17);end
   xline(ax,0,':','Color',[.45 .45 .45]);xlim(ax,[-2 3]);ylim(ax,lim);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
   if s==1,ylabel(ax,labs{row},'FontSize',21);end
   if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',17);else,ax.XTickLabel=[];end
   if row==1,title(ax,names{s},'FontSize',23,'FontWeight','bold');subtitle(ax,sprintf('n=%d events',n),'FontSize',15);end
  end
 end
 subtitle(tl,{strrep(char(S.ids(k)),'_',' '),'Original turn criteria | Thin: events | Thick: mean | Shading: SEM across events'},'FontSize',15);
 assert(numel(findall(fig,'Type','axes'))==6&&numel(findall(fig,'Tag','eventTrace'))==3*sum(S.counts(k,:)));
 stem=fullfile(out,char(S.ids(k)+"_original_criteria_control_heat_turns"));exportgraphics(fig,[stem '.png'],'Resolution',180);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
 fprintf('%s: control %d, heat %d\n',S.ids(k),S.counts(k,:));
end
end
