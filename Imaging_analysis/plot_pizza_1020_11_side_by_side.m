function plot_pizza_1020_11_side_by_side(folder)
% Re-layout saved data: peak control/heat, then slowdown+rotation control/heat.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if nargin<1,folder='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/20261020_11_2/presentation_control_vs_heat_20261102';end
S=load(fullfile(folder,'comparison_data.mat'));rel=S.rel;mapping=[1 1;1 2;2 1;2 2];
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};thinColors={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};groups={'VES041 peaks','Slowdowns + rotation'};conditions={'Control barriers','Heat barriers'};
fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
for row=1:3
 z=[];for c=1:2,for s=1:2,q=S.D{s,c}(:,:,row);z=[z;q(:)];end,end
 q=S.mu(:,row,:,1:2)-S.sem(:,row,:,1:2);r=S.mu(:,row,:,1:2)+S.sem(:,row,:,1:2);z=[z;q(:);r(:)];z=z(isfinite(z));pad=max(.01,.06*range(z));lims=[min(z)-pad max(z)+pad];
 for panel=1:4
  c=mapping(panel,1);s=mapping(panel,2);ax=nexttile(tl);hold(ax,'on');X=S.D{s,c}(:,:,row);
  plot(ax,rel,X,'Color',thinColors{row},'LineWidth',.5,'Tag','eventTrace');m=S.mu(:,row,s,c);se=S.sem(:,row,s,c);
  fill(ax,[rel;flipud(rel)],[m-se;flipud(m+se)],colors{row},'FaceAlpha',.20,'EdgeColor','none','Tag','eventSem');
  plot(ax,rel,m,'Color',colors{row},'LineWidth',3,'Tag','conditionMean');xline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);
  xlim(ax,[-2 3]);ylim(ax,lims);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if panel==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time (s)','FontSize',21);else,ax.XTickLabel=[];end
  if row==1,title(ax,{groups{c},conditions{s}},'FontSize',23,'FontWeight','bold');subtitle(ax,sprintf('n=%d events',S.counts(s,c)),'FontSize',13);end
 end
end
caption='2026-10-20 trial 11-2';if isfield(S,'decelerationCutoffRemoved')&&S.decelerationCutoffRemoved,caption=[caption ' | No minimum deceleration magnitude'];end
subtitle(tl,{caption,'Thin: individual events | Thick: mean | Shading: SEM across events | One recording, one fly'},'FontSize',16);
assert(numel(findall(fig,'Type','axes'))==12);assert(numel(findall(fig,'Tag','eventTrace'))==3*sum(S.counts(:,1:2),'all'));
stem=fullfile(folder,'20261020_11_2_control_heat_side_by_side_first_two_selections');
exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
