function export_pizza_optomotor_columns_svg()
% Split the exact saved baseline-subtracted presentation into two SVGs.
% Copy axes and their data/limits; do not redetect events or renormalize.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis';
source=fullfile(root,'ves041_noheat_imaging','gratings','optomotor_exploration', ...
 'ves041_optomotor_publication_baseline_subtracted.fig');
out=fullfile(root,'data_crunches');assert(isfolder(out));
old=openfig(source,'invisible');cleanup=onCleanup(@()close(old)); %#ok<NASGU>
axesOld=findall(old,'Type','axes');ordered=gobjects(6,1);
for a=axesOld',ordered(a.Layout.Tile)=a;end
names={'ves041_optomotor_all_visual_motion_baseline_subtracted.svg', ...
 'ves041_optomotor_matched_turns_baseline_subtracted.svg'};
colors=[.15 .43 .70;.86 .34 .17;.48 .31 .65];
for c=1:2
 fig=figure('Visible','off','Color','w','Position',[20 20 730 1080]);
 tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','loose');
 for row=1:3
  previous=ordered((row-1)*2+c);ax=copyobj(previous,fig);ax.Parent=tl;ax.Layout.Tile=row;
  ax.FontSize=18;
  labels={{'Fly yaw','(deg/s)'},{'Forward speed','(mm/s)'},{'Baseline-subtracted','dF/F'}};
  ylabel(ax,labels{row},'FontSize',21);
  assert(isequal(ax.XLim,previous.XLim)&&isequal(ax.YLim,previous.YLim));
  if row==1
   title(ax,previous.Title.String,'FontSize',23);
   hold(ax,'on');h=gobjects(4,1);
   for k=1:3,h(k)=plot(ax,NaN,NaN,'Color',colors(k,:),'LineWidth',1.7);end
   h(4)=plot(ax,NaN,NaN,'k','LineWidth',3);
   lg=legend(ax,h,{'Aug 12, fly 3','Aug 13, fly 1','Aug 19, fly 1','Equal-fly mean +/- SEM'}, ...
    'FontSize',13,'NumColumns',2,'Box','off');lg.Layout.Tile='south';
  end
 end
 subtitle(tl,{'VES041 | 3 flies', 'dF/F baseline: -1 to -0.25 s before visual motion onset'},'FontSize',14);
 print(fig,fullfile(out,names{c}),'-dsvg','-painters');close(fig);
 fprintf('Saved %s\n',fullfile(out,names{c}));
end
end
